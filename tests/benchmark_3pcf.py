"""
Benchmark and correctness comparison: gramsci (CPU) vs gramsci_gpu (GPU).

Usage:
    python benchmark_3pcf.py [--mode 3pcf|4pcf] [--rmin 1.0] [--rmax 50.0] [--nbins N]
"""
import argparse
import os
import re
import subprocess
import sys
import time

import numpy as np

BINDIR  = os.path.join('..', 'bin')
GALFILE = os.path.join('..', 'example', 'test.gal')
RANFILE = os.path.join('..', 'example', 'test.ran')

RTOL = 1e-6
ATOL = 1e-30


# ---------------------------------------------------------------------------
# Runner helpers
# ---------------------------------------------------------------------------

def run(exe, args, outfile, mode_flag, env_extra=None):
    cmd = f"{exe} {args} -out {outfile} {mode_flag}"
    env = os.environ.copy()
    if env_extra:
        env.update(env_extra)
    label = " ".join(f"{k}={v}" for k, v in (env_extra or {}).items())
    print(f"\n>>> {cmd}" + (f"  [{label}]" if label else ""))
    t0 = time.perf_counter()
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True, env=env)
    wall = time.perf_counter() - t0
    if result.returncode != 0:
        print("STDOUT:", result.stdout[-3000:])
        print("STDERR:", result.stderr[-2000:])
        raise RuntimeError(f"{exe} exited with code {result.returncode}")
    return result.stdout, wall


def parse_query_time(stdout):
    m = re.search(r'Querying graph took\s+([\d.]+)\s+seconds', stdout)
    return float(m.group(1)) if m else None


def parse_graph_time(stdout):
    m = re.search(r'Creating the graph took\s+([\d.]+)\s+seconds', stdout)
    return float(m.group(1)) if m else None


# ---------------------------------------------------------------------------
# Output loaders
# ---------------------------------------------------------------------------

def load_output(path):
    # Count header lines: skip lines that don't start with a number.
    # nvfortran may wrap long headers across two lines; gfortran doesn't.
    n_skip = 0
    header_lines = []
    with open(path) as f:
        for line in f:
            stripped = line.strip()
            if stripped and (stripped[0].isdigit() or stripped[0] in '+-'):
                break
            header_lines.append(line)
            n_skip += 1
    header = ''.join(header_lines)
    data = np.genfromtxt(path, skip_header=n_skip)
    return header, data


# ---------------------------------------------------------------------------
# Correctness comparison
# ---------------------------------------------------------------------------

def compare(cpu_data, gpu_data, count_cols, label_cols):
    """
    count_cols: list of (col_index, name) for the raw counts to compare.
    label_cols: number of leading bin-edge columns whose equality is required.
    """
    if cpu_data.shape != gpu_data.shape:
        raise AssertionError(
            f"Shape mismatch: CPU {cpu_data.shape} vs GPU {gpu_data.shape}"
        )

    # Bin edges must be identical
    edge_diff = np.max(np.abs(cpu_data[:, :label_cols] - gpu_data[:, :label_cols]))
    assert edge_diff < 1e-10, f"Bin edges differ by {edge_diff}"

    results = {}
    for col_idx, name in count_cols:
        c = cpu_data[:, col_idx]
        g = gpu_data[:, col_idx]
        cpu_nan = np.isnan(c)
        gpu_nan = np.isnan(g)
        n_nan_mismatch = int(np.sum(cpu_nan ^ gpu_nan))
        mask = ~cpu_nan & ~gpu_nan
        if np.any(mask):
            denom = np.maximum(np.abs(c[mask]), ATOL)
            max_rel = float(np.max(np.abs(c[mask] - g[mask]) / denom))
        else:
            max_rel = 0.0
        results[name] = {'max_rel_diff': max_rel, 'n_nan_disagree': n_nan_mismatch}

    return results


# ---------------------------------------------------------------------------
# Formatting
# ---------------------------------------------------------------------------

def print_banner(title):
    w = 62
    print("\n" + "=" * w)
    print(f"  {title}")
    print("=" * w)


def print_timing(cpuN, gpu, n_cpu, query_label):
    cpuN_wall, cpuN_graph, cpuN_query = cpuN
    gpu_wall,  gpu_graph,  gpu_query  = gpu

    print(f"  NOTE: Graph build is CPU-only in all cases.")
    print(f"        gramsci uses gfortran -Ofast, gramsci_gpu uses nvfortran -O3 -mp.")
    print()
    hdr = f"{'':30s}  {'CPU×'+str(n_cpu):>8s}  {'GPU':>8s}  {'GPU/CPU×N':>9s}"
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))

    def row(label, vN, vG):
        if vN and vG:
            print(f"  {label:28s}  {vN:8.2f}  {vG:8.2f}  {vN/vG:8.2f}x")
        else:
            print(f"  {label:28s}  {'N/A':>8s}  {'N/A':>8s}")

    row("Wall-clock (s)", cpuN_wall, gpu_wall)
    row("Graph build (s)", cpuN_graph, gpu_graph)
    row(f"{query_label} (s)", cpuN_query, gpu_query)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode',  choices=['3pcf', '4pcf'], default='3pcf')
    parser.add_argument('--rmin',  type=float, default=1.0)
    parser.add_argument('--rmax',  type=float, default=50.0)
    parser.add_argument('--nbins', type=int,   default=None,
                        help='number of bins (default: 6 for 3pcf, 3 for 4pcf)')
    args = parser.parse_args()

    if args.nbins is None:
        args.nbins = 6 if args.mode == '3pcf' else 3

    mode_flag   = f'-{args.mode}'
    query_label = f'{args.mode.upper()} query'

    base_args = (
        f"-gal {GALFILE} -ran {RANFILE} "
        f"-rmin {args.rmin} -rmax {args.rmax} -nbins {args.nbins} -nmu 1"
    )

    cpu_exe  = os.path.join(BINDIR, 'gramsci')
    gpu_exe  = os.path.join(BINDIR, 'gramsci_gpu')
    cpuN_out = f'bench_cpuN_{args.mode}.out'
    gpu_out  = f'bench_gpu_{args.mode}.out'
    n_cpu    = os.cpu_count() or 1

    # ---- Run CPU ×N ----
    print_banner(f"CPU run: gramsci ({n_cpu} threads)  [{args.mode}, rmax={args.rmax}, nbins={args.nbins}]")
    cpuN_stdout, cpuN_wall = run(cpu_exe, base_args, cpuN_out, mode_flag,
                                 env_extra={"OMP_NUM_THREADS": str(n_cpu)})
    cpuN_query = parse_query_time(cpuN_stdout)
    cpuN_graph = parse_graph_time(cpuN_stdout)

    # ---- Run GPU ----
    print_banner("GPU run: gramsci_gpu")
    gpu_stdout, gpu_wall = run(gpu_exe, base_args, gpu_out, mode_flag)
    gpu_query = parse_query_time(gpu_stdout)
    gpu_graph = parse_graph_time(gpu_stdout)

    # ---- Timing ----
    print_banner("Timing summary")
    print_timing(
        (cpuN_wall, cpuN_graph, cpuN_query),
        (gpu_wall,  gpu_graph,  gpu_query),
        n_cpu, query_label
    )

    # ---- Correctness ----
    print_banner(f"Correctness check (CPU×{n_cpu} vs GPU)")
    _, cpu_data = load_output(cpuN_out)
    _, gpu_data = load_output(gpu_out)

    nrows = cpu_data.shape[0]
    print(f"  Output rows: {nrows}")

    if args.mode == '3pcf':
        # cols: 6 bin edges | NNN(6) RRR(7) zeta(8)
        count_cols = [(6, 'NNN'), (7, 'RRR'), (8, 'zeta')]
        label_cols = 6
        count_name = 'NNN'
        count_idx  = 6
    else:
        # cols: 12 bin edges | NNNN(12) RRRR(13) zeta(14) zeta_disc(15) zeta_conn(16)
        count_cols = [(12, 'NNNN'), (13, 'RRRR'), (14, 'zeta')]
        label_cols = 12
        count_name = 'NNNN'
        count_idx  = 12

    diffs = compare(cpu_data, gpu_data, count_cols, label_cols)
    all_pass = True
    for name, d in diffs.items():
        ok = d['max_rel_diff'] <= RTOL and d['n_nan_disagree'] == 0
        if not ok:
            all_pass = False
        status = "PASS" if ok else "FAIL"
        print(f"  {name:8s}  max_rel_diff={d['max_rel_diff']:.2e}  "
              f"nan_mismatch={d['n_nan_disagree']}  [{status}]")

    # Print worst row
    c_col = cpu_data[:, count_idx]
    g_col = gpu_data[:, count_idx]
    mask = ~np.isnan(c_col) & ~np.isnan(g_col) & (c_col != 0)
    if np.any(mask):
        rel = np.abs(c_col[mask] - g_col[mask]) / np.maximum(np.abs(c_col[mask]), ATOL)
        worst = np.argmax(rel)
        row_idx = np.where(mask)[0][worst]
        print(f"\n  Worst {count_name} row (index {row_idx}):")
        print(f"    CPU: {cpu_data[row_idx]}")
        print(f"    GPU: {gpu_data[row_idx]}")

    print()
    if all_pass:
        print("  RESULT: CPU and GPU results agree within tolerance.")
    else:
        print("  RESULT: MISMATCH detected — see above.")
        for f in [cpuN_out, gpu_out]:
            if os.path.exists(f): os.remove(f)
        sys.exit(1)

    for f in [cpuN_out, gpu_out]:
        os.remove(f)


if __name__ == '__main__':
    main()
