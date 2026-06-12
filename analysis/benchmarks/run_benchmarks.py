"""Paper benchmark runner: scaling with rmax and out-of-core overhead.

Produces benchmarks.json with all raw timings.  Run on an otherwise idle
machine (GPU and all CPU cores are used).

Sweeps (test catalog: 207k galaxies + 500k randoms):
  A. 3PCF query time vs rmax:  CPU x 64 threads and GPU
  B. 4PCF query time vs rmax:  CPU x 64 threads (to rmax=70) and GPU (to 80)
  C. Out-of-core overhead: 3PCF rmax=80 with 1..8 forced edge windows
"""
import json
import os
import re
import subprocess
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
CPU = os.path.join(ROOT, 'bin', 'gramsci')
GPU = os.path.join(ROOT, 'bin', 'gramsci_gpu')
GAL = os.path.join(ROOT, 'example', 'test.gal')
RAN = os.path.join(ROOT, 'example', 'test.ran')
OUT = os.path.join(HERE, 'benchmarks.json')

N_THREADS = os.cpu_count() or 1


def run(exe, rmax, nbins, mode, env_extra=None, rmin=1.0):
    env = os.environ.copy()
    env.update(env_extra or {})
    cmd = [exe, '-gal', GAL, '-ran', RAN, '-rmin', str(rmin),
           '-rmax', str(rmax), '-nbins', str(nbins), '-nmu', '1',
           '-out', '/tmp/bench_paper.out', f'-{mode}']
    t0 = time.perf_counter()
    res = subprocess.run(cmd, capture_output=True, text=True, env=env)
    wall = time.perf_counter() - t0
    if res.returncode != 0:
        print(f'    FAILED: {" ".join(cmd)}')
        return None
    q = re.search(r'Querying graph took\s+([\d.]+)', res.stdout)
    g = re.search(r'Creating the graph took\s+([\d.]+)', res.stdout)
    e = re.search(r'CSR graph: \d+ nodes, (\d+) edges', res.stdout)
    return {
        'query_s': float(q.group(1)) if q else None,
        'graph_s': float(g.group(1)) if g else None,
        'wall_s': wall,
        'edges': int(e.group(1)) if e else None,
    }


def main():
    results = {'n_threads': N_THREADS, 'sweeps': {}}

    # ---- Sweep A: 3PCF vs rmax --------------------------------------
    sweep = []
    for rmax in [30, 40, 50, 60, 70, 80, 90, 100]:
        row = {'rmax': rmax}
        print(f'3PCF rmax={rmax}: GPU...', flush=True)
        row['gpu'] = run(GPU, rmax, 6, '3pcf')
        print(f'3PCF rmax={rmax}: CPU x{N_THREADS}...', flush=True)
        row['cpu'] = run(CPU, rmax, 6, '3pcf',
                         {'OMP_NUM_THREADS': str(N_THREADS)})
        sweep.append(row)
        print(f'  -> cpu {row["cpu"]["query_s"]}s  gpu {row["gpu"]["query_s"]}s')
    results['sweeps']['3pcf_rmax'] = sweep

    # ---- Sweep B: 4PCF vs rmax --------------------------------------
    sweep = []
    for rmax in [30, 40, 50, 60, 65, 70, 80]:
        row = {'rmax': rmax}
        print(f'4PCF rmax={rmax}: GPU...', flush=True)
        row['gpu'] = run(GPU, rmax, 3, '4pcf')
        if rmax <= 70:
            print(f'4PCF rmax={rmax}: CPU x{N_THREADS}...', flush=True)
            row['cpu'] = run(CPU, rmax, 3, '4pcf',
                             {'OMP_NUM_THREADS': str(N_THREADS)})
        else:
            row['cpu'] = None
        sweep.append(row)
        c = row['cpu']['query_s'] if row['cpu'] else '--'
        print(f'  -> cpu {c}s  gpu {row["gpu"]["query_s"]}s')
    results['sweeps']['4pcf_rmax'] = sweep

    # ---- Sweep C: out-of-core overhead ------------------------------
    # Fixed problem: 3PCF rmax=80.  First learn the edge count, then force
    # nwin windows via GRAMSCI_GPU_WIN_EDGES.
    probe = run(GPU, 80, 6, '3pcf')
    edges = probe['edges']
    print(f'chunk sweep: {edges} edges at rmax=80')
    sweep = [{'nwin': 1, 'res': probe}]
    for nwin in [2, 3, 4, 6, 8]:
        we = edges // nwin + 1
        print(f'3PCF rmax=80, {nwin} windows (win_edges={we})...', flush=True)
        r = run(GPU, 80, 6, '3pcf', {'GRAMSCI_GPU_WIN_EDGES': str(we)})
        sweep.append({'nwin': nwin, 'res': r})
        print(f'  -> {r["query_s"]}s')
    results['sweeps']['3pcf_chunks'] = sweep

    with open(OUT, 'w') as f:
        json.dump(results, f, indent=1)
    print(f'wrote {OUT}')


if __name__ == '__main__':
    main()
