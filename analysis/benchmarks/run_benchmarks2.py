"""Paper benchmark runner, phase 2:

  D. Legacy binary-search (v1) vs merge-walk (v2) CPU kernels vs rmax
     (adds the lineage curve to the scaling figure).
  E. OpenMP strong scaling: query time vs thread count, plus GPU reference.
  F. Density scaling: subsampled catalogs at fixed rmax, CPU x N vs GPU.

Writes benchmarks2.json.  Requires the GRAMSCI_BSEARCH env switch in the
CPU binary.  Run on an idle machine after run_benchmarks.py.
"""
import json
import os
import re
import subprocess
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
CPU = os.path.join(ROOT, 'bin', 'gramsci')
GPU = os.path.join(ROOT, 'bin', 'gramsci_gpu')
GAL = os.path.join(ROOT, 'example', 'test.gal')
RAN = os.path.join(ROOT, 'example', 'test.ran')
OUT = os.path.join(HERE, 'benchmarks2.json')

N_THREADS = os.cpu_count() or 1


def run(exe, gal, ran, rmax, nbins, mode, env_extra=None):
    env = os.environ.copy()
    env.update(env_extra or {})
    cmd = [exe, '-gal', gal, '-ran', ran, '-rmin', '1.0',
           '-rmax', str(rmax), '-nbins', str(nbins), '-nmu', '1',
           '-out', '/tmp/bench_paper2.out', f'-{mode}']
    res = subprocess.run(cmd, capture_output=True, text=True, env=env)
    if res.returncode != 0:
        print(f'    FAILED: {" ".join(cmd)}')
        return None
    q = re.search(r'Querying graph took\s+([\d.]+)', res.stdout)
    return float(q.group(1)) if q else None


def subsample(frac, tag):
    """Subsample gal+ran by the same fraction; cache in /tmp."""
    g = f'/tmp/bench_sub_{tag}.gal'
    r = f'/tmp/bench_sub_{tag}.ran'
    if not (os.path.exists(g) and os.path.exists(r)):
        rng = np.random.default_rng(12345)
        for src, dst in [(GAL, g), (RAN, r)]:
            d = np.loadtxt(src)
            idx = rng.choice(len(d), size=int(frac * len(d)), replace=False)
            np.savetxt(dst, d[np.sort(idx)], fmt='%16.8e')
    return g, r


def main():
    results = {'n_threads': N_THREADS, 'sweeps': {}}

    # ---- Sweep D: bsearch (v1) vs rmax, CPU x max threads ------------
    sweep = []
    for rmax in [30, 40, 50, 60, 70, 80, 90, 100]:
        print(f'v1 bsearch 3PCF rmax={rmax}...', flush=True)
        t3 = run(CPU, GAL, RAN, rmax, 6, '3pcf',
                 {'OMP_NUM_THREADS': str(N_THREADS), 'GRAMSCI_BSEARCH': '1'})
        t4 = None
        if rmax <= 70:
            print(f'v1 bsearch 4PCF rmax={rmax}...', flush=True)
            t4 = run(CPU, GAL, RAN, rmax, 3, '4pcf',
                     {'OMP_NUM_THREADS': str(N_THREADS), 'GRAMSCI_BSEARCH': '1'})
        sweep.append({'rmax': rmax, 'cpu_bsearch_3pcf': t3,
                      'cpu_bsearch_4pcf': t4})
        print(f'  -> 3pcf {t3}s  4pcf {t4}s')
    results['sweeps']['bsearch_rmax'] = sweep

    # ---- Sweep E: strong scaling with thread count -------------------
    # 3PCF rmax=70, 4PCF rmax=40 (sized so 1-thread runs stay manageable)
    sweep = []
    threads = [1, 2, 4, 8, 16, 32, 64]
    for nt in [t for t in threads if t <= N_THREADS]:
        print(f'threads={nt}: 3PCF rmax=70...', flush=True)
        t3 = run(CPU, GAL, RAN, 70, 6, '3pcf', {'OMP_NUM_THREADS': str(nt)})
        print(f'threads={nt}: 4PCF rmax=40...', flush=True)
        t4 = run(CPU, GAL, RAN, 40, 3, '4pcf', {'OMP_NUM_THREADS': str(nt)})
        sweep.append({'threads': nt, 'cpu_3pcf': t3, 'cpu_4pcf': t4})
        print(f'  -> 3pcf {t3}s  4pcf {t4}s')
    results['sweeps']['threads'] = {
        'cpu': sweep,
        'gpu_3pcf_rmax70': run(GPU, GAL, RAN, 70, 6, '3pcf'),
        'gpu_4pcf_rmax40': run(GPU, GAL, RAN, 40, 3, '4pcf'),
    }
    print(f"GPU refs: {results['sweeps']['threads']['gpu_3pcf_rmax70']}s / "
          f"{results['sweeps']['threads']['gpu_4pcf_rmax40']}s")

    # ---- Sweep F: density scaling at fixed rmax ----------------------
    # 3PCF rmax=80, 4PCF rmax=50; subsample fractions of the test catalog
    sweep = []
    for frac in [0.125, 0.25, 0.5, 1.0]:
        if frac == 1.0:
            g, r = GAL, RAN
        else:
            g, r = subsample(frac, f'{frac:.3f}'.replace('.', 'p'))
        row = {'frac': frac}
        print(f'frac={frac}: GPU 3PCF rmax=80...', flush=True)
        row['gpu_3pcf'] = run(GPU, g, r, 80, 6, '3pcf')
        print(f'frac={frac}: CPU 3PCF rmax=80...', flush=True)
        row['cpu_3pcf'] = run(CPU, g, r, 80, 6, '3pcf',
                              {'OMP_NUM_THREADS': str(N_THREADS)})
        print(f'frac={frac}: GPU 4PCF rmax=50...', flush=True)
        row['gpu_4pcf'] = run(GPU, g, r, 50, 3, '4pcf')
        print(f'frac={frac}: CPU 4PCF rmax=50...', flush=True)
        row['cpu_4pcf'] = run(CPU, g, r, 50, 3, '4pcf',
                              {'OMP_NUM_THREADS': str(N_THREADS)})
        sweep.append(row)
        print(f'  -> {row}')
    results['sweeps']['density'] = sweep

    with open(OUT, 'w') as f:
        json.dump(results, f, indent=1)
    print(f'wrote {OUT}')


if __name__ == '__main__':
    main()
