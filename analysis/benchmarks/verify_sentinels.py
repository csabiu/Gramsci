"""Re-verify sentinel benchmark timings on a quiescent machine.

Re-runs a representative subset of the paper's benchmark points (3 repeats
each) and compares against the values stored in benchmarks{,2}.json.
Run this on an otherwise idle node before submission; deviations beyond
TOLERANCE are flagged and indicate the stored sweep should be re-run.

Usage:  python3 verify_sentinels.py
"""
import json
import os
import re
import subprocess

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
CPU = os.path.join(ROOT, 'bin', 'gramsci')
GPU = os.path.join(ROOT, 'bin', 'gramsci_gpu')
GAL = os.path.join(ROOT, 'example', 'test.gal')
RAN = os.path.join(ROOT, 'example', 'test.ran')
NTH = os.cpu_count() or 1
TOLERANCE = 0.05      # 5%
REPEATS = 3


def run(exe, rmax, nbins, mode, env_extra=None):
    env = os.environ.copy()
    env.update(env_extra or {})
    cmd = [exe, '-gal', GAL, '-ran', RAN, '-rmin', '1.0',
           '-rmax', str(rmax), '-nbins', str(nbins), '-nmu', '1',
           '-out', '/tmp/sentinel.out', f'-{mode}']
    res = subprocess.run(cmd, capture_output=True, text=True, env=env)
    q = re.search(r'Querying graph took\s+([\d.]+)', res.stdout)
    return float(q.group(1))


def main():
    b1 = json.load(open(os.path.join(HERE, 'benchmarks.json')))
    b2 = json.load(open(os.path.join(HERE, 'benchmarks2.json')))

    def ref1(sweep, rmax, side):
        return [r for r in b1['sweeps'][sweep] if r['rmax'] == rmax][0][side]['query_s']

    cpu_env = {'OMP_NUM_THREADS': str(NTH)}
    bs_env = {'OMP_NUM_THREADS': str(NTH), 'GRAMSCI_BSEARCH': '1'}

    sentinels = [
        ('3PCF rmax=80 CPUx64', ref1('3pcf_rmax', 80, 'cpu'),
         lambda: run(CPU, 80, 6, '3pcf', cpu_env)),
        ('3PCF rmax=80 GPU', ref1('3pcf_rmax', 80, 'gpu'),
         lambda: run(GPU, 80, 6, '3pcf')),
        ('4PCF rmax=60 CPUx64', ref1('4pcf_rmax', 60, 'cpu'),
         lambda: run(CPU, 60, 3, '4pcf', cpu_env)),
        ('4PCF rmax=60 GPU', ref1('4pcf_rmax', 60, 'gpu'),
         lambda: run(GPU, 60, 3, '4pcf')),
        ('v1 bsearch 3PCF rmax=80',
         [r for r in b2['sweeps']['bsearch_rmax'] if r['rmax'] == 80][0]['cpu_bsearch_3pcf'],
         lambda: run(CPU, 80, 6, '3pcf', bs_env)),
        ('1-thread 3PCF rmax=70',
         [r for r in b2['sweeps']['threads']['cpu'] if r['threads'] == 1][0]['cpu_3pcf'],
         lambda: run(CPU, 70, 6, '3pcf', {'OMP_NUM_THREADS': '1'})),
    ]

    print(f'{REPEATS} repeats per sentinel, tolerance {TOLERANCE:.0%}\n')
    all_ok = True
    for label, ref, fn in sentinels:
        vals = [fn() for _ in range(REPEATS)]
        med = float(np.median(vals))
        dev = abs(med - ref) / ref
        status = 'OK  ' if dev <= TOLERANCE else 'FLAG'
        if dev > TOLERANCE:
            all_ok = False
        print(f'[{status}] {label:28s} stored {ref:9.2f}s  now {med:9.2f}s '
              f'({dev:+.1%}, spread {np.std(vals)/med:.1%})')

    print('\nAll sentinels within tolerance — stored sweeps stand.' if all_ok
          else '\nFLAGGED sentinels — re-run the affected sweeps on this node.')


if __name__ == '__main__':
    main()
