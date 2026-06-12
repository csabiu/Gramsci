"""Parity-estimator signal-recovery test (paper figure).

Catalogs of isolated chiral tetrahedra with left-handed fraction x are
generated; each left unit contributes +1 and each right unit -1 to the
parity-odd channel at the tetrahedron's configuration, so the estimator
must recover

    zeta_odd / zeta_even  =  2x - 1

exactly at the injected configuration (edge bins (1,2,2,2,3,3)).
Sweeps x in [0.5, 1.0], runs the GPU parity kernel, and plots measured
vs injected asymmetry.

Output: parity_recovery.{pdf,png}, parity_recovery.json
"""
import json
import os
import re
import subprocess
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'tests'))
from generate_test_fields import generate_chiral_tetra, RMIN, RMAX, NBINS, \
    BOX_SIZE, N_RAND   # noqa: E402

GPU = os.path.join(ROOT, 'bin', 'gramsci_gpu')
N_STRUCT = 500
TARGET_BINS = (1, 2, 2, 2, 3, 3)

plt.rcParams.update({
    'font.family': 'serif', 'mathtext.fontset': 'stix', 'font.size': 9,
    'axes.labelsize': 10, 'legend.fontsize': 8,
    'xtick.direction': 'in', 'ytick.direction': 'in',
    'xtick.top': True, 'ytick.right': True, 'figure.dpi': 150,
})


def measure(frac, seed):
    gal, ran, out = '/tmp/pr.gal', '/tmp/pr.ran', '/tmp/pr.4pcfp'
    generate_chiral_tetra(N_STRUCT, BOX_SIZE, N_RAND, seed, gal, ran,
                          handedness=float(frac))
    cmd = [GPU, '-gal', gal, '-ran', ran, '-rmin', str(RMIN),
           '-rmax', str(RMAX), '-nbins', str(NBINS), '-nmu', '1',
           '-out', out, '-4pcfp']
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0, res.stdout[-500:]

    # nvfortran wraps the long header over two lines; skip non-numeric rows
    n_skip = 0
    with open(out) as f:
        for line in f:
            s = line.strip()
            if s and (s[0].isdigit() or s[0] in '+-'):
                break
            n_skip += 1
    data = np.loadtxt(out, skiprows=n_skip)
    # columns: 12 bin edges | NNNN RRRR zeta_even NNNN_odd RRRR_odd zeta_odd ...
    width = (RMAX - RMIN) / NBINS
    binidx = np.rint((data[:, 0:12:2] - RMIN) / width).astype(int) + 1
    # Match the exact stored canonical tuple: several S4 orbits share the
    # same sorted multiset, so sorting would be ambiguous.
    row = np.all(binidx == TARGET_BINS, axis=1)
    assert row.sum() == 1
    nnnn_even = data[row, 12][0]
    nnnn_odd = data[row, 15][0]
    # Count-normalised odd channel: weights are 1/N_gal each, the kernel
    # counts each tetrahedron once, and mixed data-random terms cancel in
    # the odd channel in expectation, so
    #     NNNN_odd * N_gal^4 = N_left - N_right   exactly.
    # (The even channel carries a constant mixed-term offset, so the
    # zeta-/zeta+ ratio is biased at the ~7% level and is not used.)
    n_gal = 4 * N_STRUCT
    return nnnn_odd * n_gal**4 / N_STRUCT


def main():
    fracs = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    seeds = [101, 202, 303]
    rows = []
    for x in fracs:
        vals = [measure(x, s) for s in seeds]
        rows.append({'frac': x, 'injected': 2 * x - 1,  # = (N_L-N_R)/N_struct
                     'measured': float(np.mean(vals)),
                     'scatter': float(np.std(vals))})
        print(rows[-1])

    with open(os.path.join(HERE, 'parity_recovery.json'), 'w') as f:
        json.dump(rows, f, indent=1)

    inj = np.array([r['injected'] for r in rows])
    mea = np.array([r['measured'] for r in rows])
    sca = np.array([r['scatter'] for r in rows])

    fig, ax = plt.subplots(figsize=(3.6, 3.3))
    ax.plot([0, 1], [0, 1], '-', color='gray', lw=0.8, label='expected (1:1)')
    ax.errorbar(inj, mea, yerr=sca, fmt='o', color='C3', ms=5, capsize=3,
                label=f'measured ({N_STRUCT} tetrahedra, 3 seeds)')
    ax.set_xlabel(r'injected asymmetry  $(N_L - N_R)/N_{\rm struct}$')
    ax.set_ylabel(r'recovered  $\mathrm{NNNN}^{-} N_{\rm gal}^{4}/N_{\rm struct}$')
    ax.legend(loc='upper left', frameon=False)
    plt.tight_layout()
    fig.savefig(os.path.join(HERE, 'parity_recovery.pdf'))
    fig.savefig(os.path.join(HERE, 'parity_recovery.png'))
    print('wrote parity_recovery.[pdf,png]')
    print('max |measured - injected|:', np.max(np.abs(mea - inj)))


if __name__ == '__main__':
    main()
