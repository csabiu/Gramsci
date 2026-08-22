"""Exact-antisymmetry test of the parity-odd estimator.

Mirroring a catalogue (x -> L-x) reverses the handedness of every tetrahedron,
so a CORRECT parity-odd estimator must satisfy, configuration by configuration
and to floating-point precision,

    NNNN_odd(mirrored data, mirrored randoms) = - NNNN_odd(data, randoms)
    NNNN     (mirrored)                       = + NNNN     (original)

This is a property of the ESTIMATOR alone: it holds for any input catalogue,
clustered or not, and needs no ensemble and no reference implementation.  It is
the sharpest cheap statement that the odd channel is computed correctly, and it
catches sign errors, canonicalisation mistakes and mode-dependent bugs that a
null test on a parity-symmetric field cannot distinguish from noise.

The pixelized modes can satisfy it only if the direction grid is itself
mirror-symmetric under the chosen reflection (true for x -> L-x when n_phi is
even); -exactparity satisfies it by construction, up to last-bit sign flips on
tetrahedra that are coplanar to machine precision.

Usage: python tests/test_parity_mirror.py [--nreal N] [--tol T]
"""
import argparse
import os
import subprocess
import sys

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.join(ROOT, 'tests'))
GPU = os.path.join(ROOT, 'bin', 'gramsci_gpu')
CPU = os.path.join(ROOT, 'bin', 'gramsci')
WORK = os.environ.get('WORK', '/tmp/parity_mirror')
os.makedirs(WORK, exist_ok=True)

# Sized so the mean neighbour count stays near 20: the quadruplet
# enumeration grows as its cube, and this check is exact rather than
# statistical, so it needs only enough tetrahedra to be meaningful.
BOX, RMIN, RMAX, NBINS = 1000.0, 10.0, 45.0, 4
N_GAL, N_RAN = 20_000, 40_000
MODES = [('4x16', ['-ntheta', '4', '-nphi', '16']),
         ('8x32', ['-ntheta', '8', '-nphi', '32']),
         ('exact', ['-exactparity'])]
C_NNNN, C_RRRR, C_ODD = 12, 13, 15


def clustered_catalogue(seed):
    """A clumpy but statistically parity-symmetric point set: isotropic
    Gaussian clumps at uniformly random centres, plus a uniform component."""
    rng = np.random.default_rng(seed)
    n_clump = N_GAL // 60
    centres = rng.uniform(0, BOX, size=(n_clump, 3))
    per = np.full(n_clump, 40)
    pts = np.repeat(centres, per, axis=0)
    pts = pts + rng.normal(0.0, 8.0, size=pts.shape)
    pts = np.vstack([pts, rng.uniform(0, BOX, size=(N_GAL - len(pts), 3))])
    pts = np.mod(pts, BOX)
    ran = rng.uniform(0, BOX, size=(N_RAN, 3))
    return pts, ran


def write(path, pos):
    np.savetxt(path, np.column_stack([pos, np.ones(len(pos))]), fmt='%.8e')


def run(binary, gal, ran, out, extra):
    cmd = [binary, '-gal', gal, '-ran', ran, '-rmin', str(RMIN), '-rmax',
           str(RMAX), '-nbins', str(NBINS), '-nmu', '1', '-out', out,
           '-4pcfp'] + extra
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0, ' '.join(cmd) + '\n' + res.stdout[-800:]
    return np.loadtxt(out, comments=('#', 'r'), ndmin=2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--nreal', type=int, default=2)
    ap.add_argument('--tol', type=float, default=1e-5,
                    help='exact mode: max residual as a fraction of the total '
                         'even-channel weight (pixel modes must be exact)')
    ap.add_argument('--cpu', action='store_true', help='test the CPU binary')
    args = ap.parse_args()
    binary = CPU if args.cpu else GPU

    print(f'{"seed":>6s} {"mode":>7s} {"max|odd_m+odd_o|":>18s} '
          f'{"total weight":>12s} {"fraction":>10s} {"even match":>11s}')
    worst = {tag: 0.0 for tag, _ in MODES}
    even_ok_all = True
    for seed in range(7000, 7000 + args.nreal):
        pts, ran = clustered_catalogue(seed)
        mir = pts.copy(); mir[:, 0] = BOX - mir[:, 0]
        mran = ran.copy(); mran[:, 0] = BOX - mran[:, 0]
        g0, r0 = f'{WORK}/o.gal', f'{WORK}/o.ran'
        g1, r1 = f'{WORK}/m.gal', f'{WORK}/m.ran'
        write(g0, pts); write(r0, ran); write(g1, mir); write(r1, mran)

        for tag, extra in MODES:
            a = run(binary, g0, r0, f'{WORK}/o_{tag}.4pcfp', extra)
            b = run(binary, g1, r1, f'{WORK}/m_{tag}.4pcfp', extra)
            ok = a[:, C_RRRR] > 0
            resid = np.abs(a[ok, C_ODD] + b[ok, C_ODD]).max()
            # Scale by the TOTAL enumerated weight: a single misclassified
            # tetrahedron contributes at most its own weight, so this measures
            # "what fraction of the tetrahedra changed sign".
            total = a[ok, C_NNNN].sum()
            frac = resid / total if total > 0 else 0.0
            even_ok = np.array_equal(a[:, C_NNNN], b[:, C_NNNN])
            even_ok_all &= even_ok
            worst[tag] = max(worst[tag], frac)
            print(f'{seed:6d} {tag:>7s} {resid:18.3e} {total:12.3e} '
                  f'{frac:10.2e} {str(even_ok):>11s}', flush=True)
        for f in (g0, r0, g1, r1):
            os.remove(f)

    print('\n=== verdict ===')
    bad = False
    for tag, _ in MODES:
        # Pixel centres map exactly onto pixel centres under x -> L-x for even
        # n_phi, so the pixelized modes must be antisymmetric to the last bit.
        # -exactparity recomputes the spokes from mirrored coordinates, where
        # (L-x1)-(L-x2) is not bit-identical to -(x1-x2), so tetrahedra that
        # are coplanar to machine precision may be classified either way.
        limit = 0.0 if tag != 'exact' else args.tol
        good = worst[tag] <= limit
        bad |= not good
        print(f'  {tag:>7s}: worst fraction {worst[tag]:.2e} '
              f'(limit {limit:.0e})  {"PASS" if good else "FAIL"}')
    print(f'  even channel identical under mirroring: {even_ok_all}')
    if bad or not even_ok_all:
        sys.exit(1)


if __name__ == '__main__':
    main()
