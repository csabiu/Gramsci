"""Pixel-resolution convergence of the parity kernel.

The parity sign is (by default) computed from spoke directions quantized onto
an n_theta x n_phi grid on the unit sphere.  Coarse grids attenuate the
recovered odd amplitude for near-degenerate tetrahedra (small unit-direction
spoke volume): the recovered asymmetry is eta * (2x-1) with eta < 1.  This
test measures eta as a function of grid resolution and checks it against the
purely geometric prediction, and against the exact-position mode
(-exactparity), which has eta == 1 by construction.

Requires the runtime grid flags -ntheta / -nphi.

Usage:
    python tests/test_pixel_convergence.py [--quick] [--shapes a,b] [--out FILE]
"""
import argparse
import json
import os
import subprocess
import sys
import time

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.join(ROOT, 'tests'))
GPU = os.path.join(ROOT, 'bin', 'gramsci_gpu')
WORK = os.environ.get('WORK', '/tmp/pixel_conv_test')
os.makedirs(WORK, exist_ok=True)

from generate_test_fields import generate_chiral_tetra_shapes
from parity_shapes import (SHAPES, canonicalise, edge_lengths,
                            choose_binning)

N_STRUCT, N_RAND = 500, 20000
SEEDS = (1, 2, 3)
SCALE = 15.0
# (n_theta, n_phi); None == exact-position mode
GRIDS = [(4, 16), (6, 24), (8, 32), (11, 45), (16, 64), (23, 91), None]
QUICK_GRIDS = [(4, 16), (8, 32), None]


# ---------------------------------------------------------------------------
# Geometric prediction: mean pixel-snapped volume sign per left-handed
# structure over random orientations, mirroring the kernel's rule exactly.
# No reference to the code -- pure geometry.
# ---------------------------------------------------------------------------
def eta_pred(verts, nt, nphi, n=400_000, seed=7):
    rng = np.random.default_rng(seed)
    q = rng.standard_normal((n, 4))
    q /= np.linalg.norm(q, axis=1)[:, None]
    w, x, y, z = q.T
    R = np.empty((n, 3, 3))
    R[:, 0, 0] = 1 - 2 * (y * y + z * z); R[:, 0, 1] = 2 * (x * y - z * w)
    R[:, 0, 2] = 2 * (x * z + y * w);     R[:, 1, 0] = 2 * (x * y + z * w)
    R[:, 1, 1] = 1 - 2 * (x * x + z * z); R[:, 1, 2] = 2 * (y * z - x * w)
    R[:, 2, 0] = 2 * (x * z - y * w);     R[:, 2, 1] = 2 * (y * z + x * w)
    R[:, 2, 2] = 1 - 2 * (x * x + y * y)
    v = verts - verts.mean(0)
    sp = np.einsum('nij,kj->nki', R, v[1:] - v[0])
    d = sp / np.linalg.norm(sp, axis=2)[:, :, None]

    def snap(u):
        th = np.arccos(np.clip(u[:, 2], -1, 1))
        ph = np.arctan2(u[:, 1], u[:, 0])
        it = np.minimum((th / (np.pi / nt)).astype(int), nt - 1)
        ip = np.minimum(((ph + np.pi) / (2 * np.pi / nphi)).astype(int),
                        nphi - 1)
        tc = (it + 0.5) * np.pi / nt
        pc = (ip + 0.5) * 2 * np.pi / nphi - np.pi
        return np.column_stack([np.sin(tc) * np.cos(pc),
                                np.sin(tc) * np.sin(pc), np.cos(tc)])

    p = np.stack([snap(d[:, k, :]) for k in range(3)], axis=1)
    volp = np.einsum('ni,ni->n', np.cross(p[:, 0], p[:, 1]), p[:, 2])
    volt = np.einsum('ni,ni->n', np.cross(d[:, 0], d[:, 1]), d[:, 2])
    s = np.where(np.abs(volp) < 1e-9, 0.0, np.sign(volp))
    return float((s * np.sign(volt)).mean()), float(np.mean(np.abs(volp) < 1e-9))


# ---------------------------------------------------------------------------
def load(path):
    n_skip = 0
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if s and (s[0].isdigit() or s[0] in '+-'):
                break
            n_skip += 1
    return np.loadtxt(path, skiprows=n_skip)


def recovery(shape, scale, x, seed, grid):
    """eps-corrected count-level recovered asymmetry; grid None -> exact."""
    v = SHAPES[shape] * scale
    edges = edge_lengths(v)
    rmin, rmax, nb = choose_binning(edges)
    spacing = edges.max() + rmax + 10.0
    box = max(500.0, 8.0 * spacing + 1.0)
    gal = os.path.join(WORK, 'pc.gal')
    ran = os.path.join(WORK, 'pc.ran')
    out = os.path.join(WORK, 'pc.4pcfp')
    generate_chiral_tetra_shapes(N_STRUCT, box, N_RAND, seed, gal, ran, v,
                                 handedness=float(x), spacing=spacing)
    cmd = [GPU, '-gal', gal, '-ran', ran, '-rmin', str(rmin), '-rmax',
           str(rmax), '-nbins', str(nb), '-nmu', '1', '-out', out, '-4pcfp']
    if grid is None:
        cmd.append('-exactparity')
    else:
        cmd += ['-ntheta', str(grid[0]), '-nphi', str(grid[1])]
    t0 = time.time()
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0, ' '.join(cmd) + '\n' + res.stdout[-800:]
    dt = time.time() - t0

    d = load(out)
    w = (rmax - rmin) / nb
    enc = [int((e - rmin) // w) + 1 for e in edges]
    canon, eps = canonicalise(enc)
    row_bins = np.rint((d[:, 0:12:2] - rmin) / w).astype(int) + 1
    match = [i for i in range(len(d))
             if canonicalise(list(row_bins[i]))[0] == canon]
    assert len(match) == 1, f'config match failed ({len(match)} rows)'
    n_gal = 4 * N_STRUCT
    return eps * d[match[0], 15] * n_gal ** 4 / N_STRUCT, dt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--shapes', default='regular,scalene,elongated,flattened,sheared')
    ap.add_argument('--out', default=os.path.join(WORK, 'pixel_convergence.json'))
    args = ap.parse_args()

    grids = QUICK_GRIDS if args.quick else GRIDS
    shapes = args.shapes.split(',')
    seeds = SEEDS[:1] if args.quick else SEEDS

    print(f'{"shape":10s} {"grid":>9s} {"npix":>5s} '
          f'{"measured":>16s} {"eta_pred":>9s} {"zerofrac":>8s} {"t[s]":>6s}')
    results = []
    for shape in shapes:
        v = SHAPES[shape] * SCALE
        for grid in grids:
            if grid is None:
                ep, zf, npix, tag = 1.0, 0.0, 0, 'exact'
            else:
                ep, zf = eta_pred(v, grid[0], grid[1])
                npix = grid[0] * grid[1]
                tag = f'{grid[0]}x{grid[1]}'
            vals, ts = [], []
            for s in seeds:
                r, dt = recovery(shape, SCALE, 1.0, s, grid)
                vals.append(r); ts.append(dt)
            m, sd = float(np.mean(vals)), float(np.std(vals))
            print(f'{shape:10s} {tag:>9s} {npix:5d} '
                  f'{m:+8.4f} +- {sd:.4f} {ep:+9.3f} {zf:8.3f} '
                  f'{np.mean(ts):6.1f}', flush=True)
            results.append(dict(shape=shape, grid=tag, npix=npix,
                                measured=m, measured_err=sd,
                                eta_pred=ep, zero_frac=zf,
                                t_mean=float(np.mean(ts))))
    with open(args.out, 'w') as fh:
        json.dump(results, fh, indent=1)
    print(f'\nwrote {args.out}')

    # --- summary: does measured track the prediction, and converge to 1? ---
    print('\n=== convergence summary (x=1) ===')
    for shape in shapes:
        rs = [r for r in results if r['shape'] == shape]
        coarse = next(r for r in rs if r['grid'] == '4x16')
        exact = next((r for r in rs if r['grid'] == 'exact'), None)
        fine = [r for r in rs if r['npix'] >= 256]
        worst = max((abs(r['measured'] - r['eta_pred']) for r in rs
                     if r['npix'] > 0), default=float('nan'))
        line = (f'  {shape:10s} 4x16 {coarse["measured"]:+.3f} '
                f'(pred {coarse["eta_pred"]:+.3f})')
        if fine:
            line += f'   >=256px {min(r["measured"] for r in fine):+.3f}..' \
                    f'{max(r["measured"] for r in fine):+.3f}'
        if exact:
            line += f'   exact {exact["measured"]:+.3f}'
        line += f'   max|meas-pred| {worst:.3f}'
        print(line)


if __name__ == '__main__':
    main()
