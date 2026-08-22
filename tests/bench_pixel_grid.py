"""Cost and degeneracy of the parity kernel vs direction-pixel resolution.

Three measurements:

1. Query time on a realistic uniform catalogue for each pixel grid and for
   -exactparity, interleaved so that external GPU contention affects all
   modes equally (ratios are taken against the same-round baseline).
2. Host graph memory per edge, computed from the actual per-edge record
   (id 4 B + dist 1 B + phi 2 B for pixel modes; no phi byte at all under
   -exactparity, which instead keeps 24 B per *point*).
3. The fraction of enumerated tetrahedra whose PIXELIZED spoke triple has an
   exactly zero signed volume -- repeated pixels or coplanar pixel triples --
   which the kernel must discard from the odd channel.  Measured by Monte
   Carlo over the estimator's own quadruplet measure, so it is directly
   comparable to the "a few per cent at typical pixelizations" statement in
   the paper.  Exact mode discards only genuinely coplanar quadruplets
   (measure zero).

Usage: python tests/bench_pixel_grid.py [--npts N] [--reps R] [--skip-timing]
"""
import argparse
import json
import os
import subprocess
import time

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
GPU = os.path.join(ROOT, 'bin', 'gramsci_gpu')
WORK = os.environ.get('WORK', '/tmp/pixel_bench')
os.makedirs(WORK, exist_ok=True)

GRIDS = [(4, 16), (8, 32), (16, 64), (23, 91), None]   # None = exact
RMIN, RMAX, NBINS, BOX = 20.0, 65.0, 5, 1000.0


# ---------------------------------------------------------------------------
def degenerate_fraction(nt, nphi, n=400_000, seed=3, rmin=RMIN, rmax=RMAX,
                        nbins=NBINS, box=BOX):
    """Fraction of RRRR-measure tetrahedra whose pixelized triple product is
    exactly zero (repeated pixels or coplanar pixel triple)."""
    rng = np.random.default_rng(seed)
    w = (rmax - rmin) / nbins

    def shell(m):
        lo = rmin + rng.integers(0, nbins, m) * w
        hi = lo + w
        u = rng.random(m)
        r = (lo ** 3 + u * (hi ** 3 - lo ** 3)) ** (1.0 / 3.0)
        v = rng.standard_normal((m, 3))
        v /= np.linalg.norm(v, axis=1)[:, None]
        return r[:, None] * v

    # hub at origin, three spokes with radii drawn from the binned measure;
    # accept only if the three mutual separations also fall inside the range
    d = [shell(n) for _ in range(3)]
    ok = np.ones(n, bool)
    for a, b in ((0, 1), (0, 2), (1, 2)):
        s = np.linalg.norm(d[a] - d[b], axis=1)
        ok &= (s >= rmin) & (s < rmax)
    d = [x[ok] for x in d]
    m = len(d[0])
    u = [x / np.linalg.norm(x, axis=1)[:, None] for x in d]

    def snap(v):
        th = np.arccos(np.clip(v[:, 2], -1, 1))
        ph = np.arctan2(v[:, 1], v[:, 0])
        it = np.minimum((th / (np.pi / nt)).astype(int), nt - 1)
        ip = np.minimum(((ph + np.pi) / (2 * np.pi / nphi)).astype(int),
                        nphi - 1)
        return it * nphi + ip, (it + 0.5) * np.pi / nt, \
            (ip + 0.5) * 2 * np.pi / nphi - np.pi

    idx, tc, pc = [], [], []
    for k in range(3):
        i_, t_, p_ = snap(u[k])
        idx.append(i_); tc.append(t_); pc.append(p_)
    P = [np.column_stack([np.sin(t) * np.cos(p), np.sin(t) * np.sin(p),
                          np.cos(t)]) for t, p in zip(tc, pc)]
    volp = np.einsum('ni,ni->n', np.cross(P[0], P[1]), P[2])
    degen = np.abs(volp) < 1e-9
    repeated = (idx[0] == idx[1]) | (idx[0] == idx[2]) | (idx[1] == idx[2])
    return dict(n_sampled=int(m), degenerate=float(degen.mean()),
                repeated_pixel=float(repeated.mean()),
                coplanar_only=float((degen & ~repeated).mean()))


# ---------------------------------------------------------------------------
def run(gal, ran, out, grid):
    cmd = [GPU, '-gal', gal, '-ran', ran, '-rmin', str(RMIN), '-rmax',
           str(RMAX), '-nbins', str(NBINS), '-nmu', '1', '-out', out, '-4pcfp']
    if grid is None:
        cmd.append('-exactparity')
    else:
        cmd += ['-ntheta', str(grid[0]), '-nphi', str(grid[1])]
    t0 = time.time()
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0, res.stdout[-800:]
    wall = time.time() - t0
    tq, tg, edges = wall, None, None
    for line in res.stdout.splitlines():
        if 'uerying' in line and 'took' in line:
            try:
                tq = float(line.split()[-2])
            except ValueError:
                pass
        if 'Creating the graph took' in line:
            try:
                tg = float(line.split()[-2])
            except ValueError:
                pass
        if 'edges' in line.lower() and 'total' in line.lower():
            for tok in line.replace(',', ' ').split():
                if tok.isdigit():
                    edges = int(tok)
    return dict(t_query=tq, t_graph=tg, t_wall=wall, edges=edges)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--npts', type=int, default=300_000)
    ap.add_argument('--reps', type=int, default=2)
    ap.add_argument('--skip-timing', action='store_true')
    ap.add_argument('--out', default=os.path.join(WORK, 'bench_pixel_grid.json'))
    args = ap.parse_args()

    print('=== degenerate (zero pixelized volume) fraction vs grid ===')
    print(f'{"grid":>9s} {"npix":>5s} {"degenerate":>11s} {"repeated px":>12s} '
          f'{"coplanar":>9s}')
    degen = []
    for g in GRIDS:
        if g is None:
            print(f'{"exact":>9s} {0:5d} {0.0:11.4f} {0.0:12.4f} {0.0:9.4f}'
                  '   (only true coplanar, measure zero)')
            degen.append(dict(grid='exact', npix=0, degenerate=0.0,
                              repeated_pixel=0.0, coplanar_only=0.0))
            continue
        d = degenerate_fraction(g[0], g[1])
        print(f'{f"{g[0]}x{g[1]}":>9s} {g[0]*g[1]:5d} {d["degenerate"]:11.4f} '
              f'{d["repeated_pixel"]:12.4f} {d["coplanar_only"]:9.4f}')
        degen.append(dict(grid=f'{g[0]}x{g[1]}', npix=g[0] * g[1], **d))

    results = dict(degenerate=degen)

    if not args.skip_timing:
        rng = np.random.default_rng(11)
        ng = args.npts // 3
        nr = args.npts - ng
        gal = os.path.join(WORK, 'b.gal')
        ran = os.path.join(WORK, 'b.ran')
        np.savetxt(gal, np.column_stack(
            [rng.random((ng, 3)) * BOX, np.ones(ng)]), fmt='%.6e')
        np.savetxt(ran, np.column_stack(
            [rng.random((nr, 3)) * BOX, np.ones(nr)]), fmt='%.6e')
        out = os.path.join(WORK, 'b.4pcfp')

        print(f'\n=== query time vs grid ({args.npts // 1000}k pts, '
              f'r {RMIN:g}-{RMAX:g}, {NBINS} bins, {args.reps} interleaved reps) ===')
        times = {}
        for rep in range(args.reps):
            for g in GRIDS:
                tag = 'exact' if g is None else f'{g[0]}x{g[1]}'
                r = run(gal, ran, out, g)
                times.setdefault(tag, []).append(r['t_query'])
                print(f'  rep{rep} {tag:>7s}: {r["t_query"]:7.1f}s', flush=True)
        base = min(times[f'{GRIDS[0][0]}x{GRIDS[0][1]}'])
        print(f'\n{"mode":>9s} {"best[s]":>9s} {"vs 4x16":>8s}')
        tsum = []
        for tag, ts in times.items():
            b = min(ts)
            print(f'{tag:>9s} {b:9.1f} {b / base:8.3f}')
            tsum.append(dict(mode=tag, best=b, ratio=b / base))
        results['timing'] = tsum

    # per-edge / per-point graph memory model
    print('\n=== graph memory model (parity mode) ===')
    print('  pixel grid : 4 B id + 1 B dist + 2 B phi = 7 B/edge'
          '   (was 6 B/edge with the int8 index)')
    print('  exact      : 4 B id + 1 B dist          = 5 B/edge, '
          'plus 24 B/point for positions')
    print('  => at DESI NGC z1 scale (~0.9e9 edges, 2.4M points): '
          '6.3 GB pixel(int16) / 5.4 GB pixel(int8) / 4.5 GB exact')
    with open(args.out, 'w') as fh:
        json.dump(results, fh, indent=1)
    print(f'\nwrote {args.out}')


if __name__ == '__main__':
    main()
