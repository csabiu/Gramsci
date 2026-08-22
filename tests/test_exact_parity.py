"""Correctness tests + benchmark for -exactparity (exact signed volume).

1. Injection recovery: flattened / sheared / scalene chiral tetrahedra with
   random orientations.  Pixelized mode is attenuated (eta ~ 0.84 / 0.92 / 1);
   exact mode must recover ~1.000 for every shape.
2. CPU == GPU in exact mode (identical output files, unit weights).
3. GPU chunked (out-of-core) path == resident path in exact mode.
4. Benchmark pixel vs exact, GPU and CPU, interleaved runs.

Usage: python tests/test_exact_parity.py [--quick]
"""
import os
import subprocess
import sys
import time

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.join(ROOT, 'tests'))
GPU = os.path.join(ROOT, 'bin', 'gramsci_gpu')
CPU = os.path.join(ROOT, 'bin', 'gramsci')
WORK = os.environ.get('WORK', '/tmp/exactparity_test')
os.makedirs(WORK, exist_ok=True)

from generate_test_fields import generate_chiral_tetra_shapes
from parity_shapes import (SHAPES, canonicalise, edge_lengths,
                            choose_binning)

N_STRUCT, N_RAND = 500, 20000


# The coarse grid is named explicitly: it is no longer the default, and the
# point of test 1 is to contrast it with the exact evaluation.
COARSE = ['-ntheta', '4', '-nphi', '16']


def run(binary, gal, ran, out, rmin, rmax, nb, exact, env=None):
    cmd = [binary, '-gal', gal, '-ran', ran, '-rmin', str(rmin),
           '-rmax', str(rmax), '-nbins', str(nb), '-nmu', '1',
           '-out', out, '-4pcfp'] + (['-exactparity'] if exact else COARSE)
    e = dict(os.environ)
    if env:
        e.update(env)
    t0 = time.time()
    res = subprocess.run(cmd, capture_output=True, text=True, env=e)
    assert res.returncode == 0, res.stdout[-800:]
    # query time from the log if present, else wall time
    tq = time.time() - t0
    for line in res.stdout.splitlines():
        if 'Querying the graph took' in line:
            tq = float(line.split()[-2])
    return tq


def load(out):
    n_skip = 0
    with open(out) as fh:
        for line in fh:
            s = line.strip()
            if s and (s[0].isdigit() or s[0] in '+-'):
                break
            n_skip += 1
    return np.loadtxt(out, skiprows=n_skip)


def recovery(shape, scale, x, seed, exact, binary=GPU):
    v = SHAPES[shape] * scale
    edges = edge_lengths(v)
    rmin, rmax, nb = choose_binning(edges)
    diam = edges.max()
    spacing = diam + rmax + 10.0
    box = max(500.0, 8.0 * spacing + 1.0)
    gal = os.path.join(WORK, 'ep.gal')
    ran = os.path.join(WORK, 'ep.ran')
    out = os.path.join(WORK, f'ep_{int(exact)}.4pcfp')
    generate_chiral_tetra_shapes(N_STRUCT, box, N_RAND, seed, gal, ran,
                                 v, handedness=float(x), spacing=spacing)
    run(binary, gal, ran, out, rmin, rmax, nb, exact)
    d = load(out)
    w = (rmax - rmin) / nb
    enc = [int((e - rmin) // w) + 1 for e in edges]
    canon, eps = canonicalise(enc)
    row_bins = np.rint((d[:, 0:12:2] - rmin) / w).astype(int) + 1
    match = [i for i in range(len(d))
             if canonicalise(list(row_bins[i]))[0] == canon]
    assert len(match) == 1
    n_gal = 4 * N_STRUCT
    return eps * d[match[0], 15] * n_gal ** 4 / N_STRUCT, out


def main():
    quick = '--quick' in sys.argv

    print('=== 1. injection recovery: coarse 4x16 grid vs exact (GPU) ===')
    for shape in ('flattened', 'sheared', 'scalene'):
        for exact in (False, True):
            vals = [recovery(shape, 15.0, 1.0, s, exact)[0]
                    for s in (1, 2, 3)]
            print(f'  {shape:9s} {"exact" if exact else "4x16":5s}: '
                  f'{np.mean(vals):+.4f} +- {np.std(vals):.4f}')

    print('=== 2. CPU == GPU (exact mode, flattened) ===')
    _, out_g = recovery('flattened', 15.0, 1.0, 1, True, GPU)
    d_g = load(out_g)
    _, out_c = recovery('flattened', 15.0, 1.0, 1, True, CPU)
    d_c = load(out_c)
    same = np.array_equal(d_g, d_c)
    print(f'  identical output: {same}'
          + ('' if same else f'  max|diff|={np.abs(d_g - d_c).max():.3e}'))

    print('=== 3. GPU chunked == resident (exact mode) ===')
    v = SHAPES['flattened'] * 15.0
    edges = edge_lengths(v)
    rmin, rmax, nb = choose_binning(edges)
    spacing = edges.max() + rmax + 10.0
    box = max(500.0, 8.0 * spacing + 1.0)
    gal = os.path.join(WORK, 'ep.gal')
    ran = os.path.join(WORK, 'ep.ran')
    generate_chiral_tetra_shapes(N_STRUCT, box, N_RAND, 1, gal, ran, v,
                                 handedness=1.0, spacing=spacing)
    o1 = os.path.join(WORK, 'res.4pcfp')
    o2 = os.path.join(WORK, 'chk.4pcfp')
    run(GPU, gal, ran, o1, rmin, rmax, nb, True)
    run(GPU, gal, ran, o2, rmin, rmax, nb, True,
        env={'GRAMSCI_GPU_WIN_EDGES': '200000'})
    same = np.array_equal(load(o1), load(o2))
    print(f'  identical output: {same}')

    if quick:
        return

    print('=== 4. benchmark pixel vs exact (uniform box, r 20-65, 5 bins) ===')
    rng = np.random.default_rng(99)
    for tag, binary, npts in (('GPU', GPU, 200_000), ('CPU', CPU, 100_000)):
        ng, nr = npts // 3, npts - npts // 3
        gal = os.path.join(WORK, 'bench.gal')
        ran = os.path.join(WORK, 'bench.ran')
        np.savetxt(gal, np.column_stack(
            [rng.random((ng, 3)) * 1000.0, np.ones(ng)]), fmt='%.6e')
        np.savetxt(ran, np.column_stack(
            [rng.random((nr, 3)) * 1000.0, np.ones(nr)]), fmt='%.6e')
        out = os.path.join(WORK, 'bench.4pcfp')
        times = {False: [], True: []}
        for rep in range(2):
            for exact in (False, True):
                tq = run(binary, gal, ran, out, 20.0, 65.0, 5, exact)
                times[exact].append(tq)
        tp = min(times[False])
        te = min(times[True])
        print(f'  {tag} ({npts // 1000}k pts): pixel {tp:8.1f}s   '
              f'exact {te:8.1f}s   ratio {te / tp:.3f}')


if __name__ == '__main__':
    main()
