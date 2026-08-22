"""Tetrahedron shape families and the kernel's S4 canonicalisation.

Shared by the parity tests.  `canonicalise` is an exact port of the Fortran
canonicalisation in src/query_4pcf_module.F90 (edge order 12,13,14,23,24,34;
lexicographic minimum over the 24 vertex permutations, updating only on a
strictly smaller tuple so that ties keep the earlier parity), so a test can
predict the sign the kernel will attach to an injected configuration.
"""
import numpy as np

_reg = np.array([[0.0, 0.0, 0.0],
                 [1.0, 0.0, 0.0],
                 [0.5, np.sqrt(3) / 2, 0.0],
                 [0.5, np.sqrt(3) / 6, np.sqrt(2.0 / 3.0)]])
_sca = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
                 [0.5, 1.589899, 0.0],
                 [-0.833333, 0.297014, 1.412486]])
_shear = np.array([[1.0, 0.5, 0.3], [0.0, 1.0, 0.5], [0.0, 0.0, 1.0]])
SHAPES = {
    'regular':   _reg,
    'scalene':   _sca,
    'elongated': _sca * np.array([1.0, 1.0, 2.2]) / 1.5,
    'flattened': _sca * np.array([1.0, 1.0, 0.28]),
    'sheared':   _reg @ _shear.T,
}

# ---------------------------------------------------------------------------
# Exact port of the Fortran S4 canonicalisation (query_4pcf_module.F90):
# edge order (12,13,14,23,24,34); the 24 edge permutations and parities,
# lexicographic minimum with first-strictly-smaller update (ties keep the
# earlier parity).
# ---------------------------------------------------------------------------
# Verbatim transcription of S4_EDGE_PERMS / S4_PARITY (1-based edge
# indices; new(j) = old(perm(j))) in the Fortran table order, so that the
# first-strictly-smaller tie rule reproduces the kernel's convention
# exactly for every tuple, symmetric ones included.
_S4_TABLE = [
    ([1, 2, 3, 4, 5, 6], +1), ([1, 4, 5, 2, 3, 6], -1),
    ([4, 2, 6, 1, 5, 3], -1), ([5, 6, 3, 4, 1, 2], -1),
    ([2, 1, 3, 4, 6, 5], -1), ([3, 2, 1, 6, 5, 4], -1),
    ([1, 3, 2, 5, 4, 6], -1), ([1, 5, 4, 3, 2, 6], +1),
    ([6, 2, 4, 3, 5, 1], +1), ([6, 5, 3, 4, 2, 1], +1),
    ([4, 1, 5, 2, 6, 3], +1), ([5, 4, 1, 6, 3, 2], +1),
    ([2, 4, 6, 1, 3, 5], +1), ([4, 6, 2, 5, 1, 3], +1),
    ([3, 6, 5, 2, 1, 4], +1), ([5, 3, 6, 1, 4, 2], +1),
    ([2, 3, 1, 6, 4, 5], +1), ([3, 1, 2, 5, 6, 4], +1),
    ([4, 5, 1, 6, 2, 3], -1), ([5, 1, 4, 3, 6, 2], -1),
    ([6, 4, 2, 5, 3, 1], -1), ([2, 6, 4, 3, 1, 5], -1),
    ([6, 3, 5, 2, 4, 1], -1), ([3, 5, 6, 1, 2, 4], -1),
]


def canonicalise(bins6):
    """(canonical tuple, parity) the kernel assigns to this tuple's
    orbit, replicating create_4pcf_binlookup exactly (identity first,
    update only on strictly smaller, ties keep the earlier parity)."""
    bins6 = list(bins6)
    canon, best = bins6[:], +1
    for ep, sign in _S4_TABLE[1:]:
        cand = [bins6[ep[j] - 1] for j in range(6)]
        if cand < canon:
            canon, best = cand, sign
    return tuple(canon), best


def canonical_parity(bins6):
    return canonicalise(bins6)[1]


def edge_lengths(v):
    """The six pairwise distances in GRAMSCI edge order (12,13,14,23,24,34)."""
    pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    return np.array([np.linalg.norm(v[b] - v[a]) for a, b in pairs])


def choose_binning(edges):
    """(rmin, rmax, nbins) placing every edge well inside a bin."""
    emin, emax = edges.min(), edges.max()
    best = None
    for nb in range(2, 7):
        for rminf in (0.25, 0.4, 0.55, 0.7):
            rmin = max(1.0, rminf * emin)
            if rmin > 0.9 * emin:
                continue
            for rmaxf in (1.1, 1.2, 1.3, 1.45):
                rmax = rmaxf * emax
                w = (rmax - rmin) / nb
                pos = (edges - rmin) / w
                margin = np.min(np.minimum(pos - np.floor(pos),
                                           np.ceil(pos) - pos))
                if best is None or margin > best[0]:
                    best = (margin, rmin, rmax, nb)
    margin, rmin, rmax, nb = best
    if margin < 0.15:
        print(f'  WARNING: bin margin only {margin:.2f}')
    return rmin, rmax, nb


def run_one(shape, scale, x, seed):
    v = SHAPES[shape] * scale
    edges = edge_lengths(v)
    rmin, rmax, nb = choose_binning(edges)
    diam = edges.max()
    spacing = diam + rmax + 10.0
    box = max(500.0, 8.0 * spacing + 1.0)

    gal = os.path.join(WORK, 'pg.gal')
    ran = os.path.join(WORK, 'pg.ran')
    out = os.path.join(WORK, 'pg.4pcfp')
    generate_chiral_tetra_shapes(N_STRUCT, box, N_RAND, seed, gal, ran,
                                 v, handedness=float(x), spacing=spacing)
    cmd = [GPU, '-gal', gal, '-ran', ran, '-rmin', str(rmin),
           '-rmax', str(rmax), '-nbins', str(nb), '-nmu', '1',
           '-out', out, '-4pcfp',
           '-ntheta', str(N_THETA), '-nphi', str(N_PHI)]
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0, res.stdout[-500:]

    n_skip = 0
    with open(out) as fh:
        for line in fh:
            s = line.strip()
            if s and (s[0].isdigit() or s[0] in '+-'):
                break
            n_skip += 1
    d = np.loadtxt(out, skiprows=n_skip)
    w = (rmax - rmin) / nb
    enc = [int((e - rmin) // w) + 1 for e in edges]
    canon, eps = canonicalise(enc)
    # locate the injected configuration by its canonical bin tuple
    row_bins = np.rint((d[:, 0:12:2] - rmin) / w).astype(int) + 1
    match = [i for i in range(len(d))
             if canonicalise(list(row_bins[i]))[0] == canon]
    assert len(match) == 1, f'config match failed: {len(match)} rows'
    row = match[0]
    n_gal = 4 * N_STRUCT
    recovered = d[row, 15] * n_gal ** 4 / N_STRUCT  # NNNN_odd, count level
    return dict(shape=shape, scale=scale, x=x, seed=seed,
                injected=2 * x - 1, recovered=float(recovered),
                eps=eps, n_theta=N_THETA, n_phi=N_PHI,
                rmin=rmin, rmax=rmax, nbins=nb,
                edges=[round(float(e), 3) for e in edges])
