"""Delete-one jackknife tests: internal angular partition + exact
brute-force verification of every statistic's realisations.

Strategy: run gramsci with -njk (internal angular partition), read back the
region labels it wrote (<out>.jkgal/.jkran), then recompute each mode's
delete-one counts by brute force in Python with the same normalized signed
weights and the same strict rmin^2 < d^2 < rmax^2 edge window and bin
formula as create_graph.  N_m must match the .jk rows to file precision:
the jackknife identity N_m = N_total - N_touching(m) is exact for
delete-one counts, so any dedup or accumulation bug shows up directly.

The 4PCF check ports the signed bintable6 fill (seed = lexicographic
orbit minimum, entry sign = seed parity * permutation parity) rather than
canonicalising per tuple, and the chiral_4pcf mask: a canonical tuple fixed
by some odd vertex permutation is achiral at the bin level, so the kernel
zeroes its odd channel (the labelling-dependent sign it used to get was an
artefact).  Run with -exactparity so Python can reproduce sign_V from the
positions alone.
"""
import itertools
import os
import subprocess

import numpy as np

from parity_shapes import _S4_TABLE

VOL_DEGEN_TOL = 1.0e-9  # keep in sync with src/query_4pcf_module.F90


def run(cmd):
    print(cmd)
    subprocess.check_call(cmd, shell=True)


def _cleanup(*files):
    for f in files:
        if os.path.exists(f):
            os.remove(f)


def _cleanup_glob(prefix):
    import glob
    for f in glob.glob(prefix + '*'):
        os.remove(f)


# ---------------------------------------------------------------------------
# Mock catalogues: full-sky shell around the observer at the origin, so the
# angular partition exercises both sin(dec) bands and phi slices.
# ---------------------------------------------------------------------------

def _shell_uniform(rng, n, r_in, r_out):
    """n points uniform in the spherical shell [r_in, r_out]."""
    u = rng.random(n)
    r = (r_in**3 + u * (r_out**3 - r_in**3)) ** (1.0 / 3.0)
    mu = 2.0 * rng.random(n) - 1.0
    phi = 2.0 * np.pi * rng.random(n)
    s = np.sqrt(1.0 - mu * mu)
    return np.column_stack([r * s * np.cos(phi), r * s * np.sin(phi), r * mu])


def _shell_clustered(rng, n, r_in, r_out, n_blob=30, sigma=12.0):
    """Clustered points: Gaussian blobs around uniform shell centers."""
    centers = _shell_uniform(rng, n_blob, r_in + sigma, r_out - sigma)
    idx = rng.integers(0, n_blob, n)
    pts = centers[idx] + rng.normal(0.0, sigma, (n, 3))
    return pts


def _write_cat(path, pts):
    np.savetxt(path, pts, fmt='%.10e')  # x y z (weight defaults to 1)


# ---------------------------------------------------------------------------
# Brute-force machinery (replicates create_graph's window and binning)
# ---------------------------------------------------------------------------

def _adjacency(pts, rmin, rmax, nbins):
    """(adj, bins): strict rmin^2 < d^2 < rmax^2 window; clamped linear bins."""
    d2 = ((pts[:, None, :] - pts[None, :, :]) ** 2).sum(axis=2)
    adj = (d2 > rmin * rmin) & (d2 < rmax * rmax)
    np.fill_diagonal(adj, False)
    dr = (rmax - rmin) / nbins
    inv = 1.0 / dr
    with np.errstate(invalid='ignore'):
        b = np.floor((np.sqrt(d2) - rmin) * inv).astype(int) + 1
    bins = np.clip(b, 1, nbins)
    return adj, bins


def _weights(n_gal, n_ran):
    w = np.empty(n_gal + n_ran)
    w[:n_gal] = 1.0 / n_gal
    w[n_gal:] = -1.0 / n_ran
    return w


def _touch_add(acc, key, regions, val):
    """Add val to acc[key][r] for each DISTINCT region r > 0 in regions."""
    seen = set()
    for r in regions:
        if r > 0 and r not in seen:
            seen.add(r)
            acc[key][r - 1] += val


def _build_bintable6(nbins):
    """Signed 6D config table, exact port of create_4pcf_binlookup: tuples
    visited in ascending lexicographic order, so each orbit is seeded by its
    canonical (lex-min) member with seed parity +1; each orbit member's sign
    is the parity of the FIRST table permutation that reaches it.  Also
    returns chiral[cid - 1] = 0 when the canonical tuple is invariant under
    an odd permutation (achiral binned configuration, odd channel := 0)."""
    table = {}
    canon = []
    chiral = []
    for tup in itertools.product(range(1, nbins + 1), repeat=6):
        if tup in table:
            continue
        canon.append(tup)
        cid = len(canon)
        achiral = False
        for ep, sign in _S4_TABLE:
            perm = tuple(tup[ep[j] - 1] for j in range(6))
            if perm not in table:
                table[perm] = sign * cid
            if sign < 0 and perm == tup:
                achiral = True
        chiral.append(0 if achiral else 1)
    return table, canon, chiral


def _check_jkcov(path, jk_rows, est_col, njk, sigma=None):
    """Verify a .jkcov file: symmetry, agreement with the covariance built
    directly from the .jk realisations, and (optionally) that its diagonal
    matches the .jkerr sigma^2 -- all to file precision."""
    C = np.loadtxt(path, comments='#', ndmin=2)
    n = C.shape[0]
    assert C.shape == (n, n), f'{path}: not square'
    assert np.allclose(C, C.T, rtol=1e-6, atol=1e-30), f'{path}: not symmetric'
    x = jk_rows[:, est_col].reshape(n, njk)
    d = x - x.mean(axis=1, keepdims=True)
    Cref = (njk - 1) / njk * (d @ d.T)
    scale = np.abs(Cref).max()
    assert np.allclose(C, Cref, rtol=2e-6, atol=1e-6 * scale), \
        f'{path}: does not match covariance of the .jk realisations'
    if sigma is not None:
        assert np.allclose(np.sqrt(np.diag(C)), sigma,
                           rtol=3e-6, atol=1e-6 * np.abs(sigma).max() + 1e-30), \
            f'{path}: diagonal does not match .jkerr sigma^2'


def _assert_rows(name, got, want, scale):
    """Compare against file values: e14.7 rounding + summation-order slack."""
    tol = 1.0e-6 * np.abs(want) + 1.0e-9 * scale + 1.0e-30
    bad = np.abs(got - want) > tol
    assert not np.any(bad), (
        f"{name}: {bad.sum()} mismatches, worst "
        f"|diff|={np.max(np.abs(got - want)):.3e} "
        f"(first: got {got[bad][0]!r} want {want[bad][0]!r})")


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_internal_partition(bindir):
    """Internal angular partition: labels written, in range, balanced."""
    print('\n=== Test: internal angular jackknife partition ===')
    rng = np.random.default_rng(42)
    n_gal, n_ran, njk = 2000, 6000, 8
    gal, ran, out = 'tmp_jkpart.gal', 'tmp_jkpart.ran', 'tmp_jkpart.out'
    _write_cat(gal, _shell_clustered(rng, n_gal, 60.0, 130.0))
    _write_cat(ran, _shell_uniform(rng, n_ran, 60.0, 130.0))
    gramsci = os.path.join(bindir, 'gramsci')

    run(f"{gramsci} -gal {gal} -ran {ran} -rmin 10 -rmax 35 -nbins 3 "
        f"-nmu 1 -out {out} -2pcf -njk {njk}")

    lg = np.loadtxt(out + '.jkgal', dtype=int)
    lr = np.loadtxt(out + '.jkran', dtype=int)
    assert len(lg) == n_gal and len(lr) == n_ran, 'label file length mismatch'
    assert lg.min() >= 1 and lg.max() <= njk, 'galaxy label out of range'
    assert lr.min() >= 1 and lr.max() <= njk, 'random label out of range'
    # Equal-count quantiles on the randoms: perfectly balanced up to rounding
    # at the quantile boundaries.
    cnt_r = np.bincount(lr, minlength=njk + 1)[1:]
    assert cnt_r.min() > 0.9 * n_ran / njk, f'random counts unbalanced: {cnt_r}'
    assert cnt_r.max() < 1.1 * n_ran / njk, f'random counts unbalanced: {cnt_r}'
    # Purely angular: within every region the radial distribution must span
    # the full shell (a radial partition would truncate it).
    r = np.linalg.norm(np.loadtxt(ran), axis=1)
    for m in range(1, njk + 1):
        rm = r[lr == m]
        assert rm.min() < 70.0 and rm.max() > 120.0, \
            f'region {m} does not span the radial shell: angular-only violated'
    print('  partition checks passed')
    _cleanup(gal, ran)
    _cleanup_glob(out)
    print('  Test PASSED')


def test_jackknife_2pcf_3pcf_equi(bindir):
    """Brute-force delete-one identity for 2PCF, 3PCF and equilateral."""
    print('\n=== Test: jackknife realisations vs brute force (2/3PCF, equi) ===')
    rng = np.random.default_rng(7)
    n_gal, n_ran, njk = 400, 800, 5
    rmin, rmax, nbins = 10.0, 35.0, 3
    gal, ran, out = 'tmp_jk23.gal', 'tmp_jk23.ran', 'tmp_jk23.out'
    pts = np.vstack([_shell_clustered(rng, n_gal, 60.0, 130.0),
                     _shell_uniform(rng, n_ran, 60.0, 130.0)])
    _write_cat(gal, pts[:n_gal])
    _write_cat(ran, pts[n_gal:])
    gramsci = os.path.join(bindir, 'gramsci')

    run(f"{gramsci} -gal {gal} -ran {ran} -rmin {rmin} -rmax {rmax} "
        f"-nbins {nbins} -nmu 1 -out {out} -2pcf -3pcf -equi -njk {njk}")

    reg = np.concatenate([np.loadtxt(out + '.jkgal', dtype=int),
                          np.loadtxt(out + '.jkran', dtype=int)])
    w = _weights(n_gal, n_ran)
    adj, bins = _adjacency(pts, rmin, rmax, nbins)
    n = len(pts)
    is_ran = np.arange(n) >= n_gal

    # ---- brute-force accumulators ----
    NN = np.zeros(nbins + 1); RR = np.zeros(nbins + 1)
    NNjk = {b: np.zeros(njk) for b in range(1, nbins + 1)}
    RRjk = {b: np.zeros(njk) for b in range(1, nbins + 1)}
    triples = list(itertools.combinations_with_replacement(range(1, nbins + 1), 3))
    t_of = {t: i for i, t in enumerate(triples)}
    NNN = np.zeros(len(triples)); RRR = np.zeros(len(triples))
    NNNjk = {t: np.zeros(njk) for t in range(len(triples))}
    RRRjk = {t: np.zeros(njk) for t in range(len(triples))}
    eNNN = np.zeros(nbins + 1); eRRR = np.zeros(nbins + 1)
    eNNNjk = {b: np.zeros(njk) for b in range(1, nbins + 1)}
    eRRRjk = {b: np.zeros(njk) for b in range(1, nbins + 1)}
    NN_scale = np.zeros(nbins + 1)
    NNN_scale = np.zeros(len(triples))

    idx = np.arange(n)
    for i in range(n):
        nb_i = idx[adj[i] & (idx > i)]
        for j in nb_i:
            wp = w[i] * w[j]
            b = bins[i, j]
            NN[b] += wp
            NN_scale[b] += abs(wp)
            _touch_add(NNjk, b, (reg[i], reg[j]), wp)
            if is_ran[i] and is_ran[j]:
                RR[b] += wp
                _touch_add(RRjk, b, (reg[i], reg[j]), wp)
            # triangles with k > j adjacent to both
            for k in idx[(idx > j) & adj[i] & adj[j]]:
                w3 = wp * w[k]
                t = t_of[tuple(sorted((bins[i, j], bins[i, k], bins[j, k])))]
                NNN[t] += w3
                NNN_scale[t] += abs(w3)
                _touch_add(NNNjk, t, (reg[i], reg[j], reg[k]), w3)
                rrr = is_ran[i] and is_ran[j] and is_ran[k]
                if rrr:
                    RRR[t] -= w3          # 3PCF stores RRR with a minus sign
                    _touch_add(RRRjk, t, (reg[i], reg[j], reg[k]), -w3)
                if bins[i, j] == bins[i, k] == bins[j, k]:
                    b3 = bins[i, j]
                    eNNN[b3] += w3
                    _touch_add(eNNNjk, b3, (reg[i], reg[j], reg[k]), w3)
                    if rrr:
                        eRRR[b3] -= w3
                        _touch_add(eRRRjk, b3, (reg[i], reg[j], reg[k]), -w3)

    # ---- compare with the .jk files ----
    jk2 = np.loadtxt(out + '.2pcf.jk', comments='#')
    for b in range(1, nbins + 1):
        rows = jk2[(b - 1) * njk: b * njk]
        assert np.all(rows[:, 4].astype(int) == np.arange(1, njk + 1))
        _assert_rows(f'2pcf NN bin{b}', NN[b] - NNjk[b], rows[:, 5], NN_scale[b])
        _assert_rows(f'2pcf RR bin{b}', RR[b] - RRjk[b], rows[:, 6], NN_scale[b])

    jk3 = np.loadtxt(out + '.3pcf.jk', comments='#')
    for t in range(len(triples)):
        rows = jk3[t * njk: (t + 1) * njk]
        _assert_rows(f'3pcf NNN cfg{t}', NNN[t] - NNNjk[t], rows[:, 7],
                     NNN_scale[t])
        _assert_rows(f'3pcf RRR cfg{t}', RRR[t] - RRRjk[t], rows[:, 8],
                     NNN_scale[t])

    jke = np.loadtxt(out + '.equi.jk', comments='#')
    for b in range(1, nbins + 1):
        rows = jke[(b - 1) * njk: b * njk]
        _assert_rows(f'equi NNN bin{b}', eNNN[b] - eNNNjk[b], rows[:, 5],
                     NNN_scale.sum())
        _assert_rows(f'equi RRR bin{b}', eRRR[b] - eRRRjk[b], rows[:, 6],
                     NNN_scale.sum())

    # jkerr sigma must match the realisations
    err2 = np.loadtxt(out + '.2pcf.jkerr', comments='#')
    for b in range(1, nbins + 1):
        xis = jk2[(b - 1) * njk: b * njk, 7]
        sig = np.sqrt((njk - 1) / njk * np.sum((xis - xis.mean()) ** 2))
        assert abs(err2[b - 1, 5] - sig) <= 1e-5 * sig + 1e-9, \
            f'2pcf jkerr sigma mismatch bin {b}'

    # covariance files: symmetric, match the realisations, diag == sigma^2
    _check_jkcov(out + '.2pcf.jkcov', jk2, 7, njk, sigma=err2[:, 5])
    err3 = np.loadtxt(out + '.3pcf.jkerr', comments='#')
    _check_jkcov(out + '.3pcf.jkcov', jk3, 9, njk, sigma=err3[:, 7])
    erre = np.loadtxt(out + '.equi.jkerr', comments='#')
    _check_jkcov(out + '.equi.jkcov', jke, 7, njk, sigma=erre[:, 5])
    print('  .jkcov files verified against realisations and .jkerr')

    print('  2PCF/3PCF/equi delete-one identities verified')
    _cleanup(gal, ran)
    _cleanup_glob(out)
    print('  Test PASSED')


def test_jackknife_4pcf_parity(bindir):
    """Brute-force delete-one identity for the 4PCF and parity channels."""
    print('\n=== Test: jackknife realisations vs brute force (4PCF parity) ===')
    rng = np.random.default_rng(11)
    n_gal, n_ran, njk = 150, 300, 4
    # 3 bins: 65 configurations of which 21 are chiral (2 bins leave only
    # 1 of 11), at the same quadruplet-enumeration cost.
    rmin, rmax, nbins = 10.0, 40.0, 3
    gal, ran, out = 'tmp_jk4.gal', 'tmp_jk4.ran', 'tmp_jk4.out'
    pts = np.vstack([_shell_clustered(rng, n_gal, 60.0, 120.0, sigma=15.0),
                     _shell_uniform(rng, n_ran, 60.0, 120.0)])
    _write_cat(gal, pts[:n_gal])
    _write_cat(ran, pts[n_gal:])
    gramsci = os.path.join(bindir, 'gramsci')

    run(f"{gramsci} -gal {gal} -ran {ran} -rmin {rmin} -rmax {rmax} "
        f"-nbins {nbins} -nmu 1 -out {out} -4pcfp -exactparity -njk {njk}")

    reg = np.concatenate([np.loadtxt(out + '.jkgal', dtype=int),
                          np.loadtxt(out + '.jkran', dtype=int)])
    w = _weights(n_gal, n_ran)
    adj, bins = _adjacency(pts, rmin, rmax, nbins)
    n = len(pts)
    is_ran = np.arange(n) >= n_gal
    table, canon, chiral = _build_bintable6(nbins)
    ncfg = len(canon)
    print(f'  {ncfg - sum(chiral)} of {ncfg} configs achiral (odd channel 0)')

    N4 = np.zeros(ncfg + 1); N4o = np.zeros(ncfg + 1)
    R4 = np.zeros(ncfg + 1)
    N4jk = {c: np.zeros(njk) for c in range(1, ncfg + 1)}
    N4ojk = {c: np.zeros(njk) for c in range(1, ncfg + 1)}
    R4jk = {c: np.zeros(njk) for c in range(1, ncfg + 1)}
    scale = np.zeros(ncfg + 1)

    idx = np.arange(n)
    nquad = 0
    for i in range(n):
        nb_i = idx[adj[i] & (idx > i)]
        for a, j in enumerate(nb_i):
            for k in nb_i[a + 1:]:
                if not adj[j, k]:
                    continue
                common = idx[(idx > k) & adj[i] & adj[j] & adj[k]]
                for l in common:
                    tup = (bins[i, j], bins[i, k], bins[i, l],
                           bins[j, k], bins[j, l], bins[k, l])
                    raw = table[tup]
                    cid = abs(raw)
                    flip = (1 if raw > 0 else -1) * chiral[cid - 1]
                    u1 = pts[j] - pts[i]; u1 = u1 / np.linalg.norm(u1)
                    u2 = pts[k] - pts[i]; u2 = u2 / np.linalg.norm(u2)
                    u3 = pts[l] - pts[i]; u3 = u3 / np.linalg.norm(u3)
                    vol = np.dot(u1, np.cross(u2, u3))
                    if abs(vol) < VOL_DEGEN_TOL:
                        sv = 0
                    else:
                        sv = 1 if vol > 0 else -1
                    w4 = w[i] * w[j] * w[k] * w[l]
                    so = flip * sv * w4
                    nquad += 1
                    N4[cid] += w4
                    N4o[cid] += so
                    scale[cid] += abs(w4)
                    regs = (reg[i], reg[j], reg[k], reg[l])
                    _touch_add(N4jk, cid, regs, w4)
                    _touch_add(N4ojk, cid, regs, so)
                    if is_ran[i] and is_ran[j] and is_ran[k] and is_ran[l]:
                        R4[cid] += w4
                        _touch_add(R4jk, cid, regs, w4)
    print(f'  brute force enumerated {nquad} quadruplets, {ncfg} configs')
    assert nquad > 1000, 'test catalogue too sparse to be meaningful'

    # single-mode run: outputs land at <out> / <out>.jk (no .4pcfp suffix)
    jk4 = np.loadtxt(out + '.jk', comments='#')
    assert jk4.shape[0] == ncfg * njk, '4pcfp.jk row count mismatch'
    for c in range(1, ncfg + 1):
        rows = jk4[(c - 1) * njk: c * njk]
        assert np.all(rows[:, 12].astype(int) == np.arange(1, njk + 1))
        _assert_rows(f'4pcfp NNNN cfg{c}', N4[c] - N4jk[c], rows[:, 13], scale[c])
        _assert_rows(f'4pcfp RRRR cfg{c}', R4[c] - R4jk[c], rows[:, 14], scale[c])
        _assert_rows(f'4pcfp NNNN_odd cfg{c}', N4o[c] - N4ojk[c], rows[:, 16],
                     scale[c])
    # conn_odd column must equal zeta_odd (no disconnected term in the odd
    # channel), and the full-sample table row order must match canon order.
    _assert_rows('4pcfp conn_odd == zeta_odd', jk4[:, 17], jk4[:, 20],
                 np.abs(jk4[:, 17]).max())

    # covariance files: conn_even + the odd channel
    err4 = np.loadtxt(out + '.jkerr', comments='#')
    _check_jkcov(out + '.jkcov', jk4, 19, njk, sigma=err4[:, 17])
    _check_jkcov(out + '.jkcov_odd', jk4, 20, njk, sigma=err4[:, 15])
    print('  4pcfp .jkcov / .jkcov_odd verified')

    full = np.loadtxt(out, comments='#')
    for c in range(1, ncfg + 1):
        _assert_rows(f'4pcfp full NNNN cfg{c}', np.array([N4[c]]),
                     np.array([full[c - 1, 12]]), scale[c])

    print('  4PCF parity delete-one identities verified')
    _cleanup(gal, ran)
    _cleanup_glob(out)
    print('  Test PASSED')


def test_jackknife_combined_matches_single(bindir):
    """Combined-mode jackknife files match single-mode runs byte-for-byte
    in their data rows (guards the per-mode re-zeroing of N2jk/N3jk)."""
    print('\n=== Test: combined-mode jackknife == single-mode ===')
    rng = np.random.default_rng(3)
    gal, ran = 'tmp_jkc.gal', 'tmp_jkc.ran'
    _write_cat(gal, _shell_clustered(rng, 300, 60.0, 120.0))
    _write_cat(ran, _shell_uniform(rng, 600, 60.0, 120.0))
    gramsci = os.path.join(bindir, 'gramsci')
    args = f"-gal {gal} -ran {ran} -rmin 10 -rmax 35 -nbins 3 -nmu 1 -njk 4"

    run(f"{gramsci} {args} -out tmp_jkc_comb.out -2pcf -3pcf")
    run(f"{gramsci} {args} -out tmp_jkc_2.out -2pcf")
    run(f"{gramsci} {args} -out tmp_jkc_3.out -3pcf")

    def data_rows(path):
        with open(path) as f:
            return [l for l in f if not l.startswith('#')]

    assert data_rows('tmp_jkc_comb.out.2pcf.jk') == data_rows('tmp_jkc_2.out.jk'), \
        '2pcf jackknife differs between combined and single-mode runs'
    assert data_rows('tmp_jkc_comb.out.3pcf.jk') == data_rows('tmp_jkc_3.out.jk'), \
        '3pcf jackknife differs between combined and single-mode runs'
    print('  combined == single verified')
    _cleanup(gal, ran)
    for p in ('tmp_jkc_comb.out', 'tmp_jkc_2.out', 'tmp_jkc_3.out'):
        _cleanup_glob(p)
    print('  Test PASSED')


def test_jackknife_validation(bindir):
    """-njk rejects: analytic/-box-internal partitions and missing randoms."""
    print('\n=== Test: jackknife CLI validation ===')
    rng = np.random.default_rng(5)
    gal, ran = 'tmp_jkv.gal', 'tmp_jkv.ran'
    _write_cat(gal, rng.random((100, 3)) * 100.0)
    _write_cat(ran, rng.random((200, 3)) * 100.0)
    gramsci = os.path.join(bindir, 'gramsci')

    def must_fail(extra, why):
        r = subprocess.run(
            f"{gramsci} -gal {gal} -rmin 5 -rmax 20 -nbins 3 -nmu 1 "
            f"-out tmp_jkv.out -2pcf {extra}",
            shell=True, capture_output=True, text=True)
        assert r.returncode != 0, f'expected failure: {why}'

    must_fail('-njk 4', 'njk without -ran')
    must_fail(f'-ran {ran} -njk 1', 'njk = 1')
    must_fail(f'-ran {ran} -box 100 -njk 4',
              'internal angular partition in a periodic box')
    must_fail('-box 100 -njk 4', 'njk with analytic randoms')
    print('  all invalid combinations rejected')
    _cleanup(gal, ran, 'tmp_jkv.out')
    print('  Test PASSED')


def run_all(bindir):
    test_internal_partition(bindir)
    test_jackknife_2pcf_3pcf_equi(bindir)
    test_jackknife_4pcf_parity(bindir)
    test_jackknife_combined_matches_single(bindir)
    test_jackknife_validation(bindir)


if __name__ == '__main__':
    run_all(os.path.join('..', 'bin'))
