import os
import subprocess
import math
import numpy as np
from generate_test_fields import (generate_pairs, generate_regular_tetra,
                                  generate_chiral_tetra, generate_uniform_box,
                                  generate_aniso_box,
                                  RMIN, RMAX, NBINS, R0, BOX_SIZE, N_RAND, SEED)

def run(cmd):
    print(cmd)
    subprocess.check_call(cmd, shell=True)

def assert_close(val, expected, tol, msg):
    if abs(val - expected)/expected > tol:
        raise AssertionError(f"{msg}: got {val} expected {expected}")

# ---------------------------------------------------------------------------
# Helper functions for physics tests
# ---------------------------------------------------------------------------

def find_row_all_bins_equal(data, target_bin, rmin, rmax, nbins):
    """Find the output row where all 6 edge-bin pairs match target_bin."""
    dr = (rmax - rmin) / nbins
    target_lo = rmin + (target_bin - 1) * dr
    target_hi = rmin + target_bin * dr
    for row in data:
        all_match = True
        for k in range(6):
            if abs(row[2*k] - target_lo) > 0.1 or abs(row[2*k+1] - target_hi) > 0.1:
                all_match = False
                break
        if all_match:
            return row
    raise ValueError(f"No row found with all bins equal to bin {target_bin}")


def find_row_by_sorted_bins(data, sorted_bins, rmin, rmax, nbins):
    """Find the output row matching a specific sorted bin-tuple."""
    dr = (rmax - rmin) / nbins
    for row in data:
        bins = []
        for k in range(6):
            lo = row[2 * k]
            b = int(round((lo - rmin) / dr)) + 1
            bins.append(b)
        if bins == sorted_bins:
            return row
    raise ValueError(f"No row found matching sorted bins {sorted_bins}")


def _cleanup(*files):
    """Remove temporary files."""
    for f in files:
        if os.path.exists(f):
            os.remove(f)


# ---------------------------------------------------------------------------
# Physics tests
# ---------------------------------------------------------------------------

def test_pairs_2pcf(bindir):
    """Test 1: Isolated pairs validate 2PCF signal localization.

    100 pairs with separation r0=12.5 (bin 2 center) + randoms.
    The 2PCF should show a strong signal in the bin containing r0
    and be much weaker in other bins.
    """
    print('\n=== Test: Isolated pairs (2PCF) ===')

    gal_file = 'tmp_pairs.gal'
    ran_file = 'tmp_pairs.ran'
    out_2pcf = 'tmp_pairs_2pcf.out'
    gramsci = os.path.join(bindir, 'gramsci')

    generate_pairs(100, R0, BOX_SIZE, N_RAND, SEED, gal_file, ran_file)

    run(f"{gramsci} -gal {gal_file} -ran {ran_file} "
        f"-rmin {RMIN} -rmax {RMAX} -nbins {NBINS} -nmu 1 -out {out_2pcf} -2pcf")

    data_2pcf = np.loadtxt(out_2pcf, skiprows=1)
    # Columns: rmin rmax mu_min mu_max NN RR xi
    xi = data_2pcf[:, 6]
    print(f"  2PCF xi values: {xi}")

    # The signal bin (bin 2, containing r0=12.5) should have xi >> 1,
    # while off-bins should be much weaker.
    assert xi[1] > 1.0, f"Pairs: xi(bin 2) should be >> 1, got {xi[1]}"
    assert xi[1] > xi[0] + 1.0, f"Pairs: xi(bin 2) should dominate bin 1"
    assert xi[1] > xi[2] + 1.0, f"Pairs: xi(bin 2) should dominate bin 3"
    print("  2PCF checks passed")

    _cleanup(gal_file, ran_file, out_2pcf)
    print("  Test PASSED")


def test_regular_tetra_connected(bindir):
    """Test 2: Regular tetrahedra produce connected 4PCF signal.

    100 regular tetrahedra (all edges = r0 = 12.5) + 8000 randoms.
    At the all-equal-bin-2 configuration, zeta_conn should be positive
    (genuine 4-point correlations that cannot be explained by 2PCF products).
    """
    print('\n=== Test: Regular tetrahedra (connected 4PCF) ===')

    gal_file = 'tmp_regtetra.gal'
    ran_file = 'tmp_regtetra.ran'
    out_4pcf = 'tmp_regtetra_4pcf.out'
    gramsci = os.path.join(bindir, 'gramsci')

    generate_regular_tetra(100, R0, BOX_SIZE, N_RAND, SEED, gal_file, ran_file)

    run(f"{gramsci} -gal {gal_file} -ran {ran_file} "
        f"-rmin {RMIN} -rmax {RMAX} -nbins {NBINS} -nmu 1 -out {out_4pcf} -4pcf")

    data_4pcf = np.loadtxt(out_4pcf, skiprows=1)
    row = find_row_all_bins_equal(data_4pcf, target_bin=2,
                                  rmin=RMIN, rmax=RMAX, nbins=NBINS)

    zeta = row[14]
    zeta_disc = row[15]
    zeta_conn = row[16]
    print(f"  4PCF all-equal-bin-2: zeta={zeta:.6f}, "
          f"zeta_disc={zeta_disc:.6f}, zeta_conn={zeta_conn:.6f}")

    assert zeta_conn > 0.0, \
        f"Regular tetra: zeta_conn should be > 0, got {zeta_conn:.6f}"
    print("  Connected 4PCF check passed")

    _cleanup(gal_file, ran_file, out_4pcf)
    print("  Test PASSED")


def test_chiral_parity(bindir):
    """Test 3: Chiral tetrahedra validate parity-odd channel.

    3a) Left-only: 100 left-handed tetrahedra -> nonzero zeta_odd
    3b) Balanced: 50 left + 50 right -> zeta_odd cancels to ~ 0
    """
    print('\n=== Test: Chiral tetrahedra (parity-odd 4PCF) ===')

    gal_left = 'tmp_chiral_left.gal'
    ran_left = 'tmp_chiral_left.ran'
    out_left = 'tmp_chiral_left_4pcfp.out'

    gal_bal = 'tmp_chiral_bal.gal'
    ran_bal = 'tmp_chiral_bal.ran'
    out_bal = 'tmp_chiral_bal_4pcfp.out'

    gramsci = os.path.join(bindir, 'gramsci')

    # --- 3a: Left-only ---
    generate_chiral_tetra(100, BOX_SIZE, N_RAND, SEED, gal_left, ran_left,
                          handedness='left')

    run(f"{gramsci} -gal {gal_left} -ran {ran_left} "
        f"-rmin {RMIN} -rmax {RMAX} -nbins {NBINS} -nmu 1 -out {out_left} -4pcfp")

    data_left = np.loadtxt(out_left, skiprows=1)
    # Chiral tetra sorted bin-tuple: (1, 2, 2, 2, 3, 3)
    row_left = find_row_by_sorted_bins(data_left, [1, 2, 2, 2, 3, 3],
                                        rmin=RMIN, rmax=RMAX, nbins=NBINS)

    # Parity output columns: 12 edges + NNNN RRRR zeta_even NNNN_odd RRRR_odd zeta_odd
    #                         zeta_disc zeta_conn_even zeta_conn_odd
    # Use raw NNNN_odd counts for comparison (not zeta_odd = NNNN_odd/RRRR_odd
    # which is unstable when RRRR_odd ~ 0 for a uniform random catalog)
    nnnn_odd_left = row_left[15]   # NNNN_odd
    nnnn_even_left = row_left[12]  # NNNN (even)
    frac_left = abs(nnnn_odd_left / nnnn_even_left) if nnnn_even_left != 0 else 0
    print(f"  Left-only: NNNN={nnnn_even_left:.6e}, "
          f"NNNN_odd={nnnn_odd_left:.6e}, |odd/even|={frac_left:.6f}")

    assert abs(nnnn_odd_left) > 0, \
        f"Left-only: |NNNN_odd| should be > 0"
    assert frac_left > 0.01, \
        f"Left-only: |odd/even| = {frac_left:.6f} should be > 0.01"
    print("  Left-only parity check passed")

    # --- 3b: Balanced (50 left + 50 right) ---
    generate_chiral_tetra(100, BOX_SIZE, N_RAND, SEED + 1, gal_bal, ran_bal,
                          handedness='balanced')

    run(f"{gramsci} -gal {gal_bal} -ran {ran_bal} "
        f"-rmin {RMIN} -rmax {RMAX} -nbins {NBINS} -nmu 1 -out {out_bal} -4pcfp")

    data_bal = np.loadtxt(out_bal, skiprows=1)
    row_bal = find_row_by_sorted_bins(data_bal, [1, 2, 2, 2, 3, 3],
                                       rmin=RMIN, rmax=RMAX, nbins=NBINS)

    nnnn_odd_bal = row_bal[15]
    nnnn_even_bal = row_bal[12]
    frac_bal = abs(nnnn_odd_bal / nnnn_even_bal) if nnnn_even_bal != 0 else 0
    print(f"  Balanced: NNNN={nnnn_even_bal:.6e}, "
          f"NNNN_odd={nnnn_odd_bal:.6e}, |odd/even|={frac_bal:.6f}")

    # Fractional parity violation should be much smaller when balanced
    assert frac_bal < 0.5 * frac_left, \
        (f"Balanced: |odd/even| = {frac_bal:.6f} should be < "
         f"0.5 * left |odd/even| = {0.5 * frac_left:.6f}")
    print("  Balanced parity cancellation check passed")

    _cleanup(gal_left, ran_left, out_left, gal_bal, ran_bal, out_bal)
    print("  Test PASSED")


def test_analytic_randoms(bindir):
    """Test 4: periodic-box analytic randoms (-box without -ran).

    A uniform catalog in a periodic box is its own null test:
      * the RR column must equal the closed-form shell-volume expectation,
      * NN/RR, NNN/RRR, NNNN/RRRR must all be ~ 1 (xi, zeta, zeta_conn ~ 0),
      * the analytic-mode estimator must still localize a real pair signal.
    """
    print('\n=== Test: Periodic box analytic randoms ===')

    gramsci = os.path.join(bindir, 'gramsci')
    box = 200.0
    n_uni = 6000
    gal_file = 'tmp_uniform.gal'
    generate_uniform_box(n_uni, box, SEED, gal_file)

    base = (f"-rmin {RMIN} -rmax {RMAX} -nbins {NBINS} "
            f"-box {box:g}")

    # --- 4a: 2PCF, exact RR + null xi ---
    out2 = 'tmp_uniform_an.2pcf'
    run(f"{gramsci} -gal {gal_file} {base} -out {out2} -2pcf")
    d2 = np.loadtxt(out2, skiprows=1)
    # Exact analytic RR: (1 - sum w^2) / 2 * Vshell / Vbox with w_i = 1/N
    w2 = (1.0 - 1.0 / n_uni) / 2.0
    edges = np.linspace(RMIN, RMAX, NBINS + 1)
    vshell = 4.0 * np.pi / 3.0 * (edges[1:] ** 3 - edges[:-1] ** 3)
    rr_exact = w2 * vshell / box ** 3
    assert np.allclose(d2[:, 5], rr_exact, rtol=1e-6), \
        f"analytic RR mismatch: {d2[:, 5]} vs {rr_exact}"
    print(f"  2PCF analytic RR matches closed form (rel err "
          f"{np.max(np.abs(d2[:, 5] / rr_exact - 1)):.2e})")
    assert np.all(np.abs(d2[:, 6]) < 0.06), \
        f"uniform box xi should be ~0, got {d2[:, 6]}"
    print(f"  2PCF null test passed: xi = {d2[:, 6]}")

    # --- 4b: 3PCF null ---
    out3 = 'tmp_uniform_an.3pcf'
    run(f"{gramsci} -gal {gal_file} {base} -out {out3} -3pcf")
    d3 = np.loadtxt(out3, skiprows=1)
    ratio3 = d3[:, 6] / d3[:, 7]
    assert abs(np.median(ratio3) - 1.0) < 0.05, \
        f"uniform box NNN/RRR median should be ~1, got {np.median(ratio3):.4f}"
    assert np.all(np.abs(ratio3 - 1.0) < 0.25), \
        f"uniform box NNN/RRR should be ~1 per config, got {ratio3}"
    assert np.median(np.abs(d3[:, 8])) < 0.05, \
        f"uniform box zeta should be ~0, median |zeta| = {np.median(np.abs(d3[:, 8])):.4f}"
    print(f"  3PCF null test passed: median NNN/RRR = {np.median(ratio3):.4f}, "
          f"median |zeta| = {np.median(np.abs(d3[:, 8])):.4f}")

    # --- 4c: 4PCF null + disconnected identity ---
    out4 = 'tmp_uniform_an.4pcf'
    run(f"{gramsci} -gal {gal_file} {base} -out {out4} -4pcf")
    d4 = np.loadtxt(out4, skiprows=1)
    assert d4.shape[1] == 17, f'analytic 4pcf expected 17 columns, got {d4.shape[1]}'
    rrrr = d4[:, 13]
    # restrict to well-populated configs (upper half of analytic RRRR)
    pop = rrrr > np.median(rrrr[rrrr > 0])
    ratio4 = d4[pop, 12] / d4[pop, 13]
    assert abs(np.median(ratio4) - 1.0) < 0.05, \
        f"uniform box NNNN/RRRR median should be ~1, got {np.median(ratio4):.4f}"
    assert np.median(np.abs(d4[pop, 16])) < 0.05, \
        f"uniform box zeta_conn should be ~0, median |zeta_conn| = " \
        f"{np.median(np.abs(d4[pop, 16])):.4f}"
    # identity: zeta_conn = zeta - zeta_disc (same as catalog mode)
    tol = 1e-6 * (np.abs(d4[:, 14]) + np.abs(d4[:, 15])) + 1e-9
    assert np.all(np.abs(d4[:, 16] - (d4[:, 14] - d4[:, 15])) <= tol), \
        'analytic 4pcf disconnected identity failed'
    print(f"  4PCF null test passed: median NNNN/RRRR = {np.median(ratio4):.4f}, "
          f"median |zeta_conn| = {np.median(np.abs(d4[pop, 16])):.4f}")

    # --- 4d: analytic mode still detects a real signal ---
    gal_pairs = 'tmp_pairs_an.gal'
    ran_unused = 'tmp_pairs_an.ran'
    outp = 'tmp_pairs_an.2pcf'
    generate_pairs(100, R0, BOX_SIZE, 10, SEED, gal_pairs, ran_unused)
    run(f"{gramsci} -gal {gal_pairs} -rmin {RMIN} -rmax {RMAX} -nbins {NBINS} "
        f"-box {BOX_SIZE:g} -out {outp} -2pcf")
    dp = np.loadtxt(outp, skiprows=1)
    xi = dp[:, 6]
    print(f"  pairs-in-box xi: {xi}")
    assert xi[1] > 1.0 and xi[1] > xi[0] + 1.0 and xi[1] > xi[2] + 1.0, \
        f"analytic mode should localize the pair signal in bin 2, got {xi}"
    print("  analytic-mode signal localization passed")

    # --- 4e: -box combined with a random catalog (periodic + LS estimator) ---
    ran_file = 'tmp_uniform.ran'
    generate_uniform_box(20000, box, SEED + 99, ran_file)
    outc = 'tmp_uniform_cat.2pcf'
    run(f"{gramsci} -gal {gal_file} -ran {ran_file} {base} -out {outc} -2pcf")
    dc = np.loadtxt(outc, skiprows=1)
    # catalog RR (normalized) should match the analytic values within noise
    ratio = dc[:, 5] / rr_exact / ((1.0 - 1.0 / 20000) / (1.0 - 1.0 / n_uni))
    assert np.all(np.abs(ratio - 1.0) < 0.05), \
        f"periodic catalog RR should match analytic within 5%, got {ratio}"
    print(f"  periodic catalog-RR vs analytic-RR ratios: {ratio}")

    # --- 4f: invalid option combinations must be rejected ---
    r = subprocess.run(f"{gramsci} -gal {gal_file} -rmin {RMIN} -rmax 150 "
                       f"-nbins {NBINS} -box {box:g} -out tmp_bad.out -2pcf",
                       shell=True, capture_output=True, text=True)
    assert 'ERROR' in r.stdout and not os.path.exists('tmp_bad.out'), \
        "-box with rmax >= L/2 should be rejected"
    r = subprocess.run(f"{gramsci} -gal {gal_file} -rmin {RMIN} -rmax 60 "
                       f"-nbins {NBINS} -box {box:g} -out tmp_bad.out -3pcf",
                       shell=True, capture_output=True, text=True)
    assert 'ERROR' in r.stdout and not os.path.exists('tmp_bad.out'), \
        "-box 3PCF with rmax > L/4 should be rejected"
    r = subprocess.run(f"{gramsci} -gal {gal_file} -rmin {RMIN} -rmax {RMAX} "
                       f"-nbins {NBINS} -nmu 4 -box {box:g} -out tmp_bad.out -2pcf",
                       shell=True, capture_output=True, text=True)
    assert 'ERROR' in r.stdout and not os.path.exists('tmp_bad.out'), \
        "-box with -nmu > 1 should be rejected"
    print("  invalid-combination rejection passed")

    _cleanup(gal_file, out2, out3, out4, gal_pairs, ran_unused, outp,
             ran_file, outc)
    print("  Test PASSED")


def test_combined_modes(bindir):
    """Test 6: combined query modes share one graph build.

    Running -2pcf -3pcf -equi -4pcf together must produce per-mode output
    files (<out>.<mode>) that are byte-identical to the corresponding
    single-mode runs, and the single-mode runs must keep the exact -out name.
    """
    print('\n=== Test: Combined query modes ===')

    gramsci = os.path.join(bindir, 'gramsci')
    box = 200.0
    gal_file = 'tmp_uniform.gal'
    generate_uniform_box(6000, box, SEED, gal_file)
    base = f"-gal {gal_file} -rmin {RMIN} -rmax {RMAX} -nbins {NBINS} -box {box:g}"

    modes = ['2pcf', '3pcf', 'equi', '4pcf']
    for m in modes:
        run(f"{gramsci} {base} -out tmp_single.{m} -{m}")
    flags = ' '.join(f'-{m}' for m in modes)
    run(f"{gramsci} {base} -out tmp_combined {flags}")

    for m in modes:
        combined = f'tmp_combined.{m}'
        single = f'tmp_single.{m}'
        assert os.path.exists(combined), f"combined run did not write {combined}"
        with open(combined) as fc, open(single) as fs:
            assert fc.read() == fs.read(), \
                f"combined-mode {m} output differs from the single-mode run"
        print(f"  {m}: combined output identical to single-mode run")

    # a combined run must not also write the bare -out name
    assert not os.path.exists('tmp_combined'), \
        "combined run should only write suffixed outputs"

    _cleanup(gal_file, *[f'tmp_single.{m}' for m in modes],
             *[f'tmp_combined.{m}' for m in modes])
    print("  Test PASSED")


def test_aniso_disconnected(bindir):
    """Test 5: anisotropic (RSD-aware) disconnected-4PCF subtraction.

    A z-squashed blob field in a periodic box has a strong quadrupole
    xi_2(r).  The zeta_disc column must equal the orientation-averaged
    Gaussian term
        sum_pairings [ xi0*xi0 + xi2*xi2*L2(ct)/5 + xi4*xi4*L4(ct)/9 ]
    with xi_ell measured here independently by direct minimum-image pair
    counting and ct the opposite-edge angle cosine from the bin centres.
    """
    print('\n=== Test: Anisotropic disconnected 4PCF subtraction ===')

    gramsci = os.path.join(bindir, 'gramsci')
    box = 200.0
    gal_file = 'tmp_aniso.gal'
    out4 = 'tmp_aniso.4pcf'
    generate_aniso_box(300, 20, box, 8.0, 3.0, SEED + 7, gal_file)
    run(f"{gramsci} -gal {gal_file} -rmin {RMIN} -rmax {RMAX} -nbins {NBINS} "
        f"-box {box:g} -out {out4} -4pcf")

    # Independent xi_ell by direct min-image pair counting (each pair once)
    pts = np.loadtxt(gal_file)[:, :3]
    n = len(pts)
    edges = np.linspace(RMIN, RMAX, NBINS + 1)
    s0 = np.zeros(NBINS)
    s2 = np.zeros(NBINS)
    s4 = np.zeros(NBINS)
    for i in range(0, n, 500):
        d = pts[i:i + 500, None, :] - pts[None, :, :]
        d -= box * np.round(d / box)
        r = np.sqrt((d ** 2).sum(-1))
        once = np.arange(i, min(i + 500, n))[:, None] < np.arange(n)[None, :]
        sel = once & (r > RMIN) & (r < RMAX)
        mu2 = (d[..., 2][sel] / r[sel]) ** 2
        b = np.minimum(NBINS - 1,
                       ((r[sel] - RMIN) / (RMAX - RMIN) * NBINS).astype(int))
        np.add.at(s0, b, 1.0)
        np.add.at(s2, b, 1.5 * mu2 - 0.5)
        np.add.at(s4, b, 4.375 * mu2 ** 2 - 3.75 * mu2 + 0.375)
    rr = (n ** 2 - n) / 2 * (4 * np.pi / 3) \
        * (edges[1:] ** 3 - edges[:-1] ** 3) / box ** 3
    xi0 = s0 / rr - 1
    xi2 = 5 * s2 / rr
    xi4 = 9 * s4 / rr
    assert np.max(np.abs(xi2)) > 0.3, \
        f"test field should have a strong quadrupole, got xi2 = {xi2}"

    d4 = np.loadtxt(out4, skiprows=1)
    rc_all = 0.5 * (edges[:-1] + edges[1:])
    dr = (RMAX - RMIN) / NBINS
    bidx = np.rint((d4[:, 0:12:2] - RMIN) / dr).astype(int)  # 0-based b1..b6

    def leg2(x):
        return 1.5 * x * x - 0.5

    def leg4(x):
        return 4.375 * x ** 4 - 3.75 * x ** 2 + 0.375

    disc_py = np.zeros(len(d4))
    disc_iso = np.zeros(len(d4))
    for row in range(len(d4)):
        b = bidx[row]
        rc = rc_all[b]
        rc2 = rc ** 2
        cts = [(rc2[2] + rc2[3] - rc2[1] - rc2[4]) / (2 * rc[0] * rc[5]),
               (rc2[2] + rc2[3] - rc2[0] - rc2[5]) / (2 * rc[1] * rc[4]),
               (rc2[1] + rc2[4] - rc2[0] - rc2[5]) / (2 * rc[2] * rc[3])]
        for (ea, eb, ct) in [(0, 5, cts[0]), (1, 4, cts[1]), (2, 3, cts[2])]:
            ct = min(1.0, max(-1.0, ct))
            disc_py[row] += (xi0[b[ea]] * xi0[b[eb]]
                             + xi2[b[ea]] * xi2[b[eb]] * leg2(ct) / 5
                             + xi4[b[ea]] * xi4[b[eb]] * leg4(ct) / 9)
            disc_iso[row] += xi0[b[ea]] * xi0[b[eb]]

    disc_f = d4[:, 15]
    err = np.abs(disc_py - disc_f) / np.maximum(np.abs(disc_py), 1e-3)
    assert np.max(err) < 1e-4, \
        f"zeta_disc mismatch vs independent computation: max rel {np.max(err):.2e}"
    print(f"  zeta_disc matches independent computation "
          f"(max rel diff {np.max(err):.2e}, {len(d4)} configs)")
    corr = np.median(np.abs(disc_py - disc_iso))
    assert corr > 0.01, \
        "anisotropic correction should be active on this field"
    print(f"  anisotropic correction active: median |disc - disc_iso| = {corr:.4f}")

    _cleanup(gal_file, out4)
    print("  Test PASSED")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    bindir = os.path.join('..', 'bin')
    gal = os.path.join('..', 'example', 'test.gal')
    ran = os.path.join('..', 'example', 'test.ran')

    # ---- Existing regression tests ----

    out2 = 'tmp_2pcf.out'
    cmd2 = f"{os.path.join(bindir, 'gramsci')} -gal {gal} -ran {ran} -rmin 1.0 -rmax 30.0 -nbins 10 -nmu 10 -out {out2} -2pcf"
    run(cmd2)

    tmp=np.loadtxt(out2,skiprows=1)

    DD=np.mean(tmp[:,4])
    RR=np.mean(tmp[:,5])

    assert_close(DD, 1.062171939e-07, 1e-5, '2pcf DD')
    assert_close(RR, 1.681313979e-07, 1e-5, '2pcf RR')

    out3 = 'tmp_3pcf.out'
    cmd3 = f"{os.path.join(bindir, 'gramsci')} -gal {gal} -ran {ran} -rmin 1.0 -rmax 30.0 -nbins 6 -nmu 1 -out {out3} -3pcf"
    run(cmd3)
    tmp=np.loadtxt(out3,skiprows=1)

    DDD=np.mean(tmp[:,6])
    RRR=np.mean(tmp[:,7])

    assert_close(DDD, 1.986803931535715e-12, 1e-5, '3pcf DDD')
    assert_close(RRR, 1.683093412466072e-12, 1e-5, '3pcf RRR')

    # 4PCF non-parity smoke test (small nbins to keep bintable6 tiny)
    out4 = 'tmp_4pcf.out'
    cmd4 = f"{os.path.join(bindir, 'gramsci')} -gal {gal} -ran {ran} -rmin 1.0 -rmax 30.0 -nbins 3 -nmu 1 -out {out4} -4pcf"
    run(cmd4)
    with open(out4) as f:
        header = f.readline()
        assert 'NNNN' in header and 'RRRR' in header and 'zeta' in header, \
            f'4pcf header check failed: {header}'
        assert 'zeta_disc' in header and 'zeta_conn' in header, \
            f'4pcf disconnected columns missing: {header}'
        lines = f.readlines()
        assert len(lines) > 0, '4pcf produced no output rows'

    # Verify disconnected subtraction: zeta_conn = (zeta - 1) - zeta_disc
    tmp4 = np.loadtxt(out4, skiprows=1)
    assert tmp4.shape[1] == 17, f'4pcf expected 17 columns, got {tmp4.shape[1]}'
    zeta_col = tmp4[:, 14]      # zeta = N4/R4
    zeta_disc = tmp4[:, 15]     # zeta_disc
    zeta_conn = tmp4[:, 16]     # zeta_conn
    # Check zeta_conn = zeta - zeta_disc for rows where zeta != 0
    # (zeta = N4/R4 is the total 4PCF under signed weights; no "-1").
    # Tolerance scales with the inputs because the columns are written at
    # e14.7 precision, so reconstructing the identity from the rounded zeta
    # and zeta_disc columns carries ~1e-7*(|zeta|+|zeta_disc|) rounding.
    mask = zeta_col != 0.0
    if np.any(mask):
        expected_conn = zeta_col[mask] - zeta_disc[mask]
        tol = 1e-6 * (np.abs(zeta_col[mask]) + np.abs(zeta_disc[mask])) + 1e-9
        assert np.all(np.abs(zeta_conn[mask] - expected_conn) <= tol), \
            f'4pcf disconnected identity failed: zeta_conn != zeta - zeta_disc'

    # 4PCF parity smoke test
    out4p = 'tmp_4pcfp.out'
    cmd4p = f"{os.path.join(bindir, 'gramsci')} -gal {gal} -ran {ran} -rmin 1.0 -rmax 30.0 -nbins 3 -nmu 1 -out {out4p} -4pcfp"
    run(cmd4p)
    with open(out4p) as f:
        header = f.readline()
        assert 'zeta_even' in header and 'zeta_odd' in header, \
            f'4pcfp header check failed: {header}'
        assert 'zeta_disc' in header and 'zeta_conn_even' in header and 'zeta_conn_odd' in header, \
            f'4pcfp disconnected columns missing: {header}'
        lines = f.readlines()
        assert len(lines) > 0, '4pcfp produced no output rows'

    # Verify parity disconnected subtraction
    tmp4p = np.loadtxt(out4p, skiprows=1)
    assert tmp4p.shape[1] == 21, f'4pcfp expected 21 columns, got {tmp4p.shape[1]}'
    zeta_even = tmp4p[:, 14]       # zeta_even
    zeta_odd = tmp4p[:, 17]        # zeta_odd
    zeta_disc_p = tmp4p[:, 18]     # zeta_disc
    zeta_conn_even = tmp4p[:, 19]  # zeta_conn_even
    zeta_conn_odd = tmp4p[:, 20]   # zeta_conn_odd
    # Check even: zeta_conn_even = zeta_even - zeta_disc (no "-1");
    # precision-scaled tolerance for the e14.7-rounded columns (see 4pcf above).
    mask = zeta_even != 0.0
    if np.any(mask):
        expected_even = zeta_even[mask] - zeta_disc_p[mask]
        tol = 1e-6 * (np.abs(zeta_even[mask]) + np.abs(zeta_disc_p[mask])) + 1e-9
        assert np.all(np.abs(zeta_conn_even[mask] - expected_even) <= tol), \
            f'4pcfp even disconnected identity failed'
    # Check odd: zeta_conn_odd = zeta_odd (disconnected is zero for odd)
    assert np.allclose(zeta_conn_odd, zeta_odd, atol=1e-20), \
        f'4pcfp odd channel: zeta_conn_odd != zeta_odd'

    os.remove(out2)
    os.remove(out3)
    os.remove(out4)
    os.remove(out4p)

    print('Existing correlation tests passed')

    # ---- Physics regression tests with synthetic catalogs ----

    test_pairs_2pcf(bindir)
    test_regular_tetra_connected(bindir)
    test_chiral_parity(bindir)
    test_analytic_randoms(bindir)
    test_combined_modes(bindir)
    test_aniso_disconnected(bindir)

    print('\n========================================')
    print('All correlation tests passed')
    print('========================================')

if __name__ == '__main__':
    main()
