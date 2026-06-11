import os
import subprocess
import math
import numpy as np
from generate_test_fields import (generate_pairs, generate_regular_tetra,
                                  generate_chiral_tetra,
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
    # Check zeta_conn = (zeta - 1) - zeta_disc for rows where zeta != 0
    mask = zeta_col != 0.0
    if np.any(mask):
        expected_conn = (zeta_col[mask] - 1.0) - zeta_disc[mask]
        assert np.allclose(zeta_conn[mask], expected_conn, atol=1e-20), \
            f'4pcf disconnected identity failed: zeta_conn != (zeta-1) - zeta_disc'

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
    # Check even: zeta_conn_even = (zeta_even - 1) - zeta_disc
    mask = zeta_even != 0.0
    if np.any(mask):
        expected_even = (zeta_even[mask] - 1.0) - zeta_disc_p[mask]
        assert np.allclose(zeta_conn_even[mask], expected_even, atol=1e-20), \
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

    print('\n========================================')
    print('All correlation tests passed')
    print('========================================')

if __name__ == '__main__':
    main()
