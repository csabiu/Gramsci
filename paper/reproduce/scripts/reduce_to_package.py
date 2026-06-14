"""Reduce the full per-realization measurement ensembles to the compact CSV
products shipped in ../data/.  Run ONCE on the machine that holds the raw
measurement outputs (~/data/desi/ezmocks, ~/data/quijote/results); the
shipped figure scripts then read only the CSVs.

This documents the raw -> reduced step transparently.  The reductions
replicate exactly the published comparison scripts:
  - Fig 7: compare_desi_ezmock_3pcf.py   (isoceles slice, 10 Mpc rebin)
  - Fig 8: compare_desi_ezmock_4pcf.py   (S^4 cap weights, validity mask)
  - Fig 9: analyze_responses.py          (equilateral + isoceles cuts)
  - parity null: chi2/dof of fiducial parity-odd 4PCF vs zero
"""
import csv
import glob
import os
import re
import sys

import numpy as np

EZ = '/home/csabiu/data/desi/ezmocks'
DESI = '/home/csabiu/data/desi'
QJ = '/home/csabiu/data/quijote/results'
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')
sys.path.insert(0, EZ)


def loadtxt_skip(path):
    n = 0
    with open(path) as f:
        for line in f:
            s = line.strip()
            if s and (s[0].isdigit() or s[0] in '+-'):
                break
            n += 1
    return np.loadtxt(path, skiprows=n)


# ---------------------------------------------------------------------------
# Fig 7: DESI 3PCF vs EZmock, ALL configurations, NGC+SGC and both redshift
# bins combined at the count level (S^3 effective-volume weighting), 10 Mpc
# rebin.  Also the three sorted triangle-side scales per configuration.
# ---------------------------------------------------------------------------
S_DESI_3PCF = {('NGC', 1): 1.798173e5, ('NGC', 2): 1.754226e5,
               ('SGC', 1): 9.207662e4, ('SGC', 2): 8.807297e4}


def reduce_3pcf():
    from mock_errorbars_3pcf import load_sums, combine_mock_counts, rebin_counts
    sums = load_sums()
    mock_ids = sorted({int(re.search(r'mock(\d+)', f).group(1))
                       for f in glob.glob(os.path.join(EZ, 'NGC_zbin1_mock*.3pcf'))})

    # DESI: sum S^3-weighted counts over the four (cap, zbin) sub-samples
    nnn = rrr = None
    e6 = None
    for cap in ('NGC', 'SGC'):
        for z in (1, 2):
            d = loadtxt_skip(os.path.join(DESI, f'LRG_{cap}_zbin{z}.3pcf'))
            e6 = d[:, :6]
            w = S_DESI_3PCF[(cap, z)]**3
            nnn = w * d[:, 6] if nnn is None else nnn + w * d[:, 6]
            rrr = w * d[:, 7] if rrr is None else rrr + w * d[:, 7]
    ec, zdesi = rebin_counts(e6, nnn, rrr)

    # mocks: per mock, add the two zbins' cap-combined counts, then rebin
    stack = []
    for n in mock_ids:
        r1 = combine_mock_counts(n, 1, sums)
        r2 = combine_mock_counts(n, 2, sums)
        if r1 is None or r2 is None:
            continue
        e, n1, rr1 = r1
        _, n2, rr2 = r2
        _, zc = rebin_counts(e, n1 + n2, rr1 + rr2)
        stack.append(zc)
    stack = np.array(stack)
    nm = stack.shape[0]
    good = np.all(np.isfinite(stack) & (np.abs(stack) < 1e3), axis=0)
    safe = np.where(good, stack, np.nan)
    mean = np.where(good, np.nanmean(safe, axis=0), np.nan)
    sig = np.where(good, np.nanstd(safe, axis=0, ddof=1), 0.0)

    # Keep only configurations that form a genuine triangle: the bin-centre
    # side lengths must satisfy the triangle inequality r1 + r2 > r3 (with
    # r3 the largest, since the canonical tuple is sorted).  The enumeration
    # includes geometrically impossible combinations such as (55, 55, 145),
    # which have zero counts; we drop them so every plotted configuration is
    # a real triangle.  Configurations are renumbered 1..N over the kept set.
    r1c = 0.5 * (ec[:, 0] + ec[:, 1])
    r2c = 0.5 * (ec[:, 2] + ec[:, 3])
    r3c = 0.5 * (ec[:, 4] + ec[:, 5])
    tri = (r1c + r2c) > r3c
    ndrop = int((~tri).sum())

    with open(os.path.join(OUT, 'ezmock_3pcf_allconfig.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['config', 'desi', 'mock_mean', 'mock_sigma', 'valid'])
        c = 0
        for i in np.where(tri)[0]:
            c += 1
            w.writerow([c, f'{zdesi[i]:.8e}',
                        f'{mean[i]:.8e}' if np.isfinite(mean[i]) else 'nan',
                        f'{sig[i]:.8e}', int(good[i])])

    with open(os.path.join(OUT, 'threepcf_config_scales.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['config', 'r1', 'r2', 'r3'])
        c = 0
        for i in np.where(tri)[0]:
            c += 1
            w.writerow([c, f'{r1c[i]:.3f}', f'{r2c[i]:.3f}', f'{r3c[i]:.3f}'])
    print(f'  ezmock_3pcf_allconfig.csv + threepcf_config_scales.csv  '
          f'(nmock={nm}, nconfig={int(tri.sum())}, dropped {ndrop} non-triangles)')


# ---------------------------------------------------------------------------
# Fig 8: DESI CONNECTED 4PCF vs EZmock, per redshift bin (NGC+SGC combined).
#
# zeta_conn = zeta_total - zeta_disc, where
#   zeta_total = (S_N^4 N4_N + S_S^4 N4_S)/(S_N^4 R4_N + S_S^4 R4_S)   (cols 12,13)
#   zeta_disc(config) = xi(b1)xi(b6) + xi(b2)xi(b5) + xi(b3)xi(b4)
# with xi the (count-level S^2-combined) internal 2PCF.  The internal 2PCF is
# read from the run logs and corrected for the historical "-1" estimator bug
# (xi_correct = xi_logged + 1; the bug was exactly the spurious -1, now fixed
# in the Fortran).  xi_NGC ~ xi_SGC to a few percent, so the S^2-weighted cap
# average of xi is an excellent stand-in for the count-level combined 2PCF
# (the residual error in zeta_disc is second order in xi_N - xi_S).
# This reconstruction was validated bin-for-bin against the fixed binary's
# zeta_conn column.
# ---------------------------------------------------------------------------
S_DESI_4PCF = {('NGC', 1): 1.798173e5, ('NGC', 2): 1.754226e5,
               ('SGC', 1): 9.207662e4, ('SGC', 2): 8.807297e4}
RMIN4, RMAX4, NB4 = 20.0, 65.0, 4


def _parse_internal_xi(logpath):
    """Return the corrected internal 2PCF (xi_logged + 1) from a 4pcf run log."""
    xi = []
    with open(logpath) as f:
        block = False
        for line in f:
            if 'Internal 2PCF computed' in line:
                block = True
                continue
            if block:
                m = re.search(r'xi =\s*(\S+)', line)
                if m:
                    xi.append(float(m.group(1)) + 1.0)   # undo the historical -1
                elif xi:
                    break
    return np.array(xi)


SUBSAMPLES = [('NGC', 1), ('NGC', 2), ('SGC', 1), ('SGC', 2)]


def _tetrahedron_mask():
    """True for configurations whose six bin-centre edge lengths form a real
    tetrahedron (Cayley-Menger determinant > 0).  The canonical enumeration
    includes 6-tuples that satisfy every face triangle inequality yet are not
    embeddable in 3D (e.g. five edges of 26 and one of 48); these are dropped
    so every plotted configuration is a genuine tetrahedron, mirroring the
    triangle-inequality cut applied to the 3PCF.  Geometry is data-independent.
    """
    d = loadtxt_skip(os.path.join(DESI, 'LRG_NGC_zbin1_r65.4pcf'))
    c = 0.5 * (d[:, 0:12:2] + d[:, 1:12:2])     # d12,d13,d14,d23,d24,d34
    mask = np.zeros(len(c), bool)
    for i, e in enumerate(c):
        d12, d13, d14, d23, d24, d34 = e**2
        M = np.array([[0, 1, 1, 1, 1],
                      [1, 0, d12, d13, d14],
                      [1, d12, 0, d23, d24],
                      [1, d13, d23, 0, d34],
                      [1, d14, d24, d34, 0]], float)
        mask[i] = np.linalg.det(M) > 0
    return mask


def _full_and_conn(file_for, log_for, S_for):
    """Full and connected 4PCF for ALL sub-samples combined (NGC+SGC, both
    redshift bins) at the count level.  Returns (ztot, zconn) arrays over
    configs, or (None, None) if inputs are missing.

    ztot = (sum S^4 N4) / (sum S^4 R4)                 [count-level total]
    zconn = ztot - zdisc,  zdisc from the S^2-combined internal 2PCF.
    """
    files = [file_for(c, z) for c, z in SUBSAMPLES]
    logs = [log_for(c, z) for c, z in SUBSAMPLES]
    if not all(os.path.exists(p) for p in files + logs):
        return None, None
    n4 = r4 = None
    xi_num = xi_den = None
    edges = None
    for (c, z), fp, lp in zip(SUBSAMPLES, files, logs):
        d = loadtxt_skip(fp)
        edges = d[:, 0:12:2]
        S = S_for(c, z)
        n4 = S**4 * d[:, 12] if n4 is None else n4 + S**4 * d[:, 12]
        r4 = S**4 * d[:, 13] if r4 is None else r4 + S**4 * d[:, 13]
        xi = _parse_internal_xi(lp)
        if xi.size != NB4:
            return None, None
        xi_num = S**2 * xi if xi_num is None else xi_num + S**2 * xi
        xi_den = S**2 if xi_den is None else xi_den + S**2
    with np.errstate(divide='ignore', invalid='ignore'):
        ztot = n4 / r4
    xi_all = xi_num / xi_den
    b = np.clip(np.rint((edges - RMIN4) / ((RMAX4 - RMIN4) / NB4)).astype(int),
                0, NB4 - 1)
    zdisc = (xi_all[b[:, 0]] * xi_all[b[:, 5]] + xi_all[b[:, 1]] * xi_all[b[:, 4]]
             + xi_all[b[:, 2]] * xi_all[b[:, 3]])
    return ztot, ztot - zdisc


def reduce_4pcf():
    mock_ids = sorted({int(re.search(r'mock(\d+)', f).group(1))
                       for f in glob.glob(os.path.join(EZ, 'NGC_zbin1_mock*.4pcf'))})

    # DESI: full and connected, all four sub-samples combined
    desi_full, desi_conn = _full_and_conn(
        lambda c, z: os.path.join(DESI, f'LRG_{c}_zbin{z}_r65.4pcf'),
        lambda c, z: os.path.join(DESI, f'LRG_{c}_zbin{z}_r65.4pcf.log'),
        lambda c, z: S_DESI_4PCF[(c, z)])

    # mock ensemble
    full_stack, conn_stack = [], []
    for n in mock_ids:
        ft, cn = _full_and_conn(
            lambda c, z, n=n: os.path.join(EZ, f'{c}_zbin{z}_mock{n}.4pcf'),
            lambda c, z, n=n: os.path.join(EZ, f'{c}_zbin{z}_mock{n}.4pcf.log'),
            lambda c, z, n=n: _msum(c, n, z))
        if ft is not None:
            full_stack.append(ft); conn_stack.append(cn)
    full_stack = np.array(full_stack)
    conn_stack = np.array(conn_stack)
    nm = full_stack.shape[0]

    def stats(stack, desi):
        good = np.all(np.isfinite(stack) & (np.abs(stack) < 1e3), axis=0) \
            & np.isfinite(desi) & (np.abs(desi) < 1e3)
        safe = np.where(good[None, :], stack, np.nan)
        return (np.where(good, np.nanmean(safe, axis=0), np.nan),
                np.where(good, np.nanstd(safe, axis=0, ddof=1), 0.0), good)

    fmean, fsig, fgood = stats(full_stack, desi_full)
    cmean, csig, cgood = stats(conn_stack, desi_conn)
    # Keep only genuine tetrahedra (Cayley-Menger) that are also well-measured;
    # renumber 1..N over the kept set so the configuration axis has no gaps.
    keep = _tetrahedron_mask() & fgood & cgood
    ndrop = int((~keep).sum())

    rows = [['config', 'desi_full', 'full_mean', 'full_sig',
             'desi_conn', 'conn_mean', 'conn_sig', 'valid']]
    c = 0
    for i in np.where(keep)[0]:
        c += 1
        rows.append([c,
                     f'{desi_full[i]:.8e}',
                     f'{fmean[i]:.8e}' if np.isfinite(fmean[i]) else 'nan',
                     f'{fsig[i]:.8e}',
                     f'{desi_conn[i]:.8e}',
                     f'{cmean[i]:.8e}' if np.isfinite(cmean[i]) else 'nan',
                     f'{csig[i]:.8e}', 1])
    with open(os.path.join(OUT, 'ezmock_4pcf_allcomb.csv'), 'w', newline='') as f:
        csv.writer(f).writerows(rows)
    print(f'  ezmock_4pcf_allcomb.csv (full + connected, all combined, '
          f'{int(keep.sum())} tetrahedra, dropped {ndrop} non-tetrahedra/'
          f'invalid; {nm} mocks)')


# mock weight sums for the 4PCF (zbin_sums_4pcf.txt)
_MSUMS = None


def _msum(cap, n, zbin):
    global _MSUMS
    if _MSUMS is None:
        _MSUMS = {}
        with open(os.path.join(EZ, 'zbin_sums_4pcf.txt')) as f:
            for line in f:
                label, zb, ngal, sw = line.split()
                c, mk = label.split('_mock')
                _MSUMS[(c, int(mk), int(zb[-1]))] = float(sw)
    return _MSUMS[(cap, n, zbin)]


# ---------------------------------------------------------------------------
# Fig 9: Quijote fiducial 3PCF, equilateral + isoceles(48,76) cuts
# ---------------------------------------------------------------------------
def reduce_quijote():
    files = sorted(glob.glob(os.path.join(QJ, 'fiducial_r*.3pcf')))
    stack = np.array([loadtxt_skip(f)[:, 8] for f in files])   # zeta
    edges = loadtxt_skip(files[0])[:, :6]
    r1 = 0.5 * (edges[:, 0] + edges[:, 1])
    r2 = 0.5 * (edges[:, 2] + edges[:, 3])
    r3 = 0.5 * (edges[:, 4] + edges[:, 5])
    mean = np.nanmean(stack, axis=0)
    sig = np.nanstd(stack, axis=0, ddof=1)

    # equilateral
    eq = (np.abs(r1 - r2) < 1e-6) & (np.abs(r2 - r3) < 1e-6)
    o = np.argsort(r1[eq])
    rows = [['cut', 'r', 'zeta_mean', 'zeta_sigma']]
    for r, m, s in zip(r1[eq][o], mean[eq][o], sig[eq][o]):
        rows.append(['equilateral', f'{r:.4f}', f'{m:.8e}', f'{s:.8e}'])

    # isoceles slice (48,76)
    centers = np.unique(np.concatenate([r1, r2, r3]))
    t1 = centers[np.argmin(np.abs(centers - 48.0))]
    t2 = centers[np.argmin(np.abs(centers - 76.0))]
    r_s, m_s, s_s = [], [], []
    for row in range(edges.shape[0]):
        trip = [r1[row], r2[row], r3[row]]
        used = [False] * 3
        ok = True
        for t in (t1, t2):
            for j in range(3):
                if not used[j] and abs(trip[j] - t) < 1e-6:
                    used[j] = True
                    break
            else:
                ok = False
                break
        if ok:
            j = used.index(False)
            r_s.append(trip[j]); m_s.append(mean[row]); s_s.append(sig[row])
    o = np.argsort(r_s)
    for i in o:
        rows.append([f'isoceles_{t1:.0f}_{t2:.0f}', f'{r_s[i]:.4f}',
                     f'{m_s[i]:.8e}', f'{s_s[i]:.8e}'])
    with open(os.path.join(OUT, 'quijote_fiducial_3pcf.csv'), 'w', newline='') as f:
        csv.writer(f).writerows(rows)
    print(f'  quijote_fiducial_3pcf.csv  (nreal={stack.shape[0]})')


# ---------------------------------------------------------------------------
# Parity null: chi2/dof of fiducial parity-odd 4PCF vs zero (text claim)
# ---------------------------------------------------------------------------
def reduce_parity_null():
    files = sorted(glob.glob(os.path.join(QJ, 'fiducial_r*.4pcfp')))
    stack = np.array([loadtxt_skip(f)[:, 15] for f in files])  # NNNN_odd
    nm = stack.shape[0]
    mean = np.nanmean(stack, axis=0)
    err = np.nanstd(stack, axis=0, ddof=1) / np.sqrt(nm)
    m = np.isfinite(mean) & np.isfinite(err) & (err > 0)
    z = mean[m] / err[m]
    chi2, ndof = float(np.sum(z**2)), int(m.sum())
    t_exp = (nm - 1) / (nm - 3)
    with open(os.path.join(OUT, 'parity_null_summary.txt'), 'w') as f:
        f.write(f'# Quijote fiducial parity-odd 4PCF null test (raw NNNN_odd vs 0)\n')
        f.write(f'n_realizations {nm}\n')
        f.write(f'n_configs {ndof}\n')
        f.write(f'chi2 {chi2:.3f}\n')
        f.write(f'chi2_per_dof {chi2/ndof:.4f}\n')
        f.write(f'student_t_expectation {t_exp:.4f}\n')
        f.write(f'frac_beyond_2sigma {(np.abs(z) > 2).mean():.4f}\n')
        f.write(f'max_abs_z {np.abs(z).max():.3f}\n')
    print(f'  parity_null_summary.txt  (chi2/dof={chi2/ndof:.3f}, nm={nm})')


if __name__ == '__main__':
    print('Reducing ensembles to ../data/ CSVs:')
    reduce_3pcf()
    reduce_4pcf()
    reduce_quijote()
    reduce_parity_null()
    print('done.')
