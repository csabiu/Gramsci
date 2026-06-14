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

    with open(os.path.join(OUT, 'ezmock_3pcf_allconfig.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['config', 'desi', 'mock_mean', 'mock_sigma', 'valid'])
        for i in range(len(mean)):
            w.writerow([i + 1, f'{zdesi[i]:.8e}',
                        f'{mean[i]:.8e}' if np.isfinite(mean[i]) else 'nan',
                        f'{sig[i]:.8e}', int(good[i])])

    r1c = 0.5 * (ec[:, 0] + ec[:, 1])
    r2c = 0.5 * (ec[:, 2] + ec[:, 3])
    r3c = 0.5 * (ec[:, 4] + ec[:, 5])
    with open(os.path.join(OUT, 'threepcf_config_scales.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['config', 'r1', 'r2', 'r3'])
        for i in range(len(r1c)):
            w.writerow([i + 1, f'{r1c[i]:.3f}', f'{r2c[i]:.3f}', f'{r3c[i]:.3f}'])
    print(f'  ezmock_3pcf_allconfig.csv + threepcf_config_scales.csv  '
          f'(nmock={nm}, nconfig={len(mean)})')


# ---------------------------------------------------------------------------
# Fig 8: DESI 4PCF vs EZmock, all 276 configs (S^4 weights, validity mask)
# ---------------------------------------------------------------------------
def reduce_4pcf():
    import compare_desi_ezmock_4pcf as c4
    msums = c4.load_mock_sums()
    mock_ids = sorted({int(re.search(r'mock(\d+)', f).group(1))
                       for f in glob.glob(os.path.join(EZ, 'NGC_zbin1_mock*.4pcf'))})
    rows = [['zbin', 'config', 'desi', 'mock_mean', 'mock_sigma', 'valid']]
    for zbin in (1, 2):
        stack = []
        for n in mock_ids:
            fn = [os.path.join(EZ, f'{c}_zbin{zbin}_mock{n}.4pcf') for c in ('NGC', 'SGC')]
            if not all(os.path.exists(x) for x in fn):
                continue
            try:
                wN = msums[('NGC', n, zbin)]**4
                wS = msums[('SGC', n, zbin)]**4
            except KeyError:
                continue
            stack.append(c4.combine(c4.load(fn[0]), c4.load(fn[1]), wN, wS))
        stack = np.array(stack)
        good = np.all(np.isfinite(stack) & (np.abs(stack) < 1.0), axis=0)
        stack = np.where(good[None, :], stack, np.nan)
        mean = np.nanmean(stack, axis=0)
        sig = np.nanstd(stack, axis=0, ddof=1)
        sig[~good] = 0.0
        ngc = c4.load(os.path.join(DESI, f'LRG_NGC_zbin{zbin}_r65.4pcf'))
        sgc = c4.load(os.path.join(DESI, f'LRG_SGC_zbin{zbin}_r65.4pcf'))
        zd = c4.combine(ngc, sgc, c4.S_DESI[('NGC', zbin)]**4, c4.S_DESI[('SGC', zbin)]**4)
        for i in range(len(mean)):
            rows.append([zbin, i + 1, f'{zd[i]:.8e}',
                         f'{mean[i]:.8e}' if np.isfinite(mean[i]) else 'nan',
                         f'{sig[i]:.8e}', int(good[i])])
    with open(os.path.join(OUT, 'ezmock_4pcf_band.csv'), 'w', newline='') as f:
        csv.writer(f).writerows(rows)
    print('  ezmock_4pcf_band.csv')


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
