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
# Fig 7: DESI 3PCF vs EZmock, isoceles slice (10 Mpc rebin)
# ---------------------------------------------------------------------------
def reduce_3pcf():
    from mock_errorbars_3pcf import (load_sums, combine_mock_counts,
                                     rebin_counts, slice_cut, DATA_COMB)
    sums = load_sums()
    mock_ids = sorted({int(re.search(r'mock(\d+)', f).group(1))
                       for f in glob.glob(os.path.join(EZ, 'NGC_zbin1_mock*.3pcf'))})
    rows = [['zbin', 'r3', 'desi', 'mock_mean', 'mock_sigma']]
    nmock = {}
    for zbin in (1, 2):
        slices = []
        for n in mock_ids:
            res = combine_mock_counts(n, zbin, sums)
            if res is None:
                continue
            e6, nnn, rrr = res
            ec, zc = rebin_counts(e6, nnn, rrr)
            r_s, z_s, _ = slice_cut(ec, zc)
            slices.append(z_s)
        slices = np.array(slices)
        nmock[zbin] = slices.shape[0]
        mmean = np.nanmean(slices, axis=0)
        msig = np.nanstd(slices, axis=0, ddof=1)
        data = np.loadtxt(DATA_COMB[zbin])
        e6d, zd = rebin_counts(data[:, :6], data[:, 6], data[:, 7])
        r_d, z_d, _ = slice_cut(e6d, zd)
        for r, d, m, s in zip(r_d, z_d, mmean, msig):
            rows.append([zbin, f'{r:.4f}', f'{d:.8e}', f'{m:.8e}', f'{s:.8e}'])
    with open(os.path.join(OUT, 'ezmock_3pcf_slice.csv'), 'w', newline='') as f:
        csv.writer(f).writerows(rows)
    print(f'  ezmock_3pcf_slice.csv  (nmock={nmock})')


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
