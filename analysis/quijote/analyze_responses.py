"""Stack the Quijote campaign measurements and build the paper figures:

1. Fiducial 3PCF with error bars (validation figure) — equilateral cut and
   an isoceles slice.
2. Fractional response of the 3PCF to each varied parameter:
       R_theta(config) = (mean_+ - mean_-) / (2 * mean_fid)
   with errors propagated from the realization scatter.
3. Same for the 4PCF on the equilateral-tetrahedra diagonal.

Usage: python analyze_responses.py [--results /home/csabiu/data/quijote/results]
"""
import argparse
import glob
import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

VARIATIONS = [('Om', r'$\Omega_{\rm m}$'), ('s8', r'$\sigma_8$'), ('h', r'$h$')]


def load_stack(results, cosmo, suffix):
    # column layout: 3pcf = 6 edges | NNN(6) RRR(7) zeta;  4pcf = 12 edges |
    # NNNN(12) RRRR(13) zeta ...  (nvfortran wraps long headers over multiple
    # lines, so skip all non-numeric leading rows)
    ncol_edges, c_n, c_r = (6, 6, 7) if suffix == '3pcf' else (12, 12, 13)
    files = sorted(glob.glob(os.path.join(results, f'{cosmo}_r*.{suffix}')))
    stack, edges = [], None
    for fn in files:
        n_skip = 0
        with open(fn) as f:
            for line in f:
                s = line.strip()
                if s and (s[0].isdigit() or s[0] in '+-'):
                    break
                n_skip += 1
        x = np.loadtxt(fn, skiprows=n_skip)
        if edges is None:
            edges = x[:, :ncol_edges]
        with np.errstate(divide='ignore', invalid='ignore'):
            stack.append(x[:, c_n] / x[:, c_r])
    return edges, np.array(stack)


def cuts_3pcf(edges):
    r1c = 0.5 * (edges[:, 0] + edges[:, 1])
    r2c = 0.5 * (edges[:, 2] + edges[:, 3])
    r3c = 0.5 * (edges[:, 4] + edges[:, 5])
    eq = (np.abs(r1c - r2c) < 1e-6) & (np.abs(r2c - r3c) < 1e-6)
    return r1c, r2c, r3c, eq


def slice_mask(edges, target1=50.0, target2=80.0):
    """Isoceles slice: configurations containing the bins nearest target1
    and target2; returns (r3 values, row indices) sorted by r3."""
    r1c = 0.5 * (edges[:, 0] + edges[:, 1])
    r2c = 0.5 * (edges[:, 2] + edges[:, 3])
    r3c = 0.5 * (edges[:, 4] + edges[:, 5])
    centers = np.unique(np.concatenate([r1c, r2c, r3c]))
    t1 = centers[np.argmin(np.abs(centers - target1))]
    t2 = centers[np.argmin(np.abs(centers - target2))]
    rows, rvals = [], []
    for i in range(edges.shape[0]):
        vals = [r1c[i], r2c[i], r3c[i]]
        used = [False] * 3
        ok = True
        for t in (t1, t2):
            for j in range(3):
                if not used[j] and abs(vals[j] - t) < 1e-6:
                    used[j] = True
                    break
            else:
                ok = False
                break
        if ok:
            rows.append(i)
            rvals.append(vals[used.index(False)])
    o = np.argsort(rvals)
    return np.array(rvals)[o], np.array(rows)[o], (t1, t2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--results', default='/home/csabiu/data/quijote/results')
    ap.add_argument('--outdir', default=None)
    args = ap.parse_args()
    outdir = args.outdir or args.results

    # ---------------- 3PCF ----------------
    edges, fid = load_stack(args.results, 'fiducial', '3pcf')
    nfid = fid.shape[0]
    fid_mean = np.nanmean(fid, axis=0)
    fid_err  = np.nanstd(fid, axis=0, ddof=1) / np.sqrt(nfid)   # error on mean
    fid_sig  = np.nanstd(fid, axis=0, ddof=1)                   # per-realization
    print(f'fiducial 3PCF: {nfid} realizations')

    r1c, r2c, r3c, eq = cuts_3pcf(edges)
    r_eq = r1c[eq]
    o = np.argsort(r_eq)

    # Figure 1: fiducial 3PCF — equilateral + isoceles slice
    r_sl, idx_sl, (t1, t2) = slice_mask(edges)
    fig, axes2 = plt.subplots(2, 1, figsize=(4.6, 7.0))
    ax = axes2[0]
    ax.errorbar(r_eq[o], (r_eq**2 * fid_mean[eq])[o],
                yerr=(r_eq**2 * fid_sig[eq])[o],
                fmt='o-', capsize=3, label=f'fiducial mean ({nfid} realizations)')
    ax.axhline(0, color='k', lw=0.5, alpha=0.5)
    ax.set_xlabel(r'$r$ [$h^{-1}$Mpc]')
    ax.set_ylabel(r'$r^2\,\zeta(r,r,r)$')
    ax.set_title('equilateral')
    ax.legend()
    ax = axes2[1]
    ax.errorbar(r_sl, fid_mean[idx_sl], yerr=fid_sig[idx_sl],
                fmt='s-', color='C3', capsize=3)
    ax.axvline(105, color='gray', ls=':', alpha=0.8, label='BAO ~105 Mpc/h')
    ax.axhline(0, color='k', lw=0.5, alpha=0.5)
    ax.set_xlabel(r'$r_3$ [$h^{-1}$Mpc]')
    ax.set_ylabel(rf'$\zeta({t1:.0f}, {t2:.0f}, r_3)$')
    ax.set_title(f'isoceles slice ({t1:.0f}, {t2:.0f}, $r_3$)')
    ax.legend()
    plt.tight_layout()
    fig.savefig(os.path.join(outdir, 'quijote_fiducial_3pcf.png'), dpi=130)
    fig.savefig(os.path.join(outdir, 'quijote_fiducial_3pcf.pdf'))
    print('wrote quijote_fiducial_3pcf.png')

    # Figure 2: responses on equilateral (top row) and isoceles slice (bottom)
    fig, axes = plt.subplots(2, len(VARIATIONS), figsize=(13, 7.6),
                             sharey='row')
    for col, (par, label) in enumerate(VARIATIONS):
        _, p = load_stack(args.results, f'{par}_p', '3pcf')
        _, m = load_stack(args.results, f'{par}_m', '3pcf')
        if p.size == 0 or m.size == 0:
            axes[0, col].set_title(f'{label} (missing)')
            continue
        dmean = np.nanmean(p, axis=0) - np.nanmean(m, axis=0)
        derr = np.sqrt(np.nanvar(p, axis=0, ddof=1)/p.shape[0] +
                       np.nanvar(m, axis=0, ddof=1)/m.shape[0])
        with np.errstate(divide='ignore', invalid='ignore'):
            resp = dmean / fid_sig       # response per unit measurement scatter
            resp_err = np.abs(derr / fid_sig)
        ax = axes[0, col]
        ax.errorbar(r_eq[o], resp[eq][o], yerr=resp_err[eq][o],
                    fmt='o-', capsize=3)
        ax.axhline(0, color='k', lw=0.5, alpha=0.5)
        ax.set_title(label)
        ax = axes[1, col]
        ax.errorbar(r_sl, resp[idx_sl], yerr=resp_err[idx_sl],
                    fmt='s-', color='C3', capsize=3)
        ax.axhline(0, color='k', lw=0.5, alpha=0.5)
        ax.set_xlabel(r'$r$ or $r_3$ [$h^{-1}$Mpc]')
    axes[0, 0].set_ylabel(
        r'equilateral: $(\zeta_{+} - \zeta_{-})\,/\,\sigma_{\rm fid}$')
    axes[1, 0].set_ylabel(
        rf'slice $({t1:.0f}, {t2:.0f}, r_3)$: '
        r'$(\zeta_{+} - \zeta_{-})\,/\,\sigma_{\rm fid}$')
    fig.suptitle('3PCF response to parameter variations')
    plt.tight_layout()
    fig.savefig(os.path.join(outdir, 'quijote_3pcf_responses.png'), dpi=130)
    fig.savefig(os.path.join(outdir, 'quijote_3pcf_responses.pdf'))
    print('wrote quijote_3pcf_responses.png')

    # ---------------- 4PCF ----------------
    try:
        e4, fid4 = load_stack(args.results, 'fiducial', '4pcf')
    except Exception:
        print('no 4PCF outputs yet'); return
    if fid4.size == 0:
        print('no 4PCF outputs yet'); return
    n4 = fid4.shape[0]
    f4_mean = np.nanmean(fid4, axis=0)
    f4_sig  = np.nanstd(fid4, axis=0, ddof=1)

    # equilateral tetrahedra: all six edge bins identical
    c4 = 0.5 * (e4[:, 0::2] + e4[:, 1::2])      # six edge-bin centres
    eq4 = np.all(np.abs(c4 - c4[:, :1]) < 1e-6, axis=1)
    r4 = c4[eq4, 0]
    o4 = np.argsort(r4)

    fig, axes = plt.subplots(1, len(VARIATIONS), figsize=(13, 4.4), sharey=True)
    for ax, (par, label) in zip(axes, VARIATIONS):
        _, p = load_stack(args.results, f'{par}_p', '4pcf')
        _, m = load_stack(args.results, f'{par}_m', '4pcf')
        if p.size == 0 or m.size == 0:
            ax.set_title(f'{label} (missing)')
            continue
        dmean = np.nanmean(p, axis=0) - np.nanmean(m, axis=0)
        derr = np.sqrt(np.nanvar(p, axis=0, ddof=1)/p.shape[0] +
                       np.nanvar(m, axis=0, ddof=1)/m.shape[0])
        with np.errstate(divide='ignore', invalid='ignore'):
            resp = dmean / f4_sig
            resp_err = np.abs(derr / f4_sig)
        ax.errorbar(r4[o4], resp[eq4][o4], yerr=resp_err[eq4][o4],
                    fmt='s-', capsize=3, color='C3')
        ax.axhline(0, color='k', lw=0.5, alpha=0.5)
        ax.set_xlabel(r'$r$ [$h^{-1}$Mpc]')
        ax.set_title(label)
    axes[0].set_ylabel(r'$(\zeta^{(4)}_{+} - \zeta^{(4)}_{-})\,/\,\sigma^{(4)}_{\rm fid}$')
    fig.suptitle(f'4PCF response to parameter variations '
                 f'(equilateral tetrahedra, {n4} fiducial realizations)')
    plt.tight_layout()
    fig.savefig(os.path.join(outdir, 'quijote_4pcf_responses.png'), dpi=130)
    fig.savefig(os.path.join(outdir, 'quijote_4pcf_responses.pdf'))
    print('wrote quijote_4pcf_responses.png')


if __name__ == '__main__':
    main()
