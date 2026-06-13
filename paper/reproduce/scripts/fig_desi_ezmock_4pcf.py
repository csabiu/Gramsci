"""Figure 8: DESI DR1 LRG 4PCF vs the EZmock ensemble, all configurations.
Reads ../data/ezmock_4pcf_band.csv; writes ../figures/desi_vs_ezmock_4pcf.{pdf,png}.
"""
import csv

import numpy as np
from _style import datafile, save, plt

rows = list(csv.DictReader(open(datafile('ezmock_4pcf_band.csv'))))
zbins = {1: [], 2: []}
for r in rows:
    mean = float(r['mock_mean']) if r['mock_mean'] != 'nan' else np.nan
    zbins[int(r['zbin'])].append(
        (int(r['config']), float(r['desi']), mean,
         float(r['mock_sigma']), int(r['valid'])))

fig, axes = plt.subplots(2, 1, figsize=(7.0, 4.6), sharex=True)
titles = {1: r'$0.4 < z \leq 0.8$', 2: r'$0.8 < z \leq 1.1$'}
for ax, zbin in zip(axes, (1, 2)):
    a = np.array(zbins[zbin], dtype=float)
    idx, desi, mean, sig, good = a[:, 0], a[:, 1], a[:, 2], a[:, 3], a[:, 4].astype(bool)
    ax.fill_between(idx, mean - sig, mean + sig, color='C0', alpha=0.3, lw=0,
                    step='mid', label='EZmock mean $\\pm 1\\sigma$ (25 mocks)')
    ax.plot(idx, mean, '-', color='C0', lw=0.7, drawstyle='steps-mid')
    ax.plot(idx, desi, 'k.', ms=2.5, label='DESI DR1 LRG')
    ax.axhline(0, color='gray', lw=0.5, alpha=0.6)
    ax.set_ylabel(r'$\zeta^{(4)}$')
    ax.text(0.985, 0.92, titles[zbin], transform=ax.transAxes, ha='right',
            va='top', fontsize=9)
    m = good & np.isfinite(desi) & (np.abs(desi) < 1.0) & (sig > 0)
    chi2 = float(np.nansum(((desi - mean)[m] / sig[m])**2))
    print(f'  zbin{zbin}: chi2/dof = {chi2:.1f}/{int(m.sum())} = {chi2/m.sum():.2f}')
    lo = np.nanpercentile(mean - 3*sig, 2); hi = np.nanpercentile(mean + 3*sig, 98)
    ax.set_ylim(lo, hi)
axes[0].legend(loc='lower left', frameon=False, ncol=2)
axes[1].set_xlabel('configuration index')
plt.tight_layout()
save(fig, 'desi_vs_ezmock_4pcf')
