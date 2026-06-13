"""Figure 7: DESI DR1 LRG 3PCF vs the EZmock ensemble, isoceles slice.
Reads ../data/ezmock_3pcf_slice.csv; writes ../figures/desi_vs_ezmock_3pcf.{pdf,png}.
"""
import csv

import numpy as np
from _style import datafile, save, plt

rows = list(csv.DictReader(open(datafile('ezmock_3pcf_slice.csv'))))
zbins = {1: [], 2: []}
for r in rows:
    zbins[int(r['zbin'])].append(
        (float(r['r3']), float(r['desi']), float(r['mock_mean']), float(r['mock_sigma'])))

fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.1))
titles = {1: r'$0.4 < z \leq 0.8$', 2: r'$0.8 < z \leq 1.1$'}
for ax, zbin in zip(axes, (1, 2)):
    a = np.array(zbins[zbin])
    r3, desi, mean, sig = a[:, 0], a[:, 1], a[:, 2], a[:, 3]
    ax.fill_between(r3, mean - sig, mean + sig, color='C0', alpha=0.25, lw=0,
                    label='EZmock mean $\\pm 1\\sigma$ (24 mocks)')
    ax.plot(r3, mean, '-', color='C0', lw=1.2)
    ax.errorbar(r3, desi, yerr=sig, fmt='o', color='k', ms=4, capsize=2.5, lw=1,
                label='DESI DR1 LRG')
    ax.axhline(0, color='gray', lw=0.5, alpha=0.6)
    ax.set_xlabel(r'$r_3$  [$h^{-1}$Mpc]'); ax.set_title(titles[zbin], fontsize=10)
    m = sig > 0
    chi2 = float(np.sum(((desi - mean)[m] / sig[m])**2))
    print(f'  zbin{zbin}: chi2/dof = {chi2:.1f}/{int(m.sum())}')
axes[0].set_ylabel(r'$\zeta(55, 75, r_3)$')
axes[0].legend(loc='upper left', frameon=False)
plt.tight_layout()
save(fig, 'desi_vs_ezmock_3pcf')
