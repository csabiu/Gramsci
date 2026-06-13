"""Figure 9: Quijote fiducial 3PCF (equilateral + isoceles cuts) with scatter.
Reads ../data/quijote_fiducial_3pcf.csv; writes ../figures/quijote_fiducial_3pcf.{pdf,png}.
"""
import csv

import numpy as np
from _style import datafile, save, plt

rows = list(csv.DictReader(open(datafile('quijote_fiducial_3pcf.csv'))))
eq = np.array([[float(r['r']), float(r['zeta_mean']), float(r['zeta_sigma'])]
               for r in rows if r['cut'] == 'equilateral'])
iso_rows = [r for r in rows if r['cut'].startswith('isoceles')]
iso = np.array([[float(r['r']), float(r['zeta_mean']), float(r['zeta_sigma'])]
                for r in iso_rows])
iso_label = iso_rows[0]['cut'].replace('isoceles_', '').replace('_', ', ')

fig, axes = plt.subplots(2, 1, figsize=(5.0, 6.2))

r, m, s = eq[:, 0], eq[:, 1], eq[:, 2]
axes[0].errorbar(r, r**2 * m, yerr=r**2 * s, fmt='o-', ms=4, capsize=3,
                 label='fiducial mean (100 realizations)')
axes[0].axhline(0, color='k', lw=0.5, alpha=0.5)
axes[0].set_xlabel(r'$r$ [$h^{-1}$Mpc]')
axes[0].set_ylabel(r'$r^2\,\zeta(r,r,r)$')
axes[0].set_title('equilateral')
axes[0].legend(frameon=False)

r, m, s = iso[:, 0], iso[:, 1], iso[:, 2]
axes[1].errorbar(r, m, yerr=s, fmt='s-', color='C3', ms=4, capsize=3)
axes[1].axvline(105, color='gray', ls=':', alpha=0.7, label='BAO ~105 Mpc/h')
axes[1].axhline(0, color='k', lw=0.5, alpha=0.5)
axes[1].set_xlabel(r'$r_3$ [$h^{-1}$Mpc]')
axes[1].set_ylabel(r'$\zeta(%s, r_3)$' % iso_label)
axes[1].set_title(r'isoceles slice (%s, $r_3$)' % iso_label)
axes[1].legend(frameon=False)

plt.tight_layout()
save(fig, 'quijote_fiducial_3pcf')
