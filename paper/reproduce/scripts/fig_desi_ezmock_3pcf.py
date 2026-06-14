"""Figure 7: DESI DR1 LRG 3PCF vs the EZmock ensemble, all configurations
(NGC+SGC and both redshift bins combined at the count level).

Top panel:    zeta vs configuration index, DESI points over the EZmock
              mean +/- 1 sigma band (symlog y-axis).
Middle panel: the residual zeta_DESI - mean(zeta_mock) with the +/-1 sigma
              mock band.
Bottom panel: a key to the configuration axis --- the three sorted
              triangle side-length bins (r1 <= r2 <= r3) for each
              configuration, so a vertical read-off gives the triangle
              shape.

Reads ../data/ezmock_3pcf_allconfig.csv and ../data/threepcf_config_scales.csv;
writes ../figures/desi_vs_ezmock_3pcf.{pdf,png}.
"""
import csv

import numpy as np
from _style import datafile, save, plt

rows = list(csv.DictReader(open(datafile('ezmock_3pcf_allconfig.csv'))))
idx = np.array([int(r['config']) for r in rows])
desi = np.array([float(r['desi']) for r in rows])
mean = np.array([float(r['mock_mean']) if r['mock_mean'] != 'nan' else np.nan for r in rows])
sig = np.array([float(r['mock_sigma']) for r in rows])
good = np.array([int(r['valid']) for r in rows]).astype(bool)

scl = list(csv.DictReader(open(datafile('threepcf_config_scales.csv'))))
r1 = np.array([float(s['r1']) for s in scl])
r2 = np.array([float(s['r2']) for s in scl])
r3 = np.array([float(s['r3']) for s in scl])
nconf = len(idx)

fig = plt.figure(figsize=(7.0, 5.2))
gs = fig.add_gridspec(3, 1, height_ratios=[3, 2, 1.8], hspace=0.13,
                      left=0.105, right=0.975, top=0.975, bottom=0.085)
ax0 = fig.add_subplot(gs[0])
ax1 = fig.add_subplot(gs[1], sharex=ax0)
axk = fig.add_subplot(gs[2], sharex=ax0)

# --- panel 1: 3PCF, DESI vs mock band ---------------------------------------
ax0.fill_between(idx, mean - sig, mean + sig, color='C0', alpha=0.3, lw=0,
                 step='mid', label='EZmock mean $\\pm 1\\sigma$ (24 mocks)')
ax0.plot(idx, mean, '-', color='C0', lw=0.7, drawstyle='steps-mid')
ax0.plot(idx, desi, 'k.', ms=3, label='DESI DR1 LRG')
ax0.axhline(0, color='gray', lw=0.5, alpha=0.6)
ax0.set_yscale('symlog', linthresh=1e-3)
ax0.set_ylabel(r'$\zeta$')
ax0.tick_params(labelbottom=False)
ax0.set_xlim(0.5, nconf + 0.5)
ax0.legend(loc='upper right', frameon=False)

# --- panel 2: residual ------------------------------------------------------
diff = desi - mean
ax1.fill_between(idx, -sig, sig, color='gray', alpha=0.25, lw=0, step='mid',
                 label=r'$\pm 1\sigma$ (mocks)')
ax1.plot(idx, diff, 'k.', ms=3)
ax1.axhline(0, color='gray', lw=0.5, alpha=0.6)
ax1.set_ylabel(r'$\zeta_{\rm DESI}-\bar\zeta_{\rm mock}$')
ax1.tick_params(labelbottom=False)
m = good & np.isfinite(diff) & (sig > 0)
ax1.set_ylim(-3 * np.nanmedian(sig[m]), 3 * np.nanmedian(sig[m]))
ax1.legend(loc='upper right', frameon=False)
chi2 = float(np.sum((diff[m] / sig[m])**2))
print(f'  combined z-range: chi2/dof = {chi2:.1f}/{int(m.sum())} = {chi2/m.sum():.2f}')

# --- panel 3: configuration key (3 sorted triangle sides) -------------------
for arr, lab, c in [(r3, r'$r_3$ (largest)', 'C3'),
                    (r2, r'$r_2$', 'C1'),
                    (r1, r'$r_1$ (smallest)', 'C2')]:
    axk.step(idx, arr, where='mid', color=c, lw=1.0, label=lab)
axk.set_ylabel(r'$r$ [$h^{-1}$Mpc]')
axk.set_xlabel('configuration index')
axk.set_ylim(48, 158)
axk.legend(loc='lower right', frameon=True, framealpha=0.9, ncol=3,
           handlelength=1.3, fontsize=7)

save(fig, 'desi_vs_ezmock_3pcf')
