"""Figure 8: DESI DR1 LRG 4PCF vs the EZmock ensemble, all configurations,
NGC+SGC and both redshift bins combined at the count level.

Top:    full 4PCF zeta = N4/R4 (connected + disconnected), symlog axis.
Middle: connected 4PCF zeta_conn = zeta - zeta_disc, symlog axis.
Bottom: the connected residual zeta_conn(DESI) - mean(zeta_conn,mock) with
        the +/-1 sigma mock band.

Reads ../data/ezmock_4pcf_allcomb.csv; writes
../figures/desi_vs_ezmock_4pcf.{pdf,png}.  Style matches the 3PCF figure.
"""
import csv

import numpy as np
from _style import datafile, save, plt

rows = list(csv.DictReader(open(datafile('ezmock_4pcf_allcomb.csv'))))
idx = np.array([int(r['config']) for r in rows])
good = np.array([int(r['valid']) for r in rows]).astype(bool)


def col(name):
    a = np.array([float(r[name]) if r[name] != 'nan' else np.nan for r in rows])
    a[~good] = np.nan          # invalid configs -> gaps in every panel
    return a


df, fm, fs = col('desi_full'), col('full_mean'), col('full_sig')
dc, cm, cs = col('desi_conn'), col('conn_mean'), col('conn_sig')
nconf = len(idx)

fig = plt.figure(figsize=(7.0, 5.4))
gs = fig.add_gridspec(3, 1, height_ratios=[3, 3, 2], hspace=0.13,
                      left=0.115, right=0.975, top=0.975, bottom=0.085)
ax0 = fig.add_subplot(gs[0])
ax1 = fig.add_subplot(gs[1], sharex=ax0)
ax2 = fig.add_subplot(gs[2], sharex=ax0)


def band_panel(ax, desi, mean, sig, ylabel, label, linthresh):
    ax.fill_between(idx, mean - sig, mean + sig, color='C0', alpha=0.3, lw=0,
                    step='mid', label='EZmock mean $\\pm 1\\sigma$ (25 mocks)')
    ax.plot(idx, mean, '-', color='C0', lw=0.7, drawstyle='steps-mid')
    ax.plot(idx, desi, 'k.', ms=3, label='DESI DR1 LRG')
    ax.axhline(0, color='gray', lw=0.5, alpha=0.6)
    ax.set_yscale('symlog', linthresh=linthresh)
    ax.set_ylabel(ylabel)
    ax.text(0.985, 0.9, label, transform=ax.transAxes, ha='right', va='top',
            fontsize=9)
    ax.set_xlim(0.5, nconf + 0.5)


band_panel(ax0, df, fm, fs, r'$\zeta^{(4)}$', 'full', 1e-3)
ax0.tick_params(labelbottom=False)
ax0.legend(loc='lower right', frameon=False)

band_panel(ax1, dc, cm, cs, r'$\zeta^{(4)}_{\rm conn}$', 'connected', 3e-4)
ax1.tick_params(labelbottom=False)

# residual of the connected 4PCF
diff = dc - cm
ax2.fill_between(idx, -cs, cs, color='gray', alpha=0.25, lw=0, step='mid',
                 label=r'$\pm 1\sigma$ (mocks)')
ax2.plot(idx, diff, 'k.', ms=3)
ax2.axhline(0, color='gray', lw=0.5, alpha=0.6)
ax2.set_ylabel(r'$\zeta^{(4)}_{\rm conn}$ residual')
ax2.set_xlabel('configuration index')
m = good & np.isfinite(diff) & (cs > 0)
ax2.set_ylim(-3 * np.nanmedian(cs[m]), 3 * np.nanmedian(cs[m]))
ax2.legend(loc='upper right', frameon=False)

chi2 = float(np.sum((diff[m] / cs[m])**2))
print(f'  connected chi2/dof = {chi2:.1f}/{int(m.sum())} = {chi2 / m.sum():.2f}')
save(fig, 'desi_vs_ezmock_4pcf')
