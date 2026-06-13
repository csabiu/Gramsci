"""Figure 8: DESI DR1 LRG 4PCF vs the EZmock ensemble, all configurations,
with a configuration-scale key panel.

Top two panels: zeta^(4) vs configuration index for the two redshift bins
(DESI points over the EZmock mean +/- 1 sigma band).
Bottom panel: a key showing, for each configuration, the bin-centre scale
of all six tetrahedron edges (r12...r34) as a colour --- read a vertical
stripe to recover the physical scales of any configuration.

Reads ../data/ezmock_4pcf_band.csv and ../data/fourpcf_config_scales.csv;
writes ../figures/desi_vs_ezmock_4pcf.{pdf,png}.
"""
import csv

import numpy as np
from matplotlib.colors import BoundaryNorm, ListedColormap
from _style import datafile, save, plt

# --- ensemble band ----------------------------------------------------------
rows = list(csv.DictReader(open(datafile('ezmock_4pcf_band.csv'))))
zbins = {1: [], 2: []}
for r in rows:
    mean = float(r['mock_mean']) if r['mock_mean'] != 'nan' else np.nan
    zbins[int(r['zbin'])].append(
        (int(r['config']), float(r['desi']), mean,
         float(r['mock_sigma']), int(r['valid'])))

# --- per-config edge scales -------------------------------------------------
scl = list(csv.DictReader(open(datafile('fourpcf_config_scales.csv'))))
edges = ['r12', 'r13', 'r14', 'r23', 'r24', 'r34']
scales = np.array([[float(s[e]) for e in edges] for s in scl])   # (Nconfig, 6)
uniq = np.unique(np.round(scales, 1))                            # 4 discrete bins
nconf = scales.shape[0]

# --- explicit gridspec layout (no tight_layout: it is incompatible with a
#     manually-placed colorbar and would leave a blank top margin) -----------
fig = plt.figure(figsize=(7.0, 5.0))
gs = fig.add_gridspec(3, 1, height_ratios=[3, 3, 1.8], hspace=0.13,
                      left=0.085, right=0.905, top=0.975, bottom=0.085)
ax0 = fig.add_subplot(gs[0])
ax1 = fig.add_subplot(gs[1], sharex=ax0)
axk = fig.add_subplot(gs[2], sharex=ax0)
titles = {1: r'$0.4 < z \leq 0.8$', 2: r'$0.8 < z \leq 1.1$'}

for ax, zbin in zip((ax0, ax1), (1, 2)):
    a = np.array(zbins[zbin], dtype=float)
    idx, desi, mean, sig, good = a[:, 0], a[:, 1], a[:, 2], a[:, 3], a[:, 4].astype(bool)
    ax.fill_between(idx, mean - sig, mean + sig, color='C0', alpha=0.3, lw=0,
                    step='mid', label='EZmock mean $\\pm 1\\sigma$ (25 mocks)')
    ax.plot(idx, mean, '-', color='C0', lw=0.7, drawstyle='steps-mid')
    ax.plot(idx, desi, 'k.', ms=2.5, label='DESI DR1 LRG')
    ax.axhline(0, color='gray', lw=0.5, alpha=0.6)
    ax.set_ylabel(r'$\zeta^{(4)}$')
    ax.text(0.985, 0.9, titles[zbin], transform=ax.transAxes, ha='right',
            va='top', fontsize=9)
    ax.tick_params(labelbottom=False)
    m = good & np.isfinite(desi) & (np.abs(desi) < 1.0) & (sig > 0)
    chi2 = float(np.nansum(((desi - mean)[m] / sig[m])**2))
    print(f'  zbin{zbin}: chi2/dof = {chi2:.1f}/{int(m.sum())} = {chi2/m.sum():.2f}')
    lo = np.nanpercentile(mean - 3 * sig, 2); hi = np.nanpercentile(mean + 3 * sig, 98)
    ax.set_ylim(lo, hi)
    ax.set_xlim(0.5, nconf + 0.5)

# --- key panel: heatmap of the six edge scales ------------------------------
cmap = ListedColormap(plt.cm.viridis(np.linspace(0.08, 0.92, len(uniq))))
bounds = np.concatenate([[uniq[0] - 1], (uniq[:-1] + uniq[1:]) / 2, [uniq[-1] + 1]])
norm = BoundaryNorm(bounds, cmap.N)
im = axk.imshow(scales.T, aspect='auto', cmap=cmap, norm=norm, origin='lower',
                extent=[0.5, nconf + 0.5, -0.5, 5.5], interpolation='nearest')
axk.set_yticks(range(6))
axk.set_yticklabels([r'$r_{12}$', r'$r_{13}$', r'$r_{14}$',
                     r'$r_{23}$', r'$r_{24}$', r'$r_{34}$'], fontsize=8)
axk.set_xlabel('configuration index')
axk.set_ylabel('edge', fontsize=9)

# colorbar placed exactly beside the heatmap (no width stolen from panels)
pos = axk.get_position()
cax = fig.add_axes([pos.x1 + 0.012, pos.y0, 0.013, pos.height])
cb = fig.colorbar(im, cax=cax, ticks=uniq)
cb.ax.set_yticklabels([f'{u:.0f}' for u in uniq], fontsize=7)
cb.set_label(r'$r$ [$h^{-1}$Mpc]', fontsize=8)

ax0.legend(loc='lower left', frameon=False, ncol=2)
save(fig, 'desi_vs_ezmock_4pcf')
