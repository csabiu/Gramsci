"""Publication figures from benchmarks.json (AAS single-column style).

Outputs (PDF for the paper, PNG for quick viewing):
  bench_scaling.pdf   query time vs rmax, 3PCF/4PCF x CPU/GPU
  bench_chunks.pdf    out-of-core overhead vs number of edge windows
"""
import json
import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams.update({
    'font.family': 'serif',
    'mathtext.fontset': 'stix',
    'font.size': 9,
    'axes.labelsize': 10,
    'legend.fontsize': 8,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.top': True,
    'ytick.right': True,
    'figure.dpi': 150,
})

HERE = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(HERE, 'benchmarks.json')) as f:
    B = json.load(f)
NTH = B['n_threads']


def arr(sweep, side, key='query_s'):
    x, y = [], []
    for row in sweep:
        if row.get(side) and row[side].get(key) is not None:
            x.append(row['rmax'])
            y.append(row[side][key])
    return np.array(x), np.array(y)


# ----------------------------------------------------------------------
# Figure 1: scaling with rmax
# ----------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(3.6, 3.2))

s3 = B['sweeps']['3pcf_rmax']
s4 = B['sweeps']['4pcf_rmax']

for sweep, color, label in [(s3, 'C0', '3PCF'), (s4, 'C3', '4PCF')]:
    xc, yc = arr(sweep, 'cpu')
    xg, yg = arr(sweep, 'gpu')
    ax.plot(xc, yc, 'o--', color=color, mfc='white', ms=5, lw=1,
            label=f'{label}, CPU$\\times${NTH}')
    ax.plot(xg, yg, 'o-', color=color, ms=5, lw=1.3,
            label=f'{label}, GPU')

# power-law guides
r = np.array([55, 100])
ax.plot(r, 2.0 * (r / 80.0)**6, color='gray', lw=0.8, alpha=0.7)
ax.text(r[-1], 2.0 * (r[-1] / 80.0)**6, r' $\propto r_{\rm max}^{6}$',
        fontsize=8, color='gray', va='center')
r = np.array([45, 80])
ax.plot(r, 25.0 * (r / 65.0)**9, color='gray', lw=0.8, alpha=0.7)
ax.text(r[-1], 25.0 * (r[-1] / 65.0)**9, r' $\propto r_{\rm max}^{9}$',
        fontsize=8, color='gray', va='center')

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$r_{\rm max}$  [$h^{-1}$Mpc]')
ax.set_ylabel('query time  [s]')
ax.set_xticks([30, 40, 50, 60, 70, 80, 100])
ax.set_xticklabels(['30', '40', '50', '60', '70', '80', '100'])
ax.legend(loc='upper left', frameon=False)
plt.tight_layout()
fig.savefig(os.path.join(HERE, 'bench_scaling.pdf'))
fig.savefig(os.path.join(HERE, 'bench_scaling.png'))
print('wrote bench_scaling.[pdf,png]')

# speedups at largest common rmax
for sweep, name in [(s3, '3PCF'), (s4, '4PCF')]:
    xc, yc = arr(sweep, 'cpu')
    xg, yg = arr(sweep, 'gpu')
    common = set(xc) & set(xg)
    rm = max(common)
    su = yc[list(xc).index(rm)] / yg[list(xg).index(rm)]
    print(f'{name}: speedup at rmax={rm}: {su:.1f}x')

# ----------------------------------------------------------------------
# Figure 2: out-of-core overhead
# ----------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(3.6, 2.9))

ch = B['sweeps']['3pcf_chunks']
nw = np.array([c['nwin'] for c in ch])
t = np.array([c['res']['query_s'] for c in ch])
t0 = t[nw == 1][0]
edges = ch[0]['res']['edges']

ax.plot(nw, t / t0, 'o-', color='C0', ms=6, lw=1.3,
        label=f'3PCF, $r_{{\\rm max}}{{=}}80$ ({edges/1e6:.0f}M edges)')
ax.axhline(1.0, color='gray', lw=0.8, ls=':')

ax.set_xlabel('number of edge windows $W$')
ax.set_ylabel('query time relative to all-resident')
ax.set_xticks(nw)
ax.set_ylim(bottom=0.9)
ax.legend(loc='upper left', frameon=False)
plt.tight_layout()
fig.savefig(os.path.join(HERE, 'bench_chunks.pdf'))
fig.savefig(os.path.join(HERE, 'bench_chunks.png'))
print('wrote bench_chunks.[pdf,png]')
print('overheads:', {int(n): round(float(x), 3) for n, x in zip(nw, t / t0)})
