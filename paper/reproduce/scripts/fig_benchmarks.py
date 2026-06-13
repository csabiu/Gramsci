"""Figures 2-5: three-generation scaling, thread scaling, density, out-of-core.
Reads ../data/benchmarks{,2}.json; writes ../figures/bench_*.{pdf,png}.
"""
import json

import numpy as np
from _style import datafile, save, plt

B1 = json.load(open(datafile('benchmarks.json')))
B2 = json.load(open(datafile('benchmarks2.json')))
NTH = B1['n_threads']


def _arr(sweep, side, key='query_s'):
    x, y = [], []
    for row in sweep:
        if row.get(side) and row[side].get(key) is not None:
            x.append(row['rmax']); y.append(row[side][key])
    return np.array(x), np.array(y)


# --- Fig 2: three generations vs rmax ---------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.1))
s3, s4 = B1['sweeps']['3pcf_rmax'], B1['sweeps']['4pcf_rmax']
sb = B2['sweeps']['bsearch_rmax']
for ax, name, sweep, bkey, slope, rg0, amp in [
        (axes[0], '3PCF', s3, 'cpu_bsearch_3pcf', 6, 60, 1.6),
        (axes[1], '4PCF', s4, 'cpu_bsearch_4pcf', 9, 45, 12.0)]:
    xb = np.array([r['rmax'] for r in sb if r[bkey] is not None])
    yb = np.array([r[bkey] for r in sb if r[bkey] is not None])
    xc, yc = _arr(sweep, 'cpu'); xg, yg = _arr(sweep, 'gpu')
    ax.plot(xb, yb, '^:', color='gray', mfc='white', ms=5, lw=1,
            label=f'v1 binary search, CPU$\\times${NTH}')
    ax.plot(xc, yc, 'o--', color='C0', mfc='white', ms=5, lw=1,
            label=f'v2 merge-walk, CPU$\\times${NTH}')
    ax.plot(xg, yg, 'o-', color='C3', ms=5, lw=1.3, label='GPU (this work)')
    r = np.array([rg0, max(xb.max(), xg.max())])
    ax.plot(r, amp * (r / r[-1])**slope, color='gray', lw=0.8, alpha=0.6)
    ax.text(r[0]*1.02, amp*(r[0]/r[-1])**slope*0.75,
            f'$\\propto r_{{\\rm max}}^{{{slope}}}$', fontsize=8, color='gray')
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel(r'$r_{\rm max}$  [$h^{-1}$Mpc]'); ax.set_title(name, fontsize=10)
    ticks = sorted(set(list(xb) + list(xg)))
    ax.set_xticks(ticks); ax.set_xticklabels([str(int(t)) for t in ticks], fontsize=7)
    ax.minorticks_off()
axes[0].set_ylabel('query time  [s]')
axes[0].legend(loc='upper left', frameon=False)
plt.tight_layout(); save(fig, 'bench_scaling')

# --- Fig 3: thread scaling --------------------------------------------------
fig, ax = plt.subplots(figsize=(3.6, 3.1))
th = B2['sweeps']['threads']
t = np.array([r['threads'] for r in th['cpu']])
y3 = np.array([r['cpu_3pcf'] for r in th['cpu']])
y4 = np.array([r['cpu_4pcf'] for r in th['cpu']])
ax.plot(t, y3, 'o-', color='C0', ms=5, lw=1.2, label=r'3PCF, $r_{\rm max}{=}70$')
ax.plot(t, y4, 's-', color='C3', ms=5, lw=1.2, label=r'4PCF, $r_{\rm max}{=}40$')
ax.plot(t, y3[0]/t, ':', color='gray', lw=0.9)
ax.plot(t, y4[0]/t, ':', color='gray', lw=0.9)
ax.text(t[-1]*0.8, y3[0]/t[-1]*0.52, 'ideal', fontsize=8, color='gray', ha='right')
ax.axhline(th['gpu_3pcf_rmax70'], color='C0', ls='--', lw=1)
ax.axhline(th['gpu_4pcf_rmax40'], color='C3', ls='--', lw=1)
ax.text(1.05, th['gpu_3pcf_rmax70']*1.12, 'GPU', fontsize=8, color='C0')
ax.text(1.05, th['gpu_4pcf_rmax40']*1.12, 'GPU', fontsize=8, color='C3')
ax.set_xscale('log', base=2); ax.set_yscale('log')
ax.set_xticks(t); ax.set_xticklabels([str(int(x)) for x in t])
ax.set_xlabel('OpenMP threads'); ax.set_ylabel('query time  [s]')
ax.legend(loc='upper right', frameon=False)
plt.tight_layout(); save(fig, 'bench_threads')

# --- Fig 4: density ---------------------------------------------------------
fig, ax = plt.subplots(figsize=(3.6, 3.1))
d = B2['sweeps']['density']
f = np.array([r['frac'] for r in d])
for key, color, marker, lab in [
        ('cpu_3pcf', 'C0', 'o', f'3PCF, CPU$\\times${NTH}'),
        ('gpu_3pcf', 'C0', 'o', '3PCF, GPU'),
        ('cpu_4pcf', 'C3', 's', f'4PCF, CPU$\\times${NTH}'),
        ('gpu_4pcf', 'C3', 's', '4PCF, GPU')]:
    y = np.array([r[key] for r in d])
    if key.startswith('cpu'):
        ax.plot(f, y, marker+'--', color=color, mfc='white', ms=5, lw=1, label=lab)
    else:
        ax.plot(f, y, marker+'-', color=color, ms=5, lw=1.3, label=lab)
g = np.array([0.1, 1.0])
ax.plot(g, 9.5*g**3, color='gray', lw=0.9, alpha=0.8)
ax.text(g[0], 9.5*g[0]**3*1.6, r'$\propto f^{3}$', fontsize=9, color='gray')
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xticks(f); ax.set_xticklabels(['1/8', '1/4', '1/2', '1']); ax.minorticks_off()
ax.set_xlabel('catalogue sampling fraction $f$'); ax.set_ylabel('query time  [s]')
ax.legend(loc='upper left', frameon=False)
plt.tight_layout(); save(fig, 'bench_density')

# --- Fig 5: out-of-core overhead --------------------------------------------
fig, ax = plt.subplots(figsize=(3.6, 2.7))
ch = B1['sweeps']['3pcf_chunks']
nw = np.array([c['nwin'] for c in ch])
tt = np.array([c['res']['query_s'] for c in ch])
edges = ch[0]['res']['edges']
ax.plot(nw, tt/tt[0], 'o-', color='C0', ms=6, lw=1.3,
        label=f'3PCF, $r_{{\\rm max}}{{=}}80$ ({edges/1e6:.0f}M edges)')
ax.axhline(1.0, color='gray', lw=0.8, ls=':')
ax.set_xlabel('number of edge windows $W$')
ax.set_ylabel('time relative to all-resident')
ax.set_xticks(nw); ax.set_ylim(0.9, 1.3)
ax.legend(loc='upper left', frameon=False)
plt.tight_layout(); save(fig, 'bench_chunks')
