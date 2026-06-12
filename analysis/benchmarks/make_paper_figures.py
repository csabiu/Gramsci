"""Final paper figures from benchmarks.json + benchmarks2.json.

  bench_scaling.pdf   2 panels: 3PCF | 4PCF query time vs rmax,
                      three generations (v1 bsearch, v2 merge-walk, GPU)
  bench_threads.pdf   OpenMP strong scaling + GPU reference lines
  bench_density.pdf   query time vs catalog sampling fraction
  bench_chunks.pdf    out-of-core overhead vs window count
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
B1 = json.load(open(os.path.join(HERE, 'benchmarks.json')))
B2 = json.load(open(os.path.join(HERE, 'benchmarks2.json')))
NTH = B1['n_threads']


def save(fig, name):
    fig.savefig(os.path.join(HERE, name + '.pdf'))
    fig.savefig(os.path.join(HERE, name + '.png'))
    print(f'wrote {name}.[pdf,png]')


# ======================================================================
# Figure 1: three generations, query time vs rmax (two panels)
# ======================================================================
fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.1), sharex=False)

s3 = B1['sweeps']['3pcf_rmax']
s4 = B1['sweeps']['4pcf_rmax']
sb = B2['sweeps']['bsearch_rmax']

panels = [
    (axes[0], '3PCF', s3, 'cpu_bsearch_3pcf', 6, 60, 1.6),
    (axes[1], '4PCF', s4, 'cpu_bsearch_4pcf', 9, 45, 12.0),
]
for ax, name, sweep, bkey, slope, rg0, amp in panels:
    xb = np.array([r['rmax'] for r in sb if r[bkey] is not None])
    yb = np.array([r[bkey] for r in sb if r[bkey] is not None])
    xc = np.array([r['rmax'] for r in sweep if r.get('cpu')])
    yc = np.array([r['cpu']['query_s'] for r in sweep if r.get('cpu')])
    xg = np.array([r['rmax'] for r in sweep])
    yg = np.array([r['gpu']['query_s'] for r in sweep])

    ax.plot(xb, yb, '^:', color='gray', mfc='white', ms=5, lw=1,
            label=f'v1 binary search, CPU$\\times${NTH}')
    ax.plot(xc, yc, 'o--', color='C0', mfc='white', ms=5, lw=1,
            label=f'v2 merge-walk, CPU$\\times${NTH}')
    ax.plot(xg, yg, 'o-', color='C3', ms=5, lw=1.3,
            label='GPU (this work)')

    r = np.array([rg0, max(xb.max(), xg.max())])
    ax.plot(r, amp * (r / r[-1])**slope, color='gray', lw=0.8, alpha=0.6)
    ax.text(r[0] * 1.02, amp * (r[0] / r[-1])**slope * 0.75,
            f'$\\propto r_{{\\rm max}}^{{{slope}}}$',
            fontsize=8, color='gray')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$r_{\rm max}$  [$h^{-1}$Mpc]')
    ax.set_title(name, fontsize=10)
    ticks = sorted(set(list(xb) + list(xg)))
    ax.set_xticks(ticks)
    ax.set_xticklabels([str(int(t)) for t in ticks], fontsize=7)
    ax.minorticks_off()

axes[0].set_ylabel('query time  [s]')
axes[0].legend(loc='upper left', frameon=False)
plt.tight_layout()
save(fig, 'bench_scaling')

# lineage numbers
print('  3PCF rmax=80 lineage:', 29.995, '->', 9.448, '->', 3.733)
print('  4PCF rmax=70 lineage:', 1043.576, '->', 279.589, '->', 25.57,
      f'({1043.576/25.57:.0f}x v1->GPU)')

# ======================================================================
# Figure 2: OpenMP strong scaling + GPU reference
# ======================================================================
fig, ax = plt.subplots(figsize=(3.6, 3.1))

th = B2['sweeps']['threads']
t = np.array([r['threads'] for r in th['cpu']])
y3 = np.array([r['cpu_3pcf'] for r in th['cpu']])
y4 = np.array([r['cpu_4pcf'] for r in th['cpu']])

ax.plot(t, y3, 'o-', color='C0', ms=5, lw=1.2,
        label=r'3PCF, $r_{\rm max}{=}70$')
ax.plot(t, y4, 's-', color='C3', ms=5, lw=1.2,
        label=r'4PCF, $r_{\rm max}{=}40$')
ax.plot(t, y3[0] / t, ':', color='gray', lw=0.9)
ax.plot(t, y4[0] / t, ':', color='gray', lw=0.9)
ax.text(t[-1] * 0.8, y3[0] / t[-1] * 0.62, 'ideal', fontsize=8,
        color='gray', ha='right')

ax.axhline(th['gpu_3pcf_rmax70'], color='C0', ls='--', lw=1)
ax.axhline(th['gpu_4pcf_rmax40'], color='C3', ls='--', lw=1)
ax.text(1.05, th['gpu_3pcf_rmax70'] * 1.12, 'GPU', fontsize=8, color='C0')
ax.text(1.05, th['gpu_4pcf_rmax40'] * 0.72, 'GPU', fontsize=8, color='C3')

ax.set_xscale('log', base=2)
ax.set_yscale('log')
ax.set_xticks(t)
ax.set_xticklabels([str(int(x)) for x in t])
ax.set_xlabel('OpenMP threads')
ax.set_ylabel('query time  [s]')
ax.legend(loc='upper right', frameon=False)
plt.tight_layout()
save(fig, 'bench_threads')
print(f'  3PCF eff@{t[-1]}: {y3[0]/y3[-1]/t[-1]*100:.0f}%   '
      f'4PCF eff@{t[-1]}: {y4[0]/y4[-1]/t[-1]*100:.0f}%')

# ======================================================================
# Figure 3: density scaling
# ======================================================================
fig, ax = plt.subplots(figsize=(3.6, 3.1))

d = B2['sweeps']['density']
f = np.array([r['frac'] for r in d])
for key, color, marker, label in [
        ('cpu_3pcf', 'C0', 'o', f'3PCF, CPU$\\times${NTH}'),
        ('gpu_3pcf', 'C0', 'o', '3PCF, GPU'),
        ('cpu_4pcf', 'C3', 's', f'4PCF, CPU$\\times${NTH}'),
        ('gpu_4pcf', 'C3', 's', '4PCF, GPU')]:
    y = np.array([r[key] for r in d])
    if key.startswith('cpu'):
        ax.plot(f, y, marker + '--', color=color, mfc='white', ms=5, lw=1,
                label=label)
    else:
        ax.plot(f, y, marker + '-', color=color, ms=5, lw=1.3, label=label)

g = np.array([0.2, 1.0])
ax.plot(g, 9.5 * g**3, color='gray', lw=0.8, alpha=0.6)
ax.text(g[0], 9.5 * g[0]**3 * 1.4, r'$\propto f^{3}$', fontsize=8,
        color='gray')

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xticks(f)
ax.set_xticklabels(['1/8', '1/4', '1/2', '1'])
ax.minorticks_off()
ax.set_xlabel('catalogue sampling fraction $f$')
ax.set_ylabel('query time  [s]')
ax.legend(loc='upper left', frameon=False, ncol=1)
plt.tight_layout()
save(fig, 'bench_density')

# ======================================================================
# Figure 4: out-of-core overhead
# ======================================================================
fig, ax = plt.subplots(figsize=(3.6, 2.7))

ch = B1['sweeps']['3pcf_chunks']
nw = np.array([c['nwin'] for c in ch])
tt = np.array([c['res']['query_s'] for c in ch])
edges = ch[0]['res']['edges']

ax.plot(nw, tt / tt[0], 'o-', color='C0', ms=6, lw=1.3,
        label=f'3PCF, $r_{{\\rm max}}{{=}}80$ ({edges/1e6:.0f}M edges)')
ax.axhline(1.0, color='gray', lw=0.8, ls=':')

ax.set_xlabel('number of edge windows $W$')
ax.set_ylabel('time relative to all-resident')
ax.set_xticks(nw)
ax.set_ylim(0.9, 1.3)
ax.legend(loc='upper left', frameon=False)
plt.tight_layout()
save(fig, 'bench_chunks')
print('  chunk overheads:',
      {int(n): round(float(x), 3) for n, x in zip(nw, tt / tt[0])})
