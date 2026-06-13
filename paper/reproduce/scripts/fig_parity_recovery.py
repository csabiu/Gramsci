"""Figure 6: parity-estimator signal-injection recovery.
Reads ../data/parity_recovery.json; writes ../figures/parity_recovery.{pdf,png}.
"""
import json

import numpy as np
from _style import datafile, save, plt

rows = json.load(open(datafile('parity_recovery.json')))
inj = np.array([r['injected'] for r in rows])
mea = np.array([r['measured'] for r in rows])
sca = np.array([r['scatter'] for r in rows])

fig, ax = plt.subplots(figsize=(3.6, 3.3))
ax.plot([0, 1], [0, 1], '-', color='gray', lw=0.8, label='expected (1:1)')
ax.errorbar(inj, mea, yerr=sca, fmt='o', color='C3', ms=5, capsize=3,
            label='measured (500 tetrahedra, 3 seeds)')
ax.set_xlabel(r'injected asymmetry  $(N_L - N_R)/N_{\rm struct}$')
ax.set_ylabel(r'recovered  $\mathrm{NNNN}^{-} N_{\rm gal}^{4}/N_{\rm struct}$')
ax.legend(loc='upper left', frameon=False)
plt.tight_layout()
save(fig, 'parity_recovery')
print(f'  max |measured - injected| = {np.max(np.abs(mea - inj)):.4f}')
