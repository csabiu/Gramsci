"""Shared plot style and path helpers for the GRAMSCI reproducibility figures."""
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
FIGS = os.path.join(HERE, '..', 'figures')

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


def datafile(name):
    return os.path.join(DATA, name)


def save(fig, name):
    os.makedirs(FIGS, exist_ok=True)
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(FIGS, f'{name}.{ext}'))
    print(f'  wrote figures/{name}.pdf')
