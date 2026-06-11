"""Plot the equilateral 3PCF from test_bao.3pcf to look for the BAO feature (~105 Mpc/h)."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

x = np.loadtxt('test_bao.3pcf', skiprows=1)
r1c = 0.5 * (x[:, 0] + x[:, 1])
r2c = 0.5 * (x[:, 2] + x[:, 3])
r3c = 0.5 * (x[:, 4] + x[:, 5])
zeta = x[:, 8]

# Equilateral configurations: r1 = r2 = r3
eq = (np.abs(r1c - r2c) < 1e-6) & (np.abs(r2c - r3c) < 1e-6)
r_eq = r1c[eq]
z_eq = zeta[eq]
order = np.argsort(r_eq)
r_eq, z_eq = r_eq[order], z_eq[order]

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

axes[0].plot(r_eq, z_eq, 'o-')
axes[0].axvline(105, color='gray', ls='--', alpha=0.6, label='BAO ~105 Mpc/h')
axes[0].set_xlabel('r [Mpc/h]')
axes[0].set_ylabel(r'$\zeta(r,r,r)$')
axes[0].set_title('Equilateral 3PCF')
axes[0].set_yscale('symlog', linthresh=1e-4)
axes[0].legend()

axes[1].plot(r_eq, r_eq**4 * z_eq, 'o-')
axes[1].axvline(105, color='gray', ls='--', alpha=0.6, label='BAO ~105 Mpc/h')
axes[1].set_xlabel('r [Mpc/h]')
axes[1].set_ylabel(r'$r^4\,\zeta(r,r,r)$')
axes[1].set_title('Equilateral 3PCF (BAO-weighted)')
axes[1].legend()

plt.tight_layout()
plt.savefig('bao_3pcf.png', dpi=130)
print('wrote bao_3pcf.png')
print()
print(f'{"r [Mpc/h]":>10s}  {"zeta(r,r,r)":>14s}  {"r^4 zeta":>12s}')
for r, z in zip(r_eq, z_eq):
    print(f'{r:10.2f}  {z:14.5e}  {r**4 * z:12.4e}')
