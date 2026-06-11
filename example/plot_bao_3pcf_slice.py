"""Plot the 3PCF slice r1~50, r2~80, varying r3, from test_bao.3pcf.

The output file stores each configuration once with (r1,r2,r3) sorted, so we
match rows where the sorted triplet contains the two target bins and scan the
third side.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

x = np.loadtxt('test_bao.3pcf', skiprows=1)
r1c = 0.5 * (x[:, 0] + x[:, 1])
r2c = 0.5 * (x[:, 2] + x[:, 3])
r3c = 0.5 * (x[:, 4] + x[:, 5])
zeta = x[:, 8]

# Bin centers available (15 linear bins, 1..150)
centers = np.unique(np.concatenate([r1c, r2c, r3c]))
print('bin centers:', np.array2string(centers, precision=2))

# Bins closest to the requested sides
t1 = centers[np.argmin(np.abs(centers - 50.0))]   # ~45.70
t2 = centers[np.argmin(np.abs(centers - 80.0))]   # ~75.50 or 85.43
print(f'using r1 bin center = {t1:.2f}, r2 bin center = {t2:.2f}')

# Configurations are stored sorted; find rows whose triplet contains (t1, t2)
# and extract the remaining side.
r_scan, z_scan = [], []
for row in range(x.shape[0]):
    trip = [r1c[row], r2c[row], r3c[row]]
    used = [False, False, False]
    found1 = found2 = False
    for j in range(3):
        if not used[j] and abs(trip[j] - t1) < 1e-6:
            used[j] = True
            found1 = True
            break
    for j in range(3):
        if not used[j] and abs(trip[j] - t2) < 1e-6:
            used[j] = True
            found2 = True
            break
    if found1 and found2:
        j = used.index(False)
        r_scan.append(trip[j])
        z_scan.append(zeta[row])

r_scan = np.array(r_scan)
z_scan = np.array(z_scan)
order = np.argsort(r_scan)
r_scan, z_scan = r_scan[order], z_scan[order]

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(r_scan, z_scan, 'o-')
ax.axvline(105, color='gray', ls='--', alpha=0.6, label='BAO ~105 Mpc/h')
ax.axhline(0, color='k', lw=0.5, alpha=0.5)
ax.set_xlabel(r'$r_3$ [Mpc/h]')
ax.set_ylabel(r'$\zeta(r_1\!\sim\!%.0f,\ r_2\!\sim\!%.0f,\ r_3)$' % (t1, t2))
ax.set_title('3PCF slice: r1~%.0f, r2~%.0f Mpc/h, scanning r3' % (t1, t2))
ax.legend()
plt.tight_layout()
plt.savefig('bao_3pcf_slice.png', dpi=130)
print('wrote bao_3pcf_slice.png')
print()
print(f'{"r3 [Mpc/h]":>10s}  {"zeta":>14s}')
for r, z in zip(r_scan, z_scan):
    print(f'{r:10.2f}  {z:14.5e}')
