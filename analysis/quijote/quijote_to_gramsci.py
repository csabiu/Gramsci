"""Convert a Quijote FoF halo catalog to Gramsci input format.

Reads the Gadget FoF group_tab binary format directly (no Pylians
dependency), applies a fixed number-density cut (most massive halos first),
and writes <out>.gal plus a matching uniform random catalog <out>.ran.

Usage:
    python quijote_to_gramsci.py --fof-dir <dir-with-group_tab files> \
        --snapnum 2 --ndens 2e-4 --boxsize 1000.0 --out halos_r0 --seed 0

The FoF directory is e.g. .../Halos/fiducial/0/groups_002/ containing
group_tab_002.0 ... group_tab_002.7 file pieces.
"""
import argparse
import glob
import os
import struct

import numpy as np


def read_fof_groups(fof_dir, snapnum):
    """Read positions and masses from Quijote FoF group_tab file pieces.

    Format per piece (Gadget-2 FoF, no ids):
      int32 Ngroups, int32 TotNgroups, int32 Nids, int64 TotNids,
      int32 NTask, int32 Nprev? -- Quijote uses the standard P-Gadget3
      layout: Ngroups, TotNgroups, Nids, TotNids(int64), NFiles
      then arrays: GroupLen(int32), GroupOffset(int32),
      GroupMass(float32), GroupPos(float32 x 3), GroupVel(float32 x 3),
      GroupTLen(int32 x 6)?  -- we only need Len/Offset/Mass/Pos/Vel.
    """
    pieces = sorted(glob.glob(os.path.join(fof_dir, f'group_tab_{snapnum:03d}.*')),
                    key=lambda p: int(p.rsplit('.', 1)[1]))
    if not pieces:
        raise FileNotFoundError(f'no group_tab_{snapnum:03d}.* in {fof_dir}')

    pos_all, mass_all = [], []
    for p in pieces:
        with open(p, 'rb') as f:
            ngroups   = struct.unpack('<i', f.read(4))[0]
            totngroups = struct.unpack('<i', f.read(4))[0]
            nids      = struct.unpack('<i', f.read(4))[0]
            totnids   = struct.unpack('<q', f.read(8))[0]
            nfiles    = struct.unpack('<i', f.read(4))[0]
            if ngroups == 0:
                continue
            np.fromfile(f, dtype=np.int32, count=ngroups)        # GroupLen
            np.fromfile(f, dtype=np.int32, count=ngroups)        # GroupOffset
            mass = np.fromfile(f, dtype=np.float32, count=ngroups)
            pos = np.fromfile(f, dtype=np.float32, count=3 * ngroups)
            pos = pos.reshape(ngroups, 3)
            pos_all.append(pos)
            mass_all.append(mass)
    pos = np.concatenate(pos_all)
    mass = np.concatenate(mass_all)
    # Quijote FoF positions are in kpc/h -> convert to Mpc/h
    if pos.max() > 2000.0:
        pos = pos / 1e3
    return pos.astype(np.float64), mass.astype(np.float64)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--fof-dir', required=True)
    ap.add_argument('--snapnum', type=int, default=2)        # z=0.5
    ap.add_argument('--ndens', type=float, default=2e-4)     # (h/Mpc)^3
    ap.add_argument('--boxsize', type=float, default=1000.0) # Mpc/h
    ap.add_argument('--out', required=True)
    ap.add_argument('--seed', type=int, default=0)
    args = ap.parse_args()

    pos, mass = read_fof_groups(args.fof_dir, args.snapnum)

    ntarget = int(round(args.ndens * args.boxsize**3))
    if ntarget > len(pos):
        raise ValueError(f'density cut needs {ntarget} halos, catalog has {len(pos)}')
    order = np.argsort(mass)[::-1][:ntarget]   # most massive first
    pos = pos[order]

    w = np.ones(len(pos))
    np.savetxt(f'{args.out}.gal', np.column_stack([pos, w]), fmt='%16.8e')

    rng = np.random.default_rng(args.seed)
    rpos = rng.uniform(0.0, args.boxsize, size=(len(pos), 3))
    np.savetxt(f'{args.out}.ran', np.column_stack([rpos, w]), fmt='%16.8e')

    print(f'{args.out}: {len(pos)} halos (n={args.ndens:g}), '
          f'Mmin={mass[order].min():.3e}, box={args.boxsize} Mpc/h')


if __name__ == '__main__':
    main()
