#!/usr/bin/env python3
"""Generate synthetic test catalogs for GRAMSCI physics regression tests.

Each generator places known structures (pairs, tetrahedra) at well-separated
grid positions in a cubic box, alongside a uniform random catalog.  The grid
spacing (50) exceeds 2*rmax (40), so structures cannot cross-contaminate.

GRAMSCI input format: 4 columns (x y z weight), whitespace-separated.
"""

import numpy as np

# ---------------------------------------------------------------------------
# Shared constants -- must match the test runner parameters
# ---------------------------------------------------------------------------
RMIN = 5.0
RMAX = 20.0
NBINS = 3
R0 = 12.5       # target separation, center of bin 2
BOX_SIZE = 500.0
N_RAND = 50000
SPACING = 50.0   # grid spacing between structure centers
SEED = 42


def _write_catalog(filename, positions, weights):
    """Write a GRAMSCI-format catalog file (x y z weight)."""
    with open(filename, 'w') as f:
        for i in range(len(positions)):
            f.write(f"  {positions[i, 0]:14.7e}   {positions[i, 1]:14.7e}   "
                    f"{positions[i, 2]:14.7e}   {weights[i]:14.7e}\n")


def _make_grid_centers(n_struct, box_size, spacing, rng):
    """Create grid-based placement centers with small jitter.

    Returns n_struct positions on a regular grid with +-5 unit jitter.
    The grid spacing ensures minimum separation > 2*rmax between structures.
    """
    n_per_side = int(box_size / spacing)
    grid_1d = np.arange(n_per_side) * spacing + spacing / 2
    gx, gy, gz = np.meshgrid(grid_1d, grid_1d, grid_1d)
    all_centers = np.column_stack([gx.ravel(), gy.ravel(), gz.ravel()])

    if n_struct > len(all_centers):
        raise ValueError(f"Need {n_struct} grid cells but only {len(all_centers)} available")

    idx = rng.choice(len(all_centers), size=n_struct, replace=False)
    centers = all_centers[idx]
    centers += rng.uniform(-5.0, 5.0, size=centers.shape)
    return centers


def _make_randoms(n_rand, box_size, rng):
    """Generate uniform random positions in [0, box_size]^3."""
    positions = rng.uniform(0, box_size, size=(n_rand, 3))
    weights = np.ones(n_rand)
    return positions, weights


def generate_pairs(n_struct, r0, box_size, n_rand, seed, gal_file, ran_file):
    """Generate catalog of isolated pairs with separation r0.

    Each pair is oriented along the x-axis: center +/- (r0/2, 0, 0).
    """
    rng = np.random.default_rng(seed)
    centers = _make_grid_centers(n_struct, box_size, SPACING, rng)

    gal_positions = np.empty((2 * n_struct, 3))
    for i, c in enumerate(centers):
        gal_positions[2 * i] = c + np.array([r0 / 2, 0, 0])
        gal_positions[2 * i + 1] = c - np.array([r0 / 2, 0, 0])

    gal_weights = np.ones(2 * n_struct)

    ran_positions, ran_weights = _make_randoms(n_rand, box_size, rng)

    _write_catalog(gal_file, gal_positions, gal_weights)
    _write_catalog(ran_file, ran_positions, ran_weights)


def generate_regular_tetra(n_struct, r0, box_size, n_rand, seed, gal_file, ran_file):
    """Generate catalog of regular tetrahedra with edge length r0.

    Vertices are s*(+/-1, +/-1, +/-1) with s = r0 / (2*sqrt(2)),
    giving all 6 edges equal to r0.
    """
    rng = np.random.default_rng(seed)
    s = r0 / (2 * np.sqrt(2))
    local_verts = s * np.array([
        [1, 1, 1],
        [1, -1, -1],
        [-1, 1, -1],
        [-1, -1, 1],
    ], dtype=float)

    centers = _make_grid_centers(n_struct, box_size, SPACING, rng)

    gal_positions = np.empty((4 * n_struct, 3))
    for i, c in enumerate(centers):
        for j in range(4):
            gal_positions[4 * i + j] = c + local_verts[j]

    gal_weights = np.ones(4 * n_struct)

    ran_positions, ran_weights = _make_randoms(n_rand, box_size, rng)

    _write_catalog(gal_file, gal_positions, gal_weights)
    _write_catalog(ran_file, ran_positions, ran_weights)


def generate_chiral_tetra(n_struct, box_size, n_rand, seed, gal_file, ran_file,
                          handedness='left'):
    """Generate catalog of chiral (scalene) tetrahedra.

    The tetrahedron is irregular with edge lengths:
        d12=7.5 (bin 1), d13=d14=d23=12.5 (bin 2), d24=d34=17.5 (bin 3)
    Sorted bin-tuple: (1, 2, 2, 2, 3, 3)

    Left-handed: signed volume > 0  (v4 has z > 0)
    Right-handed: signed volume < 0 (v4 has z < 0)
    Balanced: first half left, second half right.

    Parameters
    ----------
    handedness : str
        'left', 'right', or 'balanced'
    """
    rng = np.random.default_rng(seed)

    # Raw vertices (before centering)
    v1_raw = np.array([0.0, 0.0, 0.0])
    v2_raw = np.array([7.5, 0.0, 0.0])
    v3_raw = np.array([3.75, 11.924240, 0.0])
    v4_left_raw = np.array([-6.25, 2.227605, 10.593643])
    v4_right_raw = np.array([-6.25, 2.227605, -10.593643])

    # Center on center-of-mass
    com_left = (v1_raw + v2_raw + v3_raw + v4_left_raw) / 4
    left_verts = np.array([v1_raw - com_left, v2_raw - com_left,
                           v3_raw - com_left, v4_left_raw - com_left])

    com_right = (v1_raw + v2_raw + v3_raw + v4_right_raw) / 4
    right_verts = np.array([v1_raw - com_right, v2_raw - com_right,
                            v3_raw - com_right, v4_right_raw - com_right])

    centers = _make_grid_centers(n_struct, box_size, SPACING, rng)

    gal_positions = np.empty((4 * n_struct, 3))
    for i, c in enumerate(centers):
        if handedness == 'left':
            verts = left_verts
        elif handedness == 'right':
            verts = right_verts
        elif handedness == 'balanced':
            verts = left_verts if i < n_struct // 2 else right_verts
        elif isinstance(handedness, float):
            # fractional: first round(x*n) structures left, rest right
            verts = left_verts if i < round(handedness * n_struct) else right_verts
        else:
            raise ValueError(f"Unknown handedness: {handedness}")

        for j in range(4):
            gal_positions[4 * i + j] = c + verts[j]

    gal_weights = np.ones(4 * n_struct)

    ran_positions, ran_weights = _make_randoms(n_rand, box_size, rng)

    _write_catalog(gal_file, gal_positions, gal_weights)
    _write_catalog(ran_file, ran_positions, ran_weights)


def generate_uniform_box(n_points, box_size, seed, gal_file):
    """Generate a uniform (unclustered) catalog filling [0, box_size)^3.

    Used by the periodic-box analytic-randoms tests: for uniform data every
    correlation function should vanish and the data tuple counts should match
    the analytic RR/RRR/RRRR expectations.
    """
    rng = np.random.default_rng(seed)
    positions = rng.uniform(0, box_size, size=(n_points, 3))
    weights = np.ones(n_points)
    _write_catalog(gal_file, positions, weights)
