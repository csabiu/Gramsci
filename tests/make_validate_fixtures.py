"""Generate the uniform test catalogues that src_gpu/validate.sh expects.

    python tests/make_validate_fixtures.py [n_gal] [n_ran] [box]

Writes tests/test.gal and tests/test.ran (deterministic, seed 5).
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
n_gal = int(sys.argv[1]) if len(sys.argv) > 1 else 40_000
n_ran = int(sys.argv[2]) if len(sys.argv) > 2 else 80_000
box = float(sys.argv[3]) if len(sys.argv) > 3 else 500.0

rng = np.random.default_rng(5)
for name, n in (('test.gal', n_gal), ('test.ran', n_ran)):
    path = os.path.join(HERE, name)
    np.savetxt(path, np.column_stack([rng.random((n, 3)) * box, np.ones(n)]),
               fmt='%.6e')
    print(f'wrote {path}  ({n} points, box {box:g})')
