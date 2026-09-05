#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# validate_shuffle.sh — the 4PCF / 4PCF-parity outputs must not depend on the
# ORDER of the input rows.  Runs the GPU binary on the test catalogues and on
# a row-permuted copy and compares every output column.  A labelling-dependent
# parity sign (achiral bin configurations before the chiral_4pcf fix) shows up
# here as O(1) changes in zeta_odd; a correct estimator changes only at the
# floating-point summation-order level.
#
# Usage: cd tests && bash ../src_gpu/validate_shuffle.sh [rel_tol]
# ---------------------------------------------------------------------------
set -euo pipefail
BINDIR="$(dirname "$0")/../bin"; TESTDIR="$(dirname "$0")/../tests"
TOL="${1:-1e-3}"
GPU_BIN="${GPU_BIN:-$BINDIR/gramsci_gpu}"
GAL="$TESTDIR/test.gal"; RAN="$TESTDIR/test.ran"
if [[ ! -f "$GAL" ]]; then GAL="$(dirname "$0")/../example/test.gal"; RAN="$(dirname "$0")/../example/test.ran"; fi
[[ -x "$GPU_BIN" ]] || { echo "ERROR: $GPU_BIN not found"; exit 1; }
W=$(mktemp -d); trap 'rm -rf "$W"' EXIT
python3 - "$GAL" "$RAN" "$W" <<'PY'
import sys, numpy as np
gal, ran, w = sys.argv[1:4]; rng = np.random.default_rng(12345)
for src, dst in ((gal, f"{w}/shuf.gal"), (ran, f"{w}/shuf.ran")):
    hdr = [l for l in open(src) if l.startswith('#')]; rows = [l for l in open(src) if not l.startswith('#')]
    open(dst, 'w').writelines(hdr + [rows[i] for i in rng.permutation(len(rows))])
PY
PASS=0; FAIL=0
for flag in -4pcf -4pcfp "-4pcfp -exactparity"; do
  echo -n "shuffle invariance, $flag ... "
  "$GPU_BIN" -gal "$GAL" -ran "$RAN" -rmin 10 -rmax 30 -nbins 3 -njk 8 $flag -out "$W/a" >/dev/null 2>&1
  "$GPU_BIN" -gal "$W/shuf.gal" -ran "$W/shuf.ran" -rmin 10 -rmax 30 -nbins 3 -njk 8 $flag -out "$W/b" >/dev/null 2>&1
  if python3 - "$W/a" "$W/b" "$TOL" <<'PY'
import sys, numpy as np
a, b, tol = sys.argv[1], sys.argv[2], float(sys.argv[3]); ok = True
load = lambda p: np.atleast_2d(np.loadtxt([l for l in open(p) if not l.startswith('#')]))
A, B = load(a), load(b); hdr = [l for l in open(a) if l.startswith('#')][-1].lstrip('#').split()
for j in range(12, A.shape[1]):
    d = np.abs(A[:, j] - B[:, j]).max(); pk = np.abs(A[:, j]).max()
    rel = d / pk if pk > 0 else 0.0
    if rel > tol: print(f"\n    {hdr[j]}: max change {d:.3e} = {rel:.2e} of peak (tol {tol:g})", end=""); ok = False
sys.exit(0 if ok else 1)
PY
  then echo "PASS"; PASS=$((PASS+1)); else echo "  FAIL"; FAIL=$((FAIL+1)); fi
done
echo "Passed: $PASS   Failed: $FAIL"; [[ $FAIL -eq 0 ]]
