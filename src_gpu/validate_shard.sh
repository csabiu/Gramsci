#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# validate_shard.sh — check the multi-GPU hub sharding of the 4PCF kernels:
# -shard k/n runs + -merge n must reproduce the unsharded output, for -4pcf
# and -4pcfp with -njk, in single-pass and forced-chunked (GRAMSCI_GPU_WIN_EDGES)
# mode.  The merge is run with no visible device (CUDA_VISIBLE_DEVICES=).
#
# Usage: cd tests && bash ../src_gpu/validate_shard.sh [nshards] [tolerance]
#
# Requirements: ../bin/gramsci_gpu (override with GPU_BIN=...), test data as
# for validate.sh, Python 3 with numpy.  Exit 0 = all passed.
# ---------------------------------------------------------------------------
set -euo pipefail
BINDIR="$(dirname "$0")/../bin"
TESTDIR="$(dirname "$0")/../tests"
NSHARD="${1:-3}"
TOL="${2:-2e-6}"      # per-column relative-to-peak; the text outputs carry 7 digits
GPU_BIN="${GPU_BIN:-$BINDIR/gramsci_gpu}"
GAL="$TESTDIR/test.gal"
RAN="$TESTDIR/test.ran"
if [[ ! -f "$GAL" ]]; then
  GAL="$(dirname "$0")/../example/test.gal"
  RAN="$(dirname "$0")/../example/test.ran"
fi
[[ -x "$GPU_BIN" ]] || { echo "ERROR: $GPU_BIN not found"; exit 1; }
[[ -f "$GAL" && -f "$RAN" ]] || { echo "ERROR: test data not found"; exit 1; }
PARAMS="-rmin 10 -rmax 30 -nbins 3 -njk 8"
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
PASS=0
FAIL=0

compare() {   # compare <ref> <merged> <mode flag>
  local sufs=" .jk .jkerr .jkcov"
  [[ "$3" == "-4pcfp" ]] && sufs="$sufs .jkcov_odd"
  python3 - "$1" "$2" "$TOL" $sufs <<'PYEOF'
import sys, numpy as np
ref, mrg, tol = sys.argv[1], sys.argv[2], float(sys.argv[3])
ok = True
for suf in [''] + sys.argv[4:]:
    load = lambda p: np.atleast_2d(np.loadtxt([l for l in open(p) if not l.startswith('#')]))
    a, b = load(mrg + suf), load(ref + suf)
    if a.shape != b.shape:
        print(f"    {suf or '.out':<11} SHAPE MISMATCH {a.shape} vs {b.shape}"); ok = False; continue
    d = np.abs(a - b)
    # per-column relative-to-peak (as validate.sh for OpenCL): cell-wise relative
    # error is meaningless for near-zero covariance entries
    peak = np.abs(b).max(axis=0); pk = peak > 0
    rel = (d.max(axis=0)[pk] / peak[pk]).max() if pk.any() else 0.0
    print(f"    {suf or '.out':<11} differing cells {int((d > 0).sum()):>5}/{d.size:<6} max rel-to-peak {rel:.2e}")
    ok &= rel <= tol
sys.exit(0 if ok else 1)
PYEOF
}

run_case() {   # run_case <label> <mode flag> [GRAMSCI_GPU_WIN_EDGES]
  local label="$1" flag="$2" win="${3:-}"
  local ref="$WORK/ref_${label// /_}" mrg="$WORK/mrg_${label// /_}"
  echo "$label ($NSHARD shards${win:+, forced-chunked})..."
  "$GPU_BIN" -gal "$GAL" -ran "$RAN" $PARAMS $flag -out "$ref" >/dev/null 2>&1
  for k in $(seq 1 "$NSHARD"); do
    GRAMSCI_GPU_WIN_EDGES=$win "$GPU_BIN" -gal "$GAL" -ran "$RAN" $PARAMS $flag \
      -shard "$k/$NSHARD" -out "$mrg" >/dev/null 2>&1
  done
  CUDA_VISIBLE_DEVICES= "$GPU_BIN" -gal "$GAL" -ran "$RAN" $PARAMS $flag \
      -merge "$NSHARD" -out "$mrg" >/dev/null 2>&1
  if compare "$ref" "$mrg" "$flag"; then
    echo "  PASS"; PASS=$((PASS + 1))
  else
    echo "  FAIL"; FAIL=$((FAIL + 1))
  fi
}

run_case "4PCF"  "-4pcf"
run_case "4PCFp" "-4pcfp"
run_case "4PCFp chunked" "-4pcfp" 3000000

echo
echo "Passed: $PASS   Failed: $FAIL"
[[ $FAIL -eq 0 ]]
