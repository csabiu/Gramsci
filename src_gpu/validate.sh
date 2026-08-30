#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# validate.sh — compare CPU and GPU GRAMSCI output for 2PCF, 3PCF, 4PCF,
# and the delete-one jackknife (-njk) outputs of the 3PCF/4PCF/4PCFp,
# in both single-pass and forced-chunked (GRAMSCI_GPU_WIN_EDGES) modes.
#
# Usage: cd tests && bash ../src_gpu/validate.sh [tolerance]
#
# The script runs each mode with both binaries and diffs the numerical output
# using Python/numpy (tolerates floating-point reordering from GPU atomics).
#
# Requirements:
#   - ../bin/gramsci and ../bin/gramsci_gpu both built
#   - test data: test.gal and test.ran (provided in tests/)
#   - Python 3 with numpy
#
# Exit code 0 = all tests passed, non-zero = one or more failed.
# ---------------------------------------------------------------------------
set -euo pipefail

BINDIR="$(dirname "$0")/../bin"
TESTDIR="$(dirname "$0")/../tests"
TOL="${1:-1e-8}"      # absolute tolerance for numerical comparison

CPU_BIN="$BINDIR/gramsci"
GPU_BIN="$BINDIR/gramsci_gpu"
GAL="$TESTDIR/test.gal"
RAN="$TESTDIR/test.ran"
if [[ ! -f "$GAL" ]]; then
  GAL="$(dirname "$0")/../example/test.gal"
  RAN="$(dirname "$0")/../example/test.ran"
fi

if [[ ! -x "$CPU_BIN" ]]; then echo "ERROR: $CPU_BIN not found"; exit 1; fi
if [[ ! -x "$GPU_BIN" ]]; then echo "ERROR: $GPU_BIN not found"; exit 1; fi
if [[ ! -f "$GAL" ]];     then echo "ERROR: $GAL not found";     exit 1; fi
if [[ ! -f "$RAN" ]];     then echo "ERROR: $RAN not found";     exit 1; fi

PASS=0
FAIL=0

run_test() {
  local mode="$1"
  local flag="$2"
  local cpu_out; cpu_out=$(mktemp)
  local gpu_out; gpu_out=$(mktemp)

  echo -n "  [$mode] ... "

  "$CPU_BIN" -gal "$GAL" -ran "$RAN" -rmin 10 -rmax 30 -nbins 3 $flag \
             -out "$cpu_out" 2>/dev/null
  "$GPU_BIN" -gal "$GAL" -ran "$RAN" -rmin 10 -rmax 30 -nbins 3 $flag \
             -out "$gpu_out" 2>/dev/null

  if python3 - <<PYEOF
import sys, numpy as np
tol = float("$TOL")
# '#' marks the header; 'r' kept for outputs written before the header fix.
cpu = np.loadtxt("$cpu_out", comments=('#', 'r'), ndmin=2)
gpu = np.loadtxt("$gpu_out", comments=('#', 'r'), ndmin=2)
if cpu.shape != gpu.shape:
    print(f"SHAPE MISMATCH: CPU {cpu.shape} vs GPU {gpu.shape}")
    sys.exit(1)
diff = np.abs(cpu - gpu)
bad = diff > tol
if bad.any():
    print(f"NUMERICAL MISMATCH (max abs diff = {diff.max():.3e}, tol = {tol:.3e})")
    sys.exit(1)
PYEOF
  then
    echo "PASS (tol=$TOL)"
    PASS=$((PASS + 1))
  else
    echo "FAIL"
    FAIL=$((FAIL + 1))
  fi

  rm -f "$cpu_out" "$gpu_out"
}

# Compare one pair of numeric files with the absolute tolerance.
compare_file() {
  python3 - "$1" "$2" <<CMPEOF
import sys, numpy as np
tol = float("$TOL")
cpu = np.loadtxt(sys.argv[1], comments=('#', 'r'), ndmin=2)
gpu = np.loadtxt(sys.argv[2], comments=('#', 'r'), ndmin=2)
if cpu.shape != gpu.shape:
    print(f"SHAPE MISMATCH: CPU {cpu.shape} vs GPU {gpu.shape}")
    sys.exit(1)
diff = np.abs(cpu - gpu)
if (diff > tol).any():
    print(f"NUMERICAL MISMATCH (max abs diff = {diff.max():.3e}, tol = {tol:.3e})")
    sys.exit(1)
CMPEOF
}

# Jackknife comparison: run one mode with -njk on both binaries and compare
# the main output plus every jackknife sidecar (.jk realisations, .jkerr,
# .jkcov, and .jkcov_odd for the parity 4PCF).  The internal angular
# partition is deterministic (quantiles of the randoms), so both binaries
# assign identical region labels — checked byte-for-byte via .jkgal.
# An optional third argument forces the chunked kernels by setting
# GRAMSCI_GPU_WIN_EDGES.
run_jk_test() {
  local mode="$1"
  local flag="$2"
  local winedges="${3:-}"
  local tmpd; tmpd=$(mktemp -d)
  local cpu_out="$tmpd/cpu.out"
  local gpu_out="$tmpd/gpu.out"

  echo -n "  [$mode] ... "

  "$CPU_BIN" -gal "$GAL" -ran "$RAN" -rmin 10 -rmax 30 -nbins 3 $flag \
             -njk 8 -out "$cpu_out" >/dev/null 2>&1
  if [[ -n "$winedges" ]]; then
    GRAMSCI_GPU_WIN_EDGES="$winedges" \
    "$GPU_BIN" -gal "$GAL" -ran "$RAN" -rmin 10 -rmax 30 -nbins 3 $flag \
               -njk 8 -out "$gpu_out" >/dev/null 2>&1
  else
    "$GPU_BIN" -gal "$GAL" -ran "$RAN" -rmin 10 -rmax 30 -nbins 3 $flag \
               -njk 8 -out "$gpu_out" >/dev/null 2>&1
  fi

  local ok=1
  if ! cmp -s "$cpu_out.jkgal" "$gpu_out.jkgal"; then
    echo -n "REGION LABELS DIFFER "
    ok=0
  fi
  local sufs="_MAIN_ .jk .jkerr .jkcov"
  [[ "$flag" == *"-4pcfp"* ]] && sufs="$sufs .jkcov_odd"
  for suf in $sufs; do
    [[ "$suf" == "_MAIN_" ]] && suf=""
    if ! compare_file "$cpu_out$suf" "$gpu_out$suf"; then
      echo -n "(in ${suf:-main output}) "
      ok=0
    fi
  done

  if [[ $ok -eq 1 ]]; then
    echo "PASS (tol=$TOL)"
    PASS=$((PASS + 1))
  else
    echo "FAIL"
    FAIL=$((FAIL + 1))
  fi
  rm -rf "$tmpd"
}

echo "=== GRAMSCI CPU vs GPU validation (tol=$TOL) ==="
echo ""
echo "Running tests..."

run_test "2PCF"  "-2pcf"
run_test "3PCF"  "-3pcf"
run_test "4PCF"  "-4pcf"
run_test "4PCFp" "-4pcfp"

echo ""
echo "Jackknife (-njk 8, single-pass)..."
run_jk_test "3PCF jk"   "-3pcf"
run_jk_test "4PCF jk"   "-4pcf"
run_jk_test "4PCFp jk"  "-4pcfp"

echo ""
echo "Jackknife (-njk 8, forced-chunked kernels)..."
run_jk_test "4PCF jk chunked"  "-4pcf"  5000000
run_jk_test "4PCFp jk chunked" "-4pcfp" 5000000

echo ""
echo "Results: $PASS passed, $FAIL failed"
[[ $FAIL -eq 0 ]] && echo "All tests PASSED." && exit 0
echo "Some tests FAILED." && exit 1
