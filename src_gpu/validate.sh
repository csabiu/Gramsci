#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# validate.sh — compare CPU and GPU GRAMSCI output for 2PCF, 3PCF, 4PCF.
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

echo "=== GRAMSCI CPU vs GPU validation (tol=$TOL) ==="
echo ""
echo "Running tests..."

run_test "2PCF"  "-2pcf"
run_test "3PCF"  "-3pcf"
run_test "4PCF"  "-4pcf"
run_test "4PCFp" "-4pcfp"

echo ""
echo "Results: $PASS passed, $FAIL failed"
[[ $FAIL -eq 0 ]] && echo "All tests PASSED." && exit 0
echo "Some tests FAILED." && exit 1
