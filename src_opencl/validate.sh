#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# validate.sh — compare CPU vs OpenCL GRAMSCI output for 2PCF/3PCF/4PCF/4PCFp.
#
# Usage: cd tests && bash ../src_opencl/validate.sh [rel_tol]
#
# The OpenCL backend runs in single precision (Apple Silicon OpenCL has no
# fp64), so this uses a per-column relative-to-peak tolerance (default 2e-3)
# rather than the tight absolute tolerance of the double-precision GPU build.
#
# Requirements:
#   - ../bin/gramsci and ../bin/gramsci_cl both built
#   - test data: test.gal and test.ran (in tests/ or example/)
#   - Python 3 with numpy
#
# Exit 0 = all passed, non-zero = one or more failed.
# ---------------------------------------------------------------------------
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
BINDIR="$HERE/../bin"
TOL="${1:-2e-3}"      # per-column relative-to-peak tolerance (single precision)

CPU_BIN="$BINDIR/gramsci"
CL_BIN="$BINDIR/gramsci_cl"

# Locate test data (tests/ preferred, else example/)
GAL=""; RAN=""
for d in "$HERE/../tests" "$HERE/../example" "."; do
  if [[ -f "$d/test.gal" && -f "$d/test.ran" ]]; then GAL="$d/test.gal"; RAN="$d/test.ran"; break; fi
done

[[ -x "$CPU_BIN" ]] || { echo "ERROR: $CPU_BIN not found (build with: cd src && make)"; exit 1; }
[[ -x "$CL_BIN"  ]] || { echo "ERROR: $CL_BIN not found (build with: cd src_opencl && make)"; exit 1; }
[[ -n "$GAL"     ]] || { echo "ERROR: test.gal/test.ran not found in tests/ or example/"; exit 1; }

PASS=0; FAIL=0

run_test() {
  local mode="$1" flag="$2"
  local cpu_out gpu_out
  cpu_out=$(mktemp); gpu_out=$(mktemp)
  echo -n "  [$mode] ... "

  "$CPU_BIN" -gal "$GAL" -ran "$RAN" -rmin 5 -rmax 30 -nbins 4 $flag -out "$cpu_out" >/dev/null 2>&1
  "$CL_BIN"  -gal "$GAL" -ran "$RAN" -rmin 5 -rmax 30 -nbins 4 $flag -out "$gpu_out" >/dev/null 2>&1

  if python3 - "$cpu_out" "$gpu_out" "$TOL" <<'PYEOF'
import sys, numpy as np
cpu = np.loadtxt(sys.argv[1], comments='r', ndmin=2)
gpu = np.loadtxt(sys.argv[2], comments='r', ndmin=2)
tol = float(sys.argv[3])
if cpu.shape != gpu.shape:
    print(f"SHAPE MISMATCH: CPU {cpu.shape} vs OpenCL {gpu.shape}"); sys.exit(1)
# per-column relative-to-peak: |dcol| <= tol * max|cpu col|
worst = 0.0
for c in range(cpu.shape[1]):
    col = cpu[:, c]
    scale = np.nanmax(np.abs(col))
    if not np.isfinite(scale) or scale == 0:  # all-zero / NaN-ratio column
        continue
    d = np.nanmax(np.abs(np.nan_to_num(cpu[:, c]) - np.nan_to_num(gpu[:, c]))) / scale
    worst = max(worst, d)
if worst > tol:
    print(f"MISMATCH (worst per-column rel = {worst:.3e}, tol = {tol:.3e})"); sys.exit(1)
print(f"max rel = {worst:.3e}", end="  ")
PYEOF
  then echo "PASS (tol=$TOL)"; PASS=$((PASS+1))
  else echo "FAIL"; FAIL=$((FAIL+1)); fi
  rm -f "$cpu_out" "$gpu_out"
}

echo "=== GRAMSCI CPU vs OpenCL validation (single precision, tol=$TOL) ==="
echo "data: $GAL / $RAN"
echo ""
run_test "2PCF"  "-2pcf"
run_test "3PCF"  "-3pcf"
run_test "equi"  "-equi"
run_test "4PCF"  "-4pcf"
run_test "4PCFp" "-4pcfp"
echo ""
echo "Results: $PASS passed, $FAIL failed"
[[ $FAIL -eq 0 ]] && echo "All tests PASSED." && exit 0
echo "Some tests FAILED." && exit 1
