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
cpu = np.loadtxt(sys.argv[1], ndmin=2)   # all header lines start with '#'
gpu = np.loadtxt(sys.argv[2], ndmin=2)
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

# Jackknife comparison: run one mode with -njk 8 on both binaries and compare
# the main output plus every sidecar (.jk, .jkerr, .jkcov, .jkcov_odd for
# -4pcfp) with the same per-column relative-to-peak test.  The jackknife
# forces the tiled launch path, whose per-window commits keep the fp32 CAS
# touching sums accurate; pass a GRAMSCI_CL_TARGET_SEC value as the third
# argument to force many tiny windows and stress the commit/re-zero cycle.
run_jk_test() {
  local mode="$1"
  local flag="$2"
  local tsec="${3:-}"
  local tmpd; tmpd=$(mktemp -d)
  local cpu_out="$tmpd/cpu.out"
  local gpu_out="$tmpd/gpu.out"
  echo -n "  [$mode] ... "

  "$CPU_BIN" -gal "$GAL" -ran "$RAN" -rmin 10 -rmax 30 -nbins 4 $flag \
             -njk 8 -out "$cpu_out" >/dev/null 2>&1
  if [[ -n "$tsec" ]]; then
    GRAMSCI_CL_TARGET_SEC="$tsec" \
    "$CL_BIN" -gal "$GAL" -ran "$RAN" -rmin 10 -rmax 30 -nbins 4 $flag \
              -njk 8 -out "$gpu_out" >/dev/null 2>&1
  else
    "$CL_BIN" -gal "$GAL" -ran "$RAN" -rmin 10 -rmax 30 -nbins 4 $flag \
              -njk 8 -out "$gpu_out" >/dev/null 2>&1
  fi

  local ok=1 worst_all=""
  if ! cmp -s "$cpu_out.jkgal" "$gpu_out.jkgal"; then
    echo -n "REGION LABELS DIFFER "
    ok=0
  fi
  local sufs="_MAIN_ .jk .jkerr .jkcov"
  [[ "$flag" == *"-4pcfp"* ]] && sufs="$sufs .jkcov_odd"
  for suf in $sufs; do
    [[ "$suf" == "_MAIN_" ]] && suf=""
    if ! python3 - "$cpu_out$suf" "$gpu_out$suf" "$TOL" <<'PYEOF'
import sys, numpy as np
cpu = np.loadtxt(sys.argv[1], comments='#', ndmin=2)
gpu = np.loadtxt(sys.argv[2], comments='#', ndmin=2)
tol = float(sys.argv[3])
if cpu.shape != gpu.shape:
    print(f"SHAPE MISMATCH: CPU {cpu.shape} vs OpenCL {gpu.shape}"); sys.exit(1)
worst = 0.0
for c in range(cpu.shape[1]):
    scale = np.nanmax(np.abs(cpu[:, c]))
    if not np.isfinite(scale) or scale == 0:
        continue
    d = np.nanmax(np.abs(np.nan_to_num(cpu[:, c]) - np.nan_to_num(gpu[:, c]))) / scale
    worst = max(worst, d)
if worst > tol:
    print(f"MISMATCH in {sys.argv[1].split('/')[-1]} (worst rel = {worst:.3e})"); sys.exit(1)
PYEOF
    then
      echo -n "(in ${suf:-main}) "
      ok=0
    fi
  done

  if [[ $ok -eq 1 ]]; then echo "PASS (tol=$TOL)"; PASS=$((PASS+1))
  else echo "FAIL"; FAIL=$((FAIL+1)); fi
  rm -rf "$tmpd"
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
echo "Jackknife (-njk 8, tiled with per-window commits)..."
run_jk_test "3PCF jk"   "-3pcf"
run_jk_test "4PCF jk"   "-4pcf"
run_jk_test "4PCFp jk"  "-4pcfp"
echo ""
echo "Jackknife commit-cycle stress (GRAMSCI_CL_TARGET_SEC=0.02)..."
run_jk_test "4PCFp jk stress" "-4pcfp" 0.02
echo ""
echo "Results: $PASS passed, $FAIL failed"
[[ $FAIL -eq 0 ]] && echo "All tests PASSED." && exit 0
echo "Some tests FAILED." && exit 1
