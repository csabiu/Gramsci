#!/bin/bash
# Quijote measurement driver: per realization, convert FoF -> gramsci
# format, run GPU 3PCF + 4PCF, clean up. Mirrors the EZmock pipeline.
#
# Expects FoF catalogs already transferred (fetch_quijote.sh / Globus) under
#   $QUIJOTE_ROOT/<cosmo>/<realization>/groups_002/
#
# Usage: ./run_quijote_3pcf.sh <cosmo> <start> <end>
#   e.g. ./run_quijote_3pcf.sh fiducial 0 99

GRAMSCI=/home/csabiu/codes/Gramsci_gpu/Gramsci/bin/gramsci_gpu
PYCONV=/home/csabiu/codes/Gramsci_gpu/Gramsci/analysis/quijote/quijote_to_gramsci.py
QUIJOTE_ROOT=${QUIJOTE_ROOT:-/home/csabiu/data/quijote}
OUTDIR=${OUTDIR:-/home/csabiu/data/quijote/results}
TMPDIR=/tmp/quijote_run
SNAP=2          # z = 0.5
NDENS=1.5e-4      # (h/Mpc)^3, matched across cosmologies

COSMO=${1:?usage: run_quijote_3pcf.sh <cosmo> <start> <end>}
START=${2:-0}
END=${3:-0}

mkdir -p "$OUTDIR" "$TMPDIR"

for R in $(seq "$START" "$END"); do
    OUT3="$OUTDIR/${COSMO}_r${R}.3pcf"
    OUT4="$OUTDIR/${COSMO}_r${R}.4pcf"
    if [ -f "$OUT3" ] && [ -f "$OUT4" ]; then
        echo "=== $COSMO r$R already done ==="
        continue
    fi

    FOFDIR="$QUIJOTE_ROOT/$COSMO/$R/groups_$(printf '%03d' $SNAP)"
    if [ ! -d "$FOFDIR" ]; then
        echo "=== $COSMO r$R: FoF catalog missing ($FOFDIR), skipping ==="
        continue
    fi

    echo "=== $COSMO r$R at $(date '+%H:%M:%S') ==="
    python3 "$PYCONV" --fof-dir "$FOFDIR" --snapnum $SNAP --ndens $NDENS \
        --out "$TMPDIR/cat" --seed "$R" \
        || { echo "  convert failed"; continue; }

    "$GRAMSCI" -gal "$TMPDIR/cat.gal" -ran "$TMPDIR/cat.ran" \
        -rmin 10.0 -rmax 150.0 -nbins 20 -nmu 1 \
        -out "$OUT3" -3pcf > "$OUTDIR/${COSMO}_r${R}.3pcf.log" 2>&1 \
        || { echo "  3PCF failed"; rm -f "$OUT3"; }

    "$GRAMSCI" -gal "$TMPDIR/cat.gal" -ran "$TMPDIR/cat.ran" \
        -rmin 10.0 -rmax 65.0 -nbins 5 -nmu 1 \
        -out "$OUT4" -4pcf > "$OUTDIR/${COSMO}_r${R}.4pcf.log" 2>&1 \
        || { echo "  4PCF failed"; rm -f "$OUT4"; }

    rm -f "$TMPDIR"/cat.gal "$TMPDIR"/cat.ran
done

echo "done: $(ls "$OUTDIR/${COSMO}"_r*.3pcf 2>/dev/null | wc -l) 3PCF outputs for $COSMO"
