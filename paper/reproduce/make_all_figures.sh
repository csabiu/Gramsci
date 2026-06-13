#!/bin/bash
# Regenerate every data-driven figure in the GRAMSCI v2 paper from the
# shipped reduced data products.  Outputs land in ./figures/.
#
# Requires only Python 3 with numpy + matplotlib (see requirements.txt).
# No GRAMSCI binary, GPU, or raw catalogues are needed: the figures are
# regenerated from the compact CSV/JSON products in ./data/.
set -e
cd "$(dirname "$0")/scripts"

echo "=== Figs 2-5: benchmarks (scaling, threads, density, out-of-core) ==="
python3 fig_benchmarks.py
echo "=== Fig 6: parity-estimator signal recovery ==="
python3 fig_parity_recovery.py
echo "=== Fig 7: DESI vs EZmock, 3PCF ==="
python3 fig_desi_ezmock_3pcf.py
echo "=== Fig 8: DESI vs EZmock, 4PCF ==="
python3 fig_desi_ezmock_4pcf.py
echo "=== Fig 9: Quijote fiducial 3PCF ==="
python3 fig_quijote_fiducial.py

echo ""
echo "All figures written to figures/ :"
ls -1 ../figures/*.pdf | sed 's#.*/#  #'
echo ""
echo "Figure 1 (chirality diagram) is pure TikZ in gramsci_v2.tex and is"
echo "produced by LaTeX directly; no script is required."
