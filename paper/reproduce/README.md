# Reproducing the GRAMSCI v2 paper figures

This directory regenerates every data-driven figure in the paper from
compact, shipped data products. It is self-contained: no GRAMSCI binary,
GPU, or multi-GB survey catalogue is required.

## Quick start

```sh
pip install -r requirements.txt     # numpy + matplotlib only
./make_all_figures.sh               # writes figures/*.pdf and *.png
```

Each figure can also be made individually, e.g.
`cd scripts && python3 fig_benchmarks.py`.

## What reproduces what

| Paper figure | Script | Data product |
|---|---|---|
| Fig. 1 — chirality diagram | (none; pure TikZ in `gramsci_v2.tex`) | — |
| Fig. 2 — three-generation scaling | `fig_benchmarks.py` | `benchmarks.json`, `benchmarks2.json` |
| Fig. 3 — OpenMP thread scaling | `fig_benchmarks.py` | `benchmarks2.json` |
| Fig. 4 — density dependence | `fig_benchmarks.py` | `benchmarks2.json` |
| Fig. 5 — out-of-core overhead | `fig_benchmarks.py` | `benchmarks.json` |
| Fig. 6 — parity signal recovery | `fig_parity_recovery.py` | `parity_recovery.json` |
| Fig. 7 — DESI vs EZmock 3PCF (all configs + key) | `fig_desi_ezmock_3pcf.py` | `ezmock_3pcf_allconfig.csv`, `threepcf_config_scales.csv` |
| Fig. 8 — DESI vs EZmock 4PCF (full + connected) | `fig_desi_ezmock_4pcf.py` | `ezmock_4pcf_allcomb.csv` |
| Fig. 9 — Quijote fiducial 3PCF | `fig_quijote_fiducial.py` | `quijote_fiducial_3pcf.csv` |

The parity null-test statistic quoted in the text
(χ²/dof = 1.10 over 638 configurations, 66 Quijote realizations) is in
`data/parity_null_summary.txt`.

## Data products

`data/` holds the reduced inputs the figure scripts consume:

- **`benchmarks.json`, `benchmarks2.json`** — raw timing measurements
  (RTX 3090 Ti vs. 64-thread CPU), produced by
  `analysis/benchmarks/run_benchmarks{,2}.py` in the main repository.
- **`parity_recovery.json`** — the chiral-tetrahedra injection sweep
  (analytic recovery test), produced by `analysis/parity/parity_recovery.py`.
- **`ezmock_3pcf_allconfig.csv`, `ezmock_4pcf_allcomb.csv`** — the DESI DR1
  LRG measurements together with the EZmock ensemble mean and scatter, over
  all binned configurations, with NGC, SGC, and both redshift bins combined
  at the count level. The 4PCF file carries both the full
  ($\zeta^{(4)} = \mathrm{NNNN}/\mathrm{RRRR}$) and connected
  ($\zeta^{(4)}_{\rm conn} = \zeta^{(4)} - \zeta^{(4)}_{\rm disc}$) columns;
  the connected piece subtracts the disconnected term built from the
  internal 2PCF (corrected for the historical estimator sign offset, now
  fixed in the Fortran). Reduced from the per-realization measurement
  outputs (24 EZmocks for the 3PCF, 25 for the 4PCF).
- **`threepcf_config_scales.csv`** — the three sorted triangle side-length
  bin-centres ($r_1 \le r_2 \le r_3$) per 3PCF configuration, for the key
  panel of Fig. 7.
- **`quijote_fiducial_3pcf.csv`** — the 100-realization Quijote fiducial
  3PCF mean and scatter along the equilateral and isoceles cuts.

## The full measurement chain (not shipped here)

The compact products above are the last reduction step of a longer
pipeline. The upstream stages require the public catalogues, the GRAMSCI
binaries, and (for the GPU timings) an NVIDIA GPU, and produce ~70 MB of
per-realization text outputs. They are reproduced by the scripts in the
main repository, not bundled here:

1. **Catalogue preparation** — `analysis/quijote/quijote_to_gramsci.py`
   (Quijote FoF halos) and the DESI/EZmock converters in
   `analysis/`/`data` tooling turn public FITS/FoF catalogues into
   GRAMSCI input format.
2. **Measurement** — `bin/gramsci_gpu` computes the 2/3/4-point functions;
   driver scripts (`run_quijote_3pcf.sh`, the EZmock campaign scripts)
   loop over realizations.
3. **Reduction** — `scripts/reduce_to_package.py` stacks the
   per-realization outputs into the CSVs in `data/`. It is included here
   for transparency; running it requires the full measurement outputs on
   disk (paths are set for the authors' machine).

The reduced products in `data/` are sufficient to regenerate every
published figure exactly, which is the purpose of this directory.

## Data provenance

DESI DR1 LRG catalogues: <https://data.desi.lbl.gov/public/dr1> (CC BY 4.0).
Quijote simulations: <https://quijote-simulations.readthedocs.io>.
See the paper's acknowledgments and data-availability statement for the
required citations.
