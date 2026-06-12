# Quijote validation & cosmology-dependence program

Measurement program for the GRAMSCI v2 paper (Sections "Validation with
Quijote" and "Cosmological dependence").

## Goal

1. **Validation**: 3PCF of fiducial-cosmology FoF halos vs. perturbation
   theory at large scales; GPU vs. CPU agreement; timing at box scale.
2. **Cosmology dependence**: fractional response of the 3PCF and 4PCF to
   single-parameter variations, with error bars from fiducial realizations.

## Data

Quijote FoF halo catalogs (Villaescusa-Navarro et al. 2020),
1 (Gpc/h)^3 boxes, 512^3 particles, snapshot z=0.5 (snapnum=2):

| Cosmology  | Realizations | Purpose                |
|------------|--------------|------------------------|
| fiducial   | 100          | mean + error bars      |
| Om_p, Om_m | 50 each      | dζ/dΩm                 |
| s8_p, s8_m | 50 each      | dζ/dσ8                 |
| h_p, h_m   | 50 each      | dζ/dh                  |

Halo selection: number density cut n = 1.5×10⁻⁴ (h/Mpc)³ (most massive
halos first) → 150,000 halos per box, matched across cosmologies so the
response isolates clustering, not abundance. Randoms: uniform in the box,
1:1.

## Access

Bulk Quijote data is served via **Globus** (endpoint: "Quijote_simulations",
see https://quijote-simulations.readthedocs.io/en/latest/access.html).
One-time setup:

```sh
pip install globus-cli
globus login                  # browser auth
./fetch_quijote.sh            # transfers the FoF catalogs listed above
```

Each FoF catalog is ~40-100 MB; the full program is ~30 GB. With the
download-process-delete driver the peak disk use stays under ~2 GB.

## Pipeline

```
fetch_quijote.sh        Globus transfer of FoF catalogs (or per-mock in driver)
quijote_to_gramsci.py   FoF binary -> .gal (density cut) + uniform randoms
run_quijote_3pcf.sh     measurement driver: 3PCF (20 bins to 150 Mpc/h)
                        and 4PCF (rmax=65, 3 bins) per realization, GPU
analyze_responses.py    stacks measurements, computes fractional responses
                        with fiducial error bars, makes paper figures
```

## Compute budget (RTX 3090 Ti)

200k halos + 200k randoms per box; estimated per realization:
- graph build ~20 s; 3PCF rmax=150 query: ~2-3 min (no chunking needed,
  ~1.5e9 edges ≈ 7.5 GB)
- 4PCF rmax=65: ~1 min

400 realizations × ~4 min ≈ **~28 GPU-hours** — about two days unattended.
