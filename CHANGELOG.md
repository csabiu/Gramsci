# Changelog

Notable changes to GRAMSCI. This project accompanies Sabiu, Hoyle, Kim & Li,
*ApJS* **242**, 29 (2019), [arXiv:1901.00296](https://arxiv.org/abs/1901.00296).

## [2.1.0] — 2026-07-04

### Added
- **Periodic-box mode with analytic randoms** (`-box L`, or `-box Lx,Ly,Lz`;
  Python: `box=`). For regular geometries such as simulation boxes, pair
  separations use the **minimum-image convention** and — when no `-ran`
  catalogue is given — the `RR`/`RRR`/`RRRR` counts are computed
  **analytically** instead of from random points: shell volumes for the 2PCF,
  the `8π² r₁r₂r₃` triangle kernel for the 3PCF, and a semi-analytic
  tetrahedron-kernel quadrature (analytic azimuthal integral, exact
  breakpoint-split Gauss–Legendre elsewhere) for the 4PCF, including the
  orbit multiplicities of the S₄-canonicalized bin 6-tuples. Estimators
  switch to the natural forms with the lower-order hierarchy terms
  subtracted internally (`ξ = DD/RR − 1`;
  `ζ = DDD/RRR − ξ(r₁) − ξ(r₂) − ξ(r₃) − 1`; the 4PCF additionally
  subtracts the six edge `ξ` and four face `ζ₃` terms via an internal 3PCF
  pass, so `zeta`/`zeta_disc`/`zeta_conn` keep their usual meanings; the
  parity-odd channel needs no subtraction and analytic `RRRR_odd = 0`).
  `-box` with `-ran` is also supported (periodic distances, catalogue
  estimators). Requires `rmax < L/2` (2PCF) or `rmax ≤ L/4` (3/4PCF), and
  `-nmu 1`. Works in all three builds (CPU, OpenCL, OpenACC) — wrapping
  happens at graph construction and the GPU kernels only see bin indices.
  Validated against Monte-Carlo counting and random-catalogue runs: the
  analytic `RR`/`RRR`/`RRRR` normalizations agree to `0.9995 ± 0.0008`,
  `0.9990 ± 0.0022` and `0.9996 ± 0.0033` respectively.
- `tests/run_correlation_tests.py`: periodic-box test suite (exact RR check,
  2/3/4PCF null tests on a uniform box, signal localization in analytic
  mode, catalogue-vs-analytic RR cross-check, invalid-flag rejection).

### Fixed
- `create_graph` no longer reads past the end of the coordinate array when
  the catalogue has fewer than 1000 points (the drivers probe the first 999
  hubs for the timing estimate; the loop is now clamped).

## [2.0.0] — 2026-06-18

### Added
- **Portable OpenCL GPU backend** (`src_opencl/` → `bin/gramsci_cl`). Offloads
  the 3PCF, equilateral 3PCF, 4PCF and parity-4PCF queries to any OpenCL 1.2
  device — including **Apple Silicon** and Intel Macs, and Linux AMD/NVIDIA/Intel
  GPUs — with the same CLI and output formats as the CPU build. Driven from
  Fortran via a hand-written `ISO_C_BINDING` interface (no third-party
  dependency). Validated against the CPU reference to ~10⁻⁶ relative.
- **Python tutorial notebook** `example/quickstart.ipynb`: catalogue → 2/3/4PCF
  → plots.
- Top-level `make gpu` (OpenCL) and `make cuda` (NVIDIA OpenACC) build targets.
- `src_opencl/validate.sh` — CPU-vs-OpenCL regression check.

### Changed
- **README rewritten in Markdown** (`README.md`): hero figure, one-command
  install, Python and command-line quickstarts, output reference, and citation.

### Notes
- The OpenCL backend runs in **single precision** (Apple Silicon GPUs report no
  fp64). For publication-grade double precision, use the CPU build (`src/`) or
  the NVIDIA OpenACC build (`src_gpu/`).
- On a small integrated GPU, the largest queries automatically fall back to a
  slower — but always correct — tiled path because of the OS GPU watchdog; see
  `src_opencl/README.md` for details.

## [1.0.0] — 2019

- Initial public release: CPU 2/3/4-point and parity-4PCF, the KD-tree
  graph-database engine, and the NVIDIA OpenACC GPU build (`src_gpu/`).
  Method paper: Sabiu, Hoyle, Kim & Li, *ApJS* **242**, 29 (2019).
