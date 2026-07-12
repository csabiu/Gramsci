# Changelog

Notable changes to GRAMSCI. This project accompanies Sabiu, Hoyle, Kim & Li,
*ApJS* **242**, 29 (2019), [arXiv:1901.00296](https://arxiv.org/abs/1901.00296).

## [2.3.0] — 2026-07-12

### Fixed
- **OpenCL: Kahan compensated accumulation in all four kernels.** The fp32
  per-column partials stopped registering increments once a column passed
  ~2/ε ≈ 1.7×10⁷ tuples (for unit weights) — a one-sided bias that silently
  undercounted the monotone RRR/RRRR denominator channels on
  production-size runs, inflating ζ. Each partial buffer now carries a
  same-size compensation buffer (`kadd` in `kernels.cl`); the host folds in
  `sum − comp` during its double reduction, keeping the error O(ε)
  independent of the tuple count. CPU-vs-OpenCL validation agreement
  tightened from ~10⁻⁴ to ~10⁻⁷ (fp32 ulp) as a side effect.
- **OpenCL: the tiled watchdog fallback now verifies every launch.**
  `cl_run_bucketed` bound the per-work-item completion flags but never read
  them back, and its first window is sized open-loop — so if a tiled launch
  itself exceeded the GPU watchdog (routine on the hour-scale 4PCF runs the
  fallback exists for), Apple's OpenCL returned partial results with no
  error and the dropped counts went undetected. The launcher now reads the
  flags after every window, periodically commits verified partials to
  double host accumulators (≤ ~5 s of GPU work ever at risk), and on
  truncation discards the uncommitted device partials, rewinds to the last
  committed window, and shrinks the window — truncation costs time, never
  counts. Set `GRAMSCI_CL_FORCE_TILED=1` to exercise the tiled path
  regardless of the single-launch outcome (used by `validate.sh` testing).
- OpenCL: the `CL_DEVICE_DOUBLE_FP_CONFIG` query no longer aborts
  initialisation on pre-1.2 devices that answer it with an error instead
  of 0.
- OpenACC driver: the graph-build ETA printed total process CPU time
  instead of the probe duration.

- `make clean` in `tests/` removes the `tmp_*`/`bench_*` scratch files
  (the recipe referenced an undefined variable and was a no-op);
  `src/Makefile` no longer installs a stale binary when the link step
  fails.

### Removed
- **The MPI/domain-decomposition path (`-cut`, `-mpi`,
  `bin/domain_decomposition`).** It was broken end-to-end: results were
  never reduced across ranks (only the master's partial counts were
  written), and the `-cut` reader took the `gal` column as the halo buffer
  flag, so region overlaps were double-counted. All three drivers are now
  single-process (OpenMP/GPU) only.

## [2.2.0] — 2026-07-06

### Added
- **Anisotropic (RSD-aware) disconnected-4PCF subtraction.** In redshift
  space the Gaussian term of the 4PCF is not the product of isotropic ξ's:
  the line-of-sight angles of the two edges in each complementary pairing
  co-vary because both edges are rigidly attached to the same tetrahedron.
  By the Legendre addition theorem the orientation-averaged term is
  `ξ₀ξ₀ + ξ₂ξ₂L₂(cosθ)/5 + ξ₄ξ₄L₄(cosθ)/9`, with θ the opposite-edge angle
  fixed by the six side lengths. The internal 2PCF pass now also measures
  ξ₂(r) and ξ₄(r) — pair Legendre sums are accumulated at graph
  construction, where the full-precision pair μ exists (plane-parallel z
  line of sight in `-box` mode; midpoint line of sight in survey mode with
  `-nmu > 1`) — and `zeta_disc`/`zeta_conn` include the multipole terms.
  For real-space data ξ₂ ≈ ξ₄ ≈ 0 and the estimator reduces exactly to the
  previous isotropic subtraction; the parity-odd channel has no
  disconnected term and is unchanged. Applies to all three builds.
  Verified against an independent pair-count + addition-theorem
  computation on a z-squashed clustered box (agreement ≤ 4×10⁻⁷ over all
  configurations; correction ~10% of `zeta_disc` on that field), and the
  1/5, 1/9 coefficients against direct Monte-Carlo orientation averaging.
- `tests/run_correlation_tests.py`: anisotropic disconnected-subtraction
  regression test (`generate_aniso_box` field).

### Fixed
- Periodic-mode graph fill now clamps the distance-bin index: the neighbor
  filter compares squared distances, and sqrt rounding could land a pair
  exactly on `rmax`, producing an out-of-range bin (found by adversarial
  review of 2.1.0).

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
