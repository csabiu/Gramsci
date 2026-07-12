# Changelog

Notable changes to GRAMSCI. This project accompanies Sabiu, Hoyle, Kim & Li,
*ApJS* **242**, 29 (2019), [arXiv:1901.00296](https://arxiv.org/abs/1901.00296).

## [2.4.0] — 2026-07-12

### Added
- **Combined query modes.** Several of `-2pcf`/`-3pcf`/`-equi`/`-4pcf`/
  `-4pcfp` may now be requested in one run: the KD-tree and neighbor graph
  (usually the dominant cost) are built once, every query starts from
  freshly zeroed accumulators, and with more than one mode each result is
  written to `<out>.<mode>` (a single mode keeps the exact `-out` name).
  This replaces the 2.3.0 rule that rejected combined flags — that rule
  fixed the cross-contamination but forced one full graph build per
  statistic. Combined outputs are verified byte-identical to the
  corresponding single-mode runs (CPU regression test added).
- **`-version` flag** (`gramsci 2.4.0 (cpu build)` — the OpenCL and
  OpenACC binaries report their backend), and a **provenance header** in
  every output file: version/build, timestamp, the exact command line,
  catalogue sizes, and box geometry, as `#` comment lines. The column-name
  header is now also a `#` comment (and the equilateral output gains one),
  so `np.loadtxt` reads every output with no `skiprows` needed — existing
  `skiprows=1` calls keep working since the skipped line is a comment.
- **Continuous integration** (GitHub Actions): every push/PR to master
  builds the CPU binary and runs the full physics regression suite on
  Linux and macOS, link-checks the OpenCL backend on both, and
  syntax-checks the OpenACC backend with `gfortran -fopenacc`. The
  GPU-dependent numerical validation (`src_opencl/validate.sh`) still
  needs local hardware.

## [2.3.0] — 2026-07-12

### Fixed
- **RSD equilateral 3PCF (`-equi` with `-nmu > 1`) was fully broken**: the
  RSD branch never accumulated the random-triple denominator, so every
  output row was `zeta = NNN/0` = Inf. It now accumulates `N3` per mu bin
  exactly like the all-configurations RSD query; the mu-summed counts
  reproduce the isotropic run to machine precision.
- **Query modes are now mutually exclusive** (exactly one of
  `-2pcf | -3pcf | -equi | -4pcf | -4pcfp`). They all accumulate into the
  shared N2/N3 arrays without resetting between queries and write to the
  same output file, so combining flags silently cross-contaminated the
  counts and clobbered the earlier output.
- **Catalogue reading is now line-based and validated.** The old
  record-based reader always consumed four values, so a 3-column file
  (weight omitted — a format the docs advertise) consumed two lines per
  point and silently left half the arrays uninitialized after a one-line
  warning; a `#` header made the line counter stop at zero points. Files
  may now be `x y z w` or `x y z` (weight = 1), blank/`#` lines are
  skipped, and any malformed line is a hard error with its line number.
- **Non-periodic graph build hardened against edge-of-bin rounding**: the
  fill pass now uses exactly the same squared-distance filter as the
  sizing pass (a sqrt-based filter could disagree at values rounding onto
  `rmax²` and overflow the per-node arrays), and the radial bin index is
  clamped into `[1, nbins]` as the periodic pass already did — an
  unclamped index could reach `nbins+1` and corrupt memory (~1e-16 per
  pair, but production runs do 10¹¹⁺ pairs).
- **`find_normal` no longer takes the square root of a negative number**
  for near-line-of-sight edge pairs (previously `int(floor(NaN))` —
  undefined behaviour); the argument is clamped and the computation
  promoted from single to double precision.
- **4PCF Burnside bound computed in int64**: `nbins**6` overflowed default
  integers at `nbins >= 36`, turning the allocation size negative; it now
  fails with a clear message if the config count exceeds the int32 range.
- **CLI parsing**: options longer than 6 characters are no longer silently
  truncated (`-nbinsfoo` used to parse as `-nbins`); a trailing option
  with no value and unparseable numeric values are clean errors; added
  validation for `-rmax <= -rmin` and `-log` with `-rmin <= 0`.
- All Fortran error paths exit with a **nonzero status** (plain `stop`
  returns 0, which made failures indistinguishable from success for
  callers); the Python wrapper accordingly treats any nonzero exit as an
  error and surfaces the captured Fortran message instead of a generic
  "did not produce output".
- **Python wrapper**: `compute_3pcf(..., nmu>=2)` no longer silently
  misparses the 11-column RSD output (mu-bin edges were returned as
  counts and raw counts as zeta); `TwoPCFResult` exposes the mu-bin
  columns; calling `compute_*` with neither `randoms_pos` nor `box`
  raises immediately instead of returning Inf/NaN; `randoms_weights`
  defaults to 1 when only `randoms_pos` is given; relative `GRAMSCI_BIN`
  paths are resolved against the caller's cwd (the subprocess runs in a
  temp dir).
- Non-short-circuit `.and.` conditions that indexed one-past-the-end
  (insertion sorts, 4PCF merge-walk seek) restructured — they aborted
  under `-fcheck=bounds`.
- `kdtree2_r_count_around_point` leaked its query-vector allocation on
  every call (once per hub when `rmin > 0`).
- The 2PCF output header mislabelled the mu-edge columns as `r2 min/max`.
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
