# Changelog

Notable changes to GRAMSCI. This project accompanies Sabiu, Hoyle, Kim & Li,
*ApJS* **242**, 29 (2019), [arXiv:1901.00296](https://arxiv.org/abs/1901.00296).

## [Unreleased]

### Added
- **Delete-one jackknife errors for every statistic** — `-2pcf`, `-3pcf`,
  `-equi`, `-4pcf` and the parity-odd `-4pcfp` (previously 3PCF-only).
  `-njk N` alone now partitions the sky **internally** into N equal-count
  angular regions (sin(dec) bands × φ slices at quantiles of the random
  catalogue, applied identically to both catalogues); the partition is
  deliberately angular-only so radial and angular systematics are never
  mixed in the error estimate. The labels used are written to
  `<out>.jkgal`/`<out>.jkran` in the exact format `-jkgal`/`-jkran` read
  back. All N realisations still come from the single normal query pass
  (`N_m = N_total − N_touching(m)`, per-tuple distinct-region dedup); each
  mode writes `<out>[.<mode>].jk` (realisations) and a new
  `<out>[.<mode>].jkerr` (jackknife mean and σ per bin). The 4PCF
  disconnected term is rebuilt per realisation from the jackknifed internal
  2PCF ξ₀ (ξ₂/ξ₄ stay full-sample); the odd channel has no disconnected
  term, so `zeta_conn_odd = zeta_odd` per realisation. The 4PCF jackknife
  accumulators use `!$OMP ATOMIC` on shared arrays instead of an OpenMP
  reduction (a per-thread copy of an (n_configs, nch, njk) array would
  multiply memory by the thread count), with a 64 GB guard. Python API:
  `njk=` on `compute_2pcf/3pcf/4pcf`, exposing `*_jk_mean`, `*_jk_sigma`
  and the raw realisation table `.jk_real`. Requirements: `-ran` and
  `N ≥ 2`; incompatible with analytic randoms; internal partition refused
  in `-box` mode (no observer direction — external labels still accepted).
  Verified by new brute-force regression tests that recompute every mode's
  delete-one realisations independently in Python (including the signed
  `bintable6` fill and `-exactparity` chirality signs) and by a
  combined-vs-single-mode identity test.

### Changed
- **Graph-build search scratch is no longer O(N) per thread.** create_graph
  allocated a kdtree2 results buffer sized num_data+num_rand for every
  OpenMP thread (~16 B x N x threads of address space; at N = 1e8 on 64
  threads that is ~100 GB of allocation, which fails outright on systems
  with strict memory accounting). The per-thread scratch is now sized from
  sampled neighbor counts and grown on demand; the search's true nfound is
  used for an exact re-allocation and retry on the rare overflow. Outputs
  verified identical on a 2M-point catalogue.

### Removed
- **Dead code left over from the 2.3.0 MPI/domain-decomposition removal**:
  the `buffer` halo-flag array (allocated, zeroed and checked at the top of
  every hot query loop on every backend -- CPU, OpenACC and OpenCL kernels
  alike -- but never set to anything but 0) and the never-read
  `loadran`/`saveran`/`ranfile` config fields. The OpenCL kernels lose
  their `buffer` argument, so every kernel argument index after it shifts
  down by one (wrappers renumbered; embedded kernel module regenerated).
  Verified: full CPU suite, CPU-vs-OpenCL agreement at file precision, and
  `src_opencl/validate.sh` on an Apple M1 -- 2PCF/3PCF/equi/4PCF pass
  identically to the pre-change baseline (the 4PCFp mismatch it reports is
  a pre-existing OpenCL bug, present at baseline, fixed separately).

### Fixed
- **The graph-RAM estimate sampled only data points**: every point (data
  and randoms alike) is a hub, and with clustered data the data hubs see
  systematically more neighbors than the uniform randoms that dominate
  the catalogue -- measured 3.7x overestimation on a clustered mock (491
  vs a true all-hub mean of 132 neighbors). The estimator now samples the
  whole catalogue, lives in one shared routine instead of two per-driver
  copies, and the OpenCL driver -- which sampled the neighbor counts and
  then never printed anything -- now reports the estimate too.
- **Segfault at multi-million-point catalogues** (silent, no backtrace,
  inside the KD-tree build): `-Ofast` implies `-fstack-arrays`, and the
  vector-subscript gather `sum(the_data(c, ind(l:u)))` in
  build_tree_for_range materialized an O(N) temporary on the default 8 MB
  stack at the tree root. The average is now an explicit loop (no
  temporary) and the gfortran builds additionally compile with
  `-fno-stack-arrays` (measured perf-neutral) so no future O(N) temporary
  can crash. A 4M-point catalogue that segfaulted now runs; regression
  test runs the KD-tree build under a 1 MB stack limit (verified to kill
  the pre-fix binary).
- kdtree2's fixed-ball searches called `kdtree2_sort_results(nfound, ...)`
  after an overflow, heapsorting `nfound > nalloc` elements past the end
  of the results array; the sort is now capped at `nalloc` so callers can
  safely detect overflow (`nfound > nalloc`) and retry with more storage
  (the grow-on-demand scratch above relies on this).
- **OpenCL 4PCF parity was silently corrupted on the default direction
  grid** -- the third instance of the same int8 pixel overflow: the OpenCL
  CSR packed the direction-pixel index as `int8` and the kernel read it as
  signed `char`, so on the 8×32 = 256-pixel default grid (v2.5.0) indices
  above 127 went negative and indexed the host-built chirality sign table
  out of bounds. `src_opencl/validate.sh` reported 4PCFp relative errors
  of ~1e2. csr_phi is now int16, the kernel takes `short`, and a
  `cl_buf_in_i16` upload helper was added. validate.sh now passes all
  five modes (4PCFp at 1.4e-7 vs the double-precision CPU, matching the
  even channel), plain and GRAMSCI_CL_FORCE_TILED=1. The OpenACC path was
  already 16-bit; `-exactparity` was unaffected.
- **Periodic-box parity was silently corrupted on the default direction
  grid**: the periodic graph pass stored the direction-pixel index through
  an `int8` conversion (`graph_module.F90`), but pixel indices reach
  `ntheta*nphi` = 256 on the default 8×32 grid — anything above 127
  wrapped to a zero/negative index, so `dir_x/y/z` were read out of
  bounds and every `-box` + `-4pcfp` chirality sign was
  orientation-dependent garbage (a parity-odd fraction of ~0.01 instead
  of ~1 on an all-left-handed test catalogue). The non-periodic pass was
  already 16-bit (v2.5.0); the periodic branch was missed, and survived
  because the historical 4×16 = 64-pixel grid fit in `int8` and no test
  covered periodic parity. `-exactparity` was unaffected (no pixels).
  New regression test: randomly rotated interior chiral tetrahedra must
  give identical NNNN/NNNN_odd in periodic and non-periodic runs, and a
  near-unity parity-odd fraction.
- **Empty bins no longer write Inf/NaN**: the 2PCF, 3PCF (isotropic and
  RSD) and equilateral writers divided `N/RR` without a zero guard, so a
  bin with no random pairs/triplets (sparse catalogues, wide `rmax`,
  narrow μ bins) produced `0/0 = NaN` rows that silently poison anything
  averaging over the output. Empty bins now carry a zero estimate — the
  convention the 4PCF and jackknife writers already used. Regression test
  added (all estimator columns finite; zero-count rows are zero).
- The OpenCL driver never called `read_jk_regions()`, so `-jkgal`/`-jkran`
  were silently ignored (and RSD 3PCF with `-njk` would crash) in the
  `gramsci_cl` build. It now reads/assigns regions like the other drivers,
  and all `-njk` queries route to the CPU kernels (the OpenCL kernels have
  no jackknife accumulation). The OpenACC build keeps its GPU jackknife
  for the isotropic 3PCF and routes equilateral/4PCF jackknife to the CPU.
- Combined query modes now also re-zero the jackknife accumulators
  (`N2jk`/`N3jk`) between modes; previously a combined `-3pcf` run after
  another jackknife mode could inherit stale touching sums (latent — only
  the 3PCF used them before this release).

## [2.5.0] — 2026-08-25

### Added
- **Exact signed-volume parity** (`-exactparity`): the odd-channel sign is
  computed from the exact galaxy positions rather than quantized spoke
  directions, eliminating direction-pixel attenuation entirely (recovered
  asymmetry consistent with unity for every tetrahedron family, including
  near-degenerate flattened and sheared shapes) at ~1.1–1.2× query cost.
- **Runtime direction-pixel grid** (`-ntheta`/`-nphi`) for the pixelized
  parity mode; pixel indices are now 16-bit (products up to 32767 pixels).
- **Delete-one jackknife resampling for the 3PCF** (CPU and GPU).
- GitHub Actions CI.

### Changed
- The default parity direction grid is now **8 × 32** (was 4 × 16),
  removing the recovered-asymmetry attenuation of near-degenerate
  configurations at <10% extra query cost.

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
