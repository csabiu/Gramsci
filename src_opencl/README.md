# GRAMSCI OpenCL (portable GPU build)

A portable, OpenCL-based GPU build of GRAMSCI driven from Fortran.  Unlike the
`src_gpu/` build (OpenACC / NVIDIA HPC SDK, NVIDIA-only), this backend runs on
**any OpenCL 1.2 device** — in particular **Apple Silicon and Intel Macs**,
where it offloads the N-point counting to the integrated GPU.

It produces `bin/gramsci_cl`, a drop-in replacement for the CPU `gramsci`
binary with identical command-line options and output formats.

## Requirements

- `gfortran` (Homebrew GCC on macOS: `brew install gcc`)
- An OpenCL runtime:
  - **macOS**: built in (`/System/Library/Frameworks/OpenCL.framework`) — nothing to install
  - **Linux**: any ICD (NVIDIA / AMD / Intel / POCL) + headers

No third-party Fortran-OpenCL library is needed: the OpenCL C API is bound
directly with `ISO_C_BINDING` in `cl_module.F90`.

## Building

```sh
cd src_opencl
make            # builds ../bin/gramsci_cl
make clean
```

The CPU modules from `../src` are recompiled here (same gfortran, so the `.mod`
files are compatible) and the OpenCL backend is linked on top.  On macOS the
Makefile auto-detects the SDK framework path; on Linux it links `-lOpenCL`.

## Usage

Identical to the CPU binary, e.g.:

```sh
gramsci_cl -gal gals.dat -ran randoms.dat \
    -rmin 1.0 -rmax 40.0 -nbins 6 -out result.3pcf -3pcf
```

Input catalogs are 4-column ASCII (`x y z weight`, comoving Mpc/h).

### What runs where

| Mode      | Execution                                            |
|-----------|------------------------------------------------------|
| `-2pcf`   | CPU (O(edges); not worth GPU offload)                |
| `-3pcf`   | GPU (isotropic, incl. `-njk`); CPU fallback when `-nmu` > 1 (RSD) |
| `-equi`   | GPU (isotropic); CPU fallback when `-nmu` > 1 or `-njk` > 0 |
| `-4pcf`   | GPU (incl. `-njk`)                                   |
| `-4pcfp`  | GPU (parity, incl. `-njk`); CPU fallback for `-exactparity` |

Graph construction (kd-tree pair finding) always runs on the CPU with OpenMP.

**Jackknife (`-njk`)** runs on the GPU for the isotropic 3PCF and both 4PCF
modes.  Apple's OpenCL has no float atomics, so the per-region touching sums
are accumulated with **CAS-emulated float atomic adds** into shared device
buffers laid out to match the host `N2jk`/`N4jk` arrays.  Those buffers get
no Kahan compensation (an atomic pair update is impossible), so with `-njk`
the query always runs through the **tiled/bucketed launcher**: every window's
fp32 partials are committed into double host accumulators and re-zeroed on
the device, which bounds the fp32 error per slot to a single window's worth
of tuples (same mechanism that protects against the watchdog).  Measured
CPU-vs-OpenCL agreement of all jackknife sidecars: ~2e-5 relative-to-peak.
`-exactparity` (any `-njk`) runs the CPU parity kernel: the OpenCL parity
chirality is looked up from the host-built pixel table and has no
exact-positions path.

## Precision (important)

**The OpenCL kernels run in single precision (fp32).**  Apple's OpenCL on
Apple Silicon reports `CL_DEVICE_DOUBLE_FP_CONFIG == 0` (no `double` in
kernels), so all device arithmetic is `float`.  To keep accuracy:

- Each work-item accumulates into its **own** partial-histogram column with
  **Kahan compensated summation** (a paired compensation buffer per partial);
  the columns are summed back into the final counts **in double on the host**
  as `sum − comp`.  The error stays O(ε) *independent of the tuple count* —
  a plain fp32 `+=` would silently stop registering increments once a column
  passed ~1.7×10⁷ tuples, one-sidedly undercounting the monotone RRR/RRRR
  channels on production-size runs.  This also avoids device atomics
  entirely (Apple OpenCL has no fp32 atomics).
- The 4PCF **parity** chirality sign is precomputed on the host in double (same
  `VOL_DEGEN_TOL` as the CPU code) into a per-pixel-triple sign table and looked
  up in the kernel — so the parity channel is bit-exact in *sign* and only the
  weight magnitude is fp32.

Measured agreement vs. the double-precision CPU reference (test catalog,
707k points):

| Query        | max relative error vs CPU |
|--------------|---------------------------|
| 3PCF, equi   | ~1e-8 (counts), RRR exact |
| 4PCF         | ~1e-7 (counts)            |
| 4PCF parity  | ~2e-7 (even & odd)        |

For full double precision (e.g. publication runs on very small signals), use
the CPU build (`src/`) or the NVIDIA OpenACC build (`src_gpu/`).

## Performance & the GPU watchdog

A single NDRange that runs too long trips the OS **GPU watchdog** (the
integrated GPU also drives the display); Apple's OpenCL then returns *partial*
results with no error.  Each query therefore runs as **one launch first**, with
per-work-item completion flags; if the watchdog truncated it, the backend
**re-zeros and re-runs the query as several shorter interleaved-bucket
launches** ("single GPU launch hit the watchdog; retrying tiled" is printed).
The tiled path re-checks the flags after *every* window (the first window is
sized open-loop, so it can trip the watchdog too), periodically commits
verified partials to double host accumulators, and on truncation discards the
uncommitted counts, rewinds to the last commit, and shrinks the window.
Results are correct either way; only the timing differs.
(`GRAMSCI_CL_FORCE_TILED=1` forces the tiled path, for testing.)

On Apple's GPU a single *sustained* launch is much faster than many short ones
(the GPU clock ramps up only under sustained load), so:

- **Small/medium queries** (whatever fits under the watchdog — on an 8-GPU-core
  M1, roughly ≤ ~1 s of device work) run in one launch and are competitive with
  or faster than the multi-threaded CPU.
- **Large queries** fall back to the tiled path, which is correct but can be
  *several times slower than the CPU* on a small integrated GPU.

This is a property of the hardware, not a bug: the M1's 8-core integrated GPU is
not dramatically faster than an 8-core CPU for this memory-bound graph traversal,
and the watchdog penalises the largest queries.  The OpenCL backend pays off
most on **discrete GPUs** (Linux AMD/NVIDIA, or higher-core Apple GPUs), where
the watchdog is generous and the core count is far higher.

`GRAMSCI_CL_TARGET_SEC` tunes the per-launch target of the tiled fallback (lower
it if you still see a watchdog message at very large `-rmax`).

## Validation

```sh
cd tests
bash ../src_opencl/validate.sh          # CPU vs OpenCL, all five modes
bash ../src_opencl/validate.sh 1e-3     # custom relative tolerance
```

## Scope / limitations (v1)

- **Single-pass only.** The whole CSR graph, the `nbins^6` 4PCF config table,
  and the accumulators must fit in device memory.  On Apple Silicon's unified
  memory this covers any laptop-scale catalog; the largest single buffer is
  `CL_DEVICE_MAX_MEM_ALLOC_SIZE` (≈2.2 GB on an 8 GB M1, so ≈5×10⁸ edges in
  `csr_id`).  The out-of-core chunking of the OpenACC build is **not** ported.
  A clear error is printed if the `nbins^6` table exceeds the max buffer size.
- **RSD (`-nmu > 1`) is not offloaded** — those queries run on the CPU, exactly
  as in the OpenACC build.
- **GPU watchdog → fast launch + tiled fallback.** See "Performance & the GPU
  watchdog" above: each query runs as one launch, with completion-flag detection
  of watchdog truncation and an automatic re-run as shorter interleaved-bucket
  launches.  Correct always; large queries on a small integrated GPU are slower.
- **4PCF uses binary search, not the `lmat` connectivity-matrix optimisation**
  of the OpenACC build.  `lmat` needs a per-work-item read-write *global*
  scratch buffer; the bsearch kernel uses only read-only inputs + a private
  accumulator column (structurally identical to the reliable 3PCF kernel), so
  it has no scratch to mismanage and is deterministic.  It does ~25× more
  searches per hub, so on a laptop's integrated GPU the 4PCF query can be
  *slower* than the multi-threaded CPU; the win is portability and offloading
  the host, and it scales to discrete OpenCL GPUs.

## Files

| File                         | Role                                              |
|------------------------------|---------------------------------------------------|
| `cl_module.F90`              | minimal `ISO_C_BINDING` interface to OpenCL 1.2    |
| `cl_env_module.F90`          | device/context/queue/program + buffer/arg helpers  |
| `kernels.cl`                 | OpenCL C compute kernels (3PCF, equi, 4PCF, 4PCFp)  |
| `cl_kernels_module.F90`      | auto-generated: `kernels.cl` embedded as a string  |
| `embed_kernels.sh`           | generator for the above (run by the Makefile)      |
| `csr_cl_module.F90`          | CSR flattening (portable; no glibc `malloc_trim`)  |
| `query_3pcf_cl_module.F90`   | host orchestration for 3PCF / equilateral          |
| `query_4pcf_cl_module.F90`   | host orchestration for 4PCF / 4PCF-parity          |
| `gramsci_cl.F90`             | main driver                                        |
| `validate.sh`                | CPU-vs-OpenCL regression check                      |
