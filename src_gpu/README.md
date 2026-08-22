# GRAMSCI GPU

GPU-accelerated build of GRAMSCI using OpenACC (NVIDIA HPC SDK / nvfortran).
Produces `bin/gramsci_gpu`, a drop-in replacement for the CPU `gramsci` binary
with identical command-line options and output formats.

## Requirements

- NVIDIA HPC SDK (nvfortran) — https://developer.nvidia.com/hpc-sdk
- An NVIDIA GPU (any recent compute capability; default build targets cc86)
- Linux + glibc (the build uses `malloc_trim` to control host memory)

## Building

```sh
cd src_gpu
make                # default GPU target: cc86 (RTX 30xx / A40 class)
make GPU=cc70       # Volta
make GPU=cc80       # Ampere A100
make GPU=cc90       # Hopper
```

The CPU modules from `../src` are recompiled here with nvfortran (gfortran
`.mod` files are incompatible), then the GPU modules are linked on top.
The binary is installed to `../bin/gramsci_gpu`.

## Usage

Identical to the CPU binary, e.g.:

```sh
gramsci_gpu -gal gals.dat -ran randoms.dat \
    -rmin 50.0 -rmax 150.0 -nbins 20 -nmu 1 \
    -out result.3pcf -3pcf
```

Input catalogs are 4-column ASCII (`x y z weight`, comoving Mpc/h).

### What runs where

| Mode      | Execution                                  |
|-----------|--------------------------------------------|
| `-2pcf`   | CPU (O(edges), not worth GPU offload)       |
| `-3pcf`   | GPU (isotropic); CPU fallback when `-nmu` > 1 (RSD) |
| `-equi`   | GPU (isotropic); CPU fallback when `-nmu` > 1       |
| `-4pcf`   | GPU                                         |
| `-4pcfp`  | GPU (parity decomposition)                  |

Graph construction (kd-tree pair finding) always runs on the CPU with OpenMP
and typically takes a small fraction of the total runtime.

## Memory handling

All device sizing is determined at runtime from free GPU memory:

- Per-hub scratch (the 4PCF connectivity matrix) is sized from the longest
  adjacency row in the actual graph.
- When the flattened graph (CSR) fits on the device, kernels run single-pass.
- When it does not, the query automatically switches to a **chunked mode**:
  the edge arrays are split into row-aligned windows and the kernels run as
  tiles over (hub-window × search-window) combinations, with only 2-3 windows
  resident on the device at a time.  Measured overhead vs. a hypothetical
  all-resident run is ~20% (PCIe re-transfers); arbitrarily large graphs are
  supported (tested to 9×10⁹ edges, 45 GB of CSR, on a 24 GB card).
- Edge indexing is 64-bit throughout; graphs may exceed 2³¹ edges.
- Host RAM is kept to roughly one live copy of the graph: jagged rows are
  freed (and returned to the OS) as the CSR is filled.

`GRAMSCI_GPU_WIN_EDGES=<n>` overrides the per-window edge budget — useful for
testing the chunked paths on small data, or constraining device usage on a
shared GPU.

## Output notes

- Counts are accumulated with floating-point atomics; results validate
  exactly against the CPU binary in practice, but bit-level reproducibility
  across runs is not formally guaranteed.
- 4PCF parity (`-4pcfp`): tetrahedra whose pixelized direction vectors give a
  (mathematically) zero signed volume have no defined chirality and contribute
  0 to the parity-odd channel.  The CPU code applies the same rule
  (`VOL_DEGEN_TOL` in `src/query_4pcf_module.F90`); without it the odd channel
  picks up compiler-dependent floating-point noise.
- The discarded fraction and the odd-channel attenuation both fall steeply
  with grid resolution: measured over the estimator's own quadruplet measure,
  6.8% of tetrahedra are degenerate at the default `4 x 16` grid, 1.4% at
  `8 x 32` and 0.12% at `23 x 91`.  `-ntheta N` / `-nphi M` set the grid
  (`N*M <= 32767`, 2 bytes per edge); `-exactparity` computes the signed
  volume from the galaxy positions instead, which discards nothing, has no
  attenuation, needs no per-edge index (24 B per point instead of 2 B per
  edge) and costs ~8% more query time.

## Validation

`tests/benchmark_3pcf.py` compares CPU and GPU outputs and timings:

```sh
cd tests
python benchmark_3pcf.py --mode 3pcf --rmax 50    # also: --mode 4pcf
```

All four kernels (3PCF, equilateral, 4PCF, 4PCF-parity) have been validated
against the CPU reference in both single-pass and chunked modes with exact
agreement.  To force the chunked paths on a small catalog:

```sh
GRAMSCI_GPU_WIN_EDGES=10000000 gramsci_gpu ... 
```

## Performance reference

DESI LRG-scale catalog (2.95M points, 1:1 randoms), RTX 3090 Ti vs. 64-thread
EPYC, query time only:

| Measurement                 | CPU×64  | GPU     | Speedup |
|-----------------------------|---------|---------|---------|
| 3PCF, rmax=80               | 9.4 s   | 3.6 s   | 2.6×    |
| 4PCF, rmax=65               | 147 s   | 16.4 s  | 9.0×    |
| 3PCF, rmax=150 (9.0×10⁹ edges, chunked) | —  | 1674 s  | —       |

## Implementation notes (for developers)

Design constraints discovered the hard way; see comments in the source for
detail:

- One gang (CUDA block) per hub, explicit `!$ACC LOOP VECTOR` on the inner
  neighbor loop.  A bare `PARALLEL LOOP GANG` lets nvfortran place a different
  hub on every vector lane — divergent memory access, ~6× slower than CPU.
- Large per-hub arrays cannot go in `PRIVATE` when `VECTOR_LENGTH > 1`
  (nvfortran allocates them per lane).  Per-gang scratch is carved from a
  global array indexed by an explicit gang-slot loop.
- Accumulators are slot-strided partial arrays updated with `ATOMIC UPDATE`
  and summed on the host.  nvfortran array reductions on multidimensional
  allocatables do not reliably propagate to the host, and gang+vector array
  reductions re-reduce per inner-loop instance.
- Scalar arguments of `!$ACC ROUTINE` device functions must be declared
  `VALUE`: by-reference scalars receive a one-shot device copy that is *not*
  refreshed on later kernel launches — host updates between launches are
  silently ignored (observed as out-of-bounds reads in chunked tiles).
