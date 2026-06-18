# Changelog

Notable changes to GRAMSCI. This project accompanies Sabiu, Hoyle, Kim & Li,
*ApJS* **242**, 29 (2019), [arXiv:1901.00296](https://arxiv.org/abs/1901.00296).

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
