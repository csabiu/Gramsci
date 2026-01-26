# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

GRAMSCI (GRAph Made Statistics for Cosmological Information) is a Fortran-based tool for calculating N-point spatial correlation functions (2PCF, 3PCF, 4PCF) of 3D point sets. It uses a KD-tree data structure for efficient neighbor searches and supports OpenMP parallelization. The code is designed for cosmological analysis of galaxy distributions.

**Paper**: [arXiv:1901.00296](https://arxiv.org/abs/1901.00296) - The Astrophysical Journal Supplement Series, Volume 242, Issue 2 (2019)

## Build Commands

```bash
# Build all executables (gramsci and domain_decomposition)
make all

# Run tests
make test

# Clean build artifacts
make clean
```

The build requires `gfortran` with OpenMP support. Executables are placed in `bin/`.

## Running the Code

```bash
# 2-point correlation function
./bin/gramsci -gal <galaxy_file> -ran <random_file> -rmin <Rmin> -rmax <Rmax> -nbins <N> -out <output> -2pcf

# 3-point correlation function
./bin/gramsci -gal <galaxy_file> -ran <random_file> -rmin <Rmin> -rmax <Rmax> -nbins <N> -out <output> -3pcf

# 4-point correlation function
./bin/gramsci -gal <galaxy_file> -ran <random_file> -rmin <Rmin> -rmax <Rmax> -nbins <N> -out <output> -4pcf
```

Key parameters:
- `-rmin`, `-rmax`: Radial separation range
- `-nbins`: Number of radial bins
- `-nmu`: Number of angular (mu) bins for anisotropic analysis
- `-gal`, `-ran`: Galaxy and random catalog files (format: x y z [weight])

## Architecture

```
src/
├── gramsci.F90              # Main program - correlation function calculations
├── kdtree2.F90              # KD-tree implementation for spatial searches
├── node_module.F90          # Node data structure
├── sorting_module.F90       # Sorting utilities
├── extension.F90            # Helper extensions
└── domain_decomposition.f90 # Domain decomposition for large datasets
```

The code uses:
- Conditional compilation for MPI support (`-DMPI`)
- OpenMP for thread parallelization
- Double precision (`-DDOUBLE` flag)

## Testing

Tests are in `tests/run_correlation_tests.py` and verify 2PCF and 3PCF calculations against expected DD/RR and DDD/RRR pair counts. Run with `make test`.

## Example Data

The `example/` directory contains test galaxy (`test.gal`) and random (`test.ran`) catalogs, plus example scripts showing typical usage patterns.
