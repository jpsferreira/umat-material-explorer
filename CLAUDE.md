# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Fortran codebase for exploring ABAQUS-compatible UMAT (User MATerial) subroutine behavior under different loading scenarios without requiring an ABAQUS installation. Ships with a GHO (Generalized Humphrey-Ogden) viscoelastic model as a working example.

## Build and Run Commands

```bash
# Build everything
make

# Build and run individual modules
make -C monotonic run    # monotonic loading (uniaxial, biaxial, pure shear, simple shear)
make -C cyclic run       # cyclic loading (frequency and amplitude sweeps)
make -C fitting run      # genetic algorithm parameter fitting

# Clean all build artifacts
make clean
```

Compiler: `gfortran` with `-O2` optimization. Visualization: `gnuplot` (optional, called automatically after runs).

## Architecture

### UMAT Interface

The core pattern: each module is a standalone Fortran program that sets up deformation gradients (DFGRD1), material properties (PROPS array), and calls `SUBROUTINE UMAT(...)` using the standard ABAQUS interface. There is no ABAQUS dependency — the drivers replicate what ABAQUS would do.

### Module Structure

- **`umat/`** — Shared UMAT source. `umat.for` is the constitutive model, `resetdfgr.for` resets deformation gradients, `param_umat.inc` defines NSDV and constants, `aba_param.inc` provides ABAQUS compatibility types.
- **`monotonic/`** — Single driver (`monotonic.f90`) that applies monotonic stretch/shear and writes `.out` files to `stress_curves/`. Links against `../umat/umat.for` and `../umat/resetdfgr.for`.
- **`cyclic/`** — Single driver (`cyclic.f90`) for frequency and amplitude sweeps. Same linking pattern as monotonic.
- **`fitting/`** — Self-contained module with its own copy of UMAT source. Uses a Simple Genetic Algorithm (`sga.f95`) driven by `main.f90`. Configuration in `ga.inp`, parameter bounds in `params.inc`. `uexternaldb.for` and `uexternald.for` handle I/O for fiber data and experimental data (`soft_tissue.csv`).

### Key Conventions

- Material properties are assigned positionally in the PROPS array: PROPS(1)=KBULK, PROPS(2)=C10, PROPS(3)=C01, PROPS(4)=K1, PROPS(5)=K2, PROPS(6)=kdisp, PROPS(7..13)=viscoelastic Maxwell parameters.
- State variable count (NSDV) is set in `param_umat.inc` — must match between UMAT and drivers.
- The fitting module has its own `param_umat.inc`, `aba_param.inc`, and `umat.for` (may differ from the shared version in `umat/`).
- Output goes to `stress_curves/` directories; gnuplot scripts in `plots/` subdirectories visualize results.
- Fiber orientation data is read from `fibers.inp` (symlinked from `umat/` at runtime via `make run`).

### Replacing the UMAT

To use a different material model: replace `umat/umat.for`, update NSDV in `umat/param_umat.inc`, and adjust PROPS assignments in the driver programs. For fitting, update the copies under `fitting/` independently.
