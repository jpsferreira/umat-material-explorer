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

Compiler: `gfortran` with `-O2 -Wall`. Visualization: `gnuplot` (optional, called automatically after runs).

## Architecture

### UMAT Interface

The core pattern: each module is a standalone Fortran program that sets up deformation gradients (DFGRD1), material properties (PROPS array), and calls `SUBROUTINE UMAT(...)` using the standard ABAQUS interface. There is no ABAQUS dependency — the drivers replicate what ABAQUS would do.

### Module Structure

- **`umat/`** — Shared UMAT source. `umat.for` is the constitutive model, `resetdfgr.for` resets deformation gradients, `param_umat.inc` defines NSDV and constants, `aba_param.inc` provides ABAQUS compatibility types.
- **`monotonic/`** — Single driver (`monotonic.f90`) that applies monotonic stretch/shear and writes `.out` files to `stress_curves/`. Links against `../umat/umat.for` and `../umat/resetdfgr.for`.
- **`cyclic/`** — Single driver (`cyclic.f90`) for frequency and amplitude sweeps. Same linking pattern as monotonic.
- **`fitting/`** — Self-contained module with its own copy of UMAT source. Uses a Simple Genetic Algorithm (`sga.f95`) driven by `main.f90`. Configuration in `ga.inp`, parameter bounds in `params.inc`. `uexternaldb.for` and `uexternald.for` handle I/O for fiber data and experimental data (`soft_tissue.csv`).

### Key Conventions

- Output goes to `stress_curves/` directories; gnuplot scripts in `plots/` subdirectories visualize results.
- Fiber orientation data is read from `fibers.inp` (symlinked from `umat/` at runtime via `make run`).

### PROPS Array Layout (shared UMAT, NPROPS=13)

| Index | Parameter | Description |
|-------|-----------|-------------|
| 1 | KBULK | Bulk modulus / penalty parameter for near-incompressibility |
| 2 | C10 | Neo-Hookean coefficient (isotropic matrix) |
| 3 | C01 | Neo-Hookean coefficient (isotropic matrix) |
| 4 | K1 | HGO fiber stiffness |
| 5 | K2 | HGO fiber nonlinearity exponent |
| 6 | kdisp | Fiber dispersion parameter (kappa) |
| 7 | V | Number of active Maxwell branches (1–3) |
| 8 | tau1 | Relaxation time, Maxwell branch 1 |
| 9 | theta1 | Viscous weight, Maxwell branch 1 |
| 10 | tau2 | Relaxation time, Maxwell branch 2 |
| 11 | theta2 | Viscous weight, Maxwell branch 2 |
| 12 | tau3 | Relaxation time, Maxwell branch 3 |
| 13 | theta3 | Viscous weight, Maxwell branch 3 |

The fitting module uses only NPROPS=6 (indices 1–6) since it fits the elastic response without viscoelasticity.

### State Variables (STATEV)

The shared UMAT uses NSDV=27. State variables store the viscous overstress tensor (HV, 3x3) for each Maxwell branch:
- STATEV(1–9): HV for Maxwell branch 1 (row-major 3x3)
- STATEV(10–18): HV for Maxwell branch 2
- STATEV(19–27): HV for Maxwell branch 3

See `HVREAD`/`HVWRITE` subroutines in `umat/umat.for` for the packing layout.

### Fitting Module Divergence

The `fitting/` module intentionally maintains its own copy of the UMAT with key differences from `umat/umat.for`:
- **No viscoelasticity** — `INITIALIZE`, `SDVREAD`, `SDVWRITE` omit the viscous parameter `VV`; no `VISCO` call.
- **Fewer state variables** — NSDV=9 (vs 27), NSTATEV=3 in the driver.
- **Fewer properties** — NPROPS=6 (elastic only, no Maxwell parameters).
- STATEV(3) stores the first Piola-Kirchhoff stress PK1(1,1) for comparing against experimental data.

When updating the shared UMAT, the fitting copy must be updated independently. The two are not expected to stay in sync.

### Replacing the UMAT

To use a different material model: replace `umat/umat.for`, update NSDV in `umat/param_umat.inc`, and adjust PROPS assignments in the driver programs. For fitting, update the copies under `fitting/` independently.
