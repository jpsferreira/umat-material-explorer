# UMAT Material Explorer

Explore the mechanical behavior of soft biological tissues under different loading scenarios and fit material parameters to experimental data. Uses a standalone GHO (Generalized Humphrey-Ogden) viscoelastic UMAT — no ABAQUS required.

## Loading Scenarios

### Monotonic Loading (`monotonic/`)
Continuous loading to study stress-strain response under:
- Uniaxial tension
- Equibiaxial tension
- Pure shear
- Simple shear

### Cyclic Loading (`cyclic/`)
Oscillatory loading to characterize viscoelastic behavior:
- **Frequency sweeps**: storage modulus, loss modulus, tan(delta) vs frequency
- **Amplitude sweeps**: same quantities vs strain amplitude

### Parameter Fitting (`fitting/`)
Genetic algorithm (SGA) to fit material parameters (C10, K1, K2, kappa) to experimental PK2 stress-strain data from uniaxial tests.

## Prerequisites

- **gfortran** 7+ (or any Fortran compiler)
- **GNU Make**
- **gnuplot** (optional, for visualization)

## Quick Start

```bash
# Build all
make

# Run monotonic tests
make -C monotonic run

# Run cyclic tests
make -C cyclic run

# Run parameter fitting
make -C fitting run
```

## Material Model

The GHO model combines:
- **Isotropic matrix**: Neo-Hookean (C10, C01)
- **Anisotropic fibers**: HGO with dispersion (K1, K2, kappa)
- **Viscoelasticity**: Generalized Maxwell (up to 3 branches)
- **Volumetric**: Penalty formulation (KBULK)

The UMAT source is in `umat/umat_gho.for` (shared by monotonic and cyclic tests). The fitting module uses its own UMAT variant in `fitting/umat_gho.for`.

## Project Structure

```
umat/                 Shared UMAT source and includes
  umat_gho.for          GHO viscoelastic constitutive model
  param_umat.inc        Material parameters and dimensions
  aba_param.inc         ABAQUS compatibility include
  resetdfgr.for         Deformation gradient reset utility
  fibers.inp            Fiber orientation data

monotonic/            Monotonic loading test driver
cyclic/               Cyclic loading test driver (freq + amplitude sweeps)
fitting/              Genetic algorithm parameter fitting
  sga.f95               Simple Genetic Algorithm implementation
  soft_tissue.csv       Experimental data (uniaxial PK2 stress-strain)
  ga.inp                GA configuration (population, generations, bounds)
```
