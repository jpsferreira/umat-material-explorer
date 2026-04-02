# UMAT Material Explorer

Explore constitutive material behavior under different loading scenarios and fit material parameters to experimental data. Works with any ABAQUS-compatible UMAT subroutine — no ABAQUS installation required.

Ships with a GHO (Generalized Humphrey-Ogden) viscoelastic UMAT as a working example. Replace it with your own UMAT to explore any material model.

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
Genetic algorithm (SGA) to fit material parameters to experimental stress-strain data from uniaxial tests.

## Prerequisites

- **gfortran** 7+ (or any Fortran compiler)
- **GNU Make**
- **gnuplot** (optional, for visualization)

## Quick Start

```bash
# Build all
make

# Run monotonic loading
make -C monotonic run

# Run cyclic loading
make -C cyclic run

# Run parameter fitting
make -C fitting run
```

## Using Your Own UMAT

The test drivers call the standard ABAQUS `subroutine umat(...)` interface. To use a different material model:

1. Replace `umat/umat.for` with your UMAT source (keep the same filename, or update the Makefiles)
2. Update `umat/param_umat.inc` with your NSDV (number of state variables)
3. Edit the PROPS assignments in the test driver (`monotonic.f90`, `cyclic.f90`) to match your material parameters
4. For the fitting module: update `fitting/umat.for` and the parameter bounds in `fitting/ga.inp`

## Included Example: GHO Viscoelastic Model

The shipped UMAT implements:
- **Isotropic matrix**: Neo-Hookean (C10, C01)
- **Anisotropic fibers**: HGO with dispersion (K1, K2, kappa)
- **Viscoelasticity**: Generalized Maxwell (up to 3 branches)
- **Volumetric**: Penalty formulation (KBULK)

## Project Structure

```
umat/                 Shared UMAT source and includes
  umat.for              Constitutive model (replace with your own)
  param_umat.inc        State variable dimensions
  aba_param.inc         ABAQUS compatibility include
  resetdfgr.for         Deformation gradient reset utility
  fibers.inp            Fiber orientation data (model-specific)

monotonic/            Monotonic loading driver
cyclic/               Cyclic loading driver (freq + amplitude sweeps)
fitting/              Genetic algorithm parameter fitting
  umat.for              Fitting-specific UMAT variant
  sga.f95               Simple Genetic Algorithm implementation
  soft_tissue.csv       Experimental data (uniaxial stress-strain)
  ga.inp                GA configuration (population, generations, bounds)
```
