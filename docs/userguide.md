# User Guide

The `aerokit` package mainly consists in `aerokit.aero` module for various gas function computations. Default gas definition is handled in `aerokit.common.defaultgas` module.

Most implemented functions use normalized properties.

## `aerokit.aero` module

- `Isentropic`: isentropic functions for compressible ideal flows (total quantities)
- `MassFlow`: mass flow normalization, choking, 1D nozzle
- `Supersonic`: supersonic specific functions (invariants, Prandtl-Meyer/Busemann function)
- `ShockWave`: 2D or conical shockwaves (Rankine-Hugoniot equations)
- `Propulsion`: flow functions for thrust
- `Fanno`: Fanno 1D flow (momentum source)
- `Rayleigh`: Rayleigh 1D flow (energy source)
- `model1D`: class for 1D compressible state
- `unsteady1D`: extended class for 1D unsteady computations

## `aerokit.instance` module

- `riemann`: solution for generalized Riemann problems (shock tube)
- `nozzle`: solution of 1D nozzle problem
