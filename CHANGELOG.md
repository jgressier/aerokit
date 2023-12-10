# Changelog & Release Notes

## Upgrading

To upgrade to the latest version of `aerokit` use `pip`:

```bash
pip install aerokit --upgrade
```

You can determine your currently installed version using this command:

```bash
pip show aerokit
```

## Versions

### [1.2.1](https://pypi.org/project/aerokit/) (2023-12-10)

#### new

- `common.numspectral` extrapolation and fitting
- `stability.Euler` with Rankine-Hugoniot conditions

#### fixed

- compilation of documentation on readthedocs
 
### [1.2.0](https://pypi.org/project/aerokit/) (2023-11-22)

#### new

- `common.numspectral` subpackage with Chebyshev derivation operator
- `stability.OrrSommerfeld` with new Orr-Sommerfeld model
- `stability.Euler` with new Linearized Euler model
- eigenvalue problem, solver, sorting functions, convergence with Rayleigh method

### [1.1.3](https://pypi.org/project/aerokit/) (2021-06-02)

#### fixed

- `instance.riemann`: fix wave speeds for right wave if it is an expansion

### [1.1.2](https://pypi.org/project/aerokit/) (2021-04-15)

#### fixed

- remove debug print in `aero.ShockWave` conical functions
- add reverse function `conical_Mach_walldeflection_sigma` (upstream Mach number from conical wall deflection and shock angle)

### [1.1.1](https://pypi.org/project/aerokit/) (2021-03-29)

#### changed

- save and restore default in `common.defaultgas`
- change Pi to Pt in `ShockWave` functions (backward compatibility handled)

#### fixed

- improve `aero.MassFlow` initialization for Sigma to Mach computations
- `instance.nozzle` is now using prescribed gamma
- some recommended fixes by LGTM code analysis
- improve test coverage

### [1.1.0](https://pypi.org/project/aerokit/) (2021-02-05)

#### changed

- changed total/stagnation notation Ti,Pi to Tt,Pt (backward compatibility handled)
- moved `nozzle` and `riemann` to `instance.*` submodule 

### [1.0.0](https://pypi.org/project/aerokit/) (2021-01-20)

- `aero.Isentropic` submodule: classical compressible isentropic flow functions
- `aero.ShockWave` submodule: local (2D) and conical shock functions (shock angle and deviation)
- `aero.Supersonic` submodule: supersonic 2D invariants, Busemann/Prandtl meyer functions
- `aero.Massflow` submodule: 1D massflow functions
- `aero.Fanno` submodule: theoretical solution of Fanno 1D flow (momentum source)
- `aero.Rayleigh` submodule: theoretical solution of Rayleigh 1D flow (energy source)
