# AeroKit

Python packages for compressible flow computations

[![PyPi Version](https://img.shields.io/pypi/v/aerokit.svg?style=flat)](https://pypi.org/project/aerokit)
[![PyPI pyversions](https://img.shields.io/pypi/pyversions/aerokit.svg?style=flat)](https://pypi.org/pypi/aerokit/)
[![Doc](https://readthedocs.org/projects/aerokit/badge/?version=latest)](https://aerokit.readthedocs.io/)
[![Slack](https://img.shields.io/static/v1?logo=slack&label=slack&message=contact&style=flat)](https://join.slack.com/t/isae-opendev/shared_invite/zt-obqywf6r-UUuHR4_hc5iTzyL5bFCwpw
)

[![GitHub stars](https://img.shields.io/github/stars/jgressier/aerokit.svg?style=flat&logo=github&label=Stars&logoColor=white)](https://github.com/jgressier/aerokit)
[![PyPi downloads](https://img.shields.io/pypi/dm/aerokit.svg?style=flat)](https://pypistats.org/packages/aerokit)
[![codecov](https://img.shields.io/codecov/c/github/jgressier/aerokit.svg?style=flat)](https://codecov.io/gh/jgressier/aerokit)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/cf3b4b76534f42a1a1fe7c913d2ee782)](https://app.codacy.com/gh/jgressier/aerokit/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_grade)

### Features

* Isentropic relations, compressible mass-flow functions, and converging-diverging nozzle flows
* Normal, oblique, and conical shock relations; shock polars and two-shock interactions
* Prandtl-Meyer expansions, Fanno flows, and Rayleigh flows
* One-dimensional unsteady compressible flows, including exact Riemann (shock-tube) solutions and inflow conditions
* Ideal two-dimensional blade rows and preliminary compressor maps
* Ideal turbojet and turbofan cycle-performance models
* Linear-stability tools based on Chebyshev collocation:
  * incompressible Orr-Sommerfeld analysis of parallel shear flows
  * linearized one-dimensional compressible Euler equations
  * local inviscid axisymmetric stability of compressible jets
* Real-gas and kerosene thermodynamic functions
* Degree-based trigonometric functions, nonlinear solvers, ODE integration, and spectral-discretization utilities

### Installation & usage

    pip install aerokit

Documentation is available on [readthedocs](https://aerokit.readthedocs.io/en/latest/). More examples are available in the repository.

### Requirements

* `numpy`, `scipy`
* examples are plotted using [`matplotlib`](http://matplotlib.org)
