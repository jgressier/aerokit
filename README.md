# aerokit

Python packages for compressible flow computations

[![PyPi Version](https://img.shields.io/pypi/v/aerokit.svg?style=flat)](https://pypi.org/project/aerokit)
[![PyPI pyversions](https://img.shields.io/pypi/pyversions/aerokit.svg?style=flat)](https://pypi.org/pypi/aerokit/)
[![Doc](https://readthedocs.org/projects/aerokit/badge/?version=latest)](https://readthedocs.org/projects/aerokit/)
[![Slack](https://img.shields.io/static/v1?logo=slack&label=slack&message=contact&style=flat)](https://join.slack.com/t/isae-opendev/shared_invite/zt-obqywf6r-UUuHR4_hc5iTzyL5bFCwpw
)

[![GitHub stars](https://img.shields.io/github/stars/jgressier/aerokit.svg?style=flat&logo=github&label=Stars&logoColor=white)](https://github.com/jgressier/aerokit)
[![PyPi downloads](https://img.shields.io/pypi/dm/aerokit.svg?style=flat)](https://pypistats.org/packages/aerokit)
[![codecov](https://img.shields.io/codecov/c/github/jgressier/aerokit.svg?style=flat)](https://codecov.io/gh/jgressier/aerokit)
[![lgtm](https://img.shields.io/lgtm/grade/python/github/jgressier/aerokit.svg?style=flat)](https://lgtm.com/projects/g/jgressier/aerokit/)

### Features

* Mach dependent functions for isentropic total pressure, temperature and mass flow
* local Rankine-Hugoniot shock wave equations (or 2d planar shocks)
* conical shock waves
* Fanno    equations for momentum losses in a duct
* Rayleigh equations for heating/cooling in a duct
* unsteady compressible 1D flows
* applications: nozzle flows and (Riemann) shock tube
* misc: degree based trigo functions, Newton iterative solve, ODE integration

### Installation & usage

    pip install aerokit

Documentation is available on [readthedocs](https://aerokit.readthedocs.io/en/latest/). More examples are available in the repository.

### Requirements

* `numpy`
* `scipy`
* examples are plotted using [`matplotlib`](http://matplotlib.org)
