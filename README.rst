aero
----

Python packages for compressible flow computations

#### Features
* Mach dependent functions for isentropic total pressure, temperature and mass flow
* local Rankine-Hugoniot shock wave equations
* conical shock waves
* real gaz thermo functions
* kerozen  thermo functions
* Fanno    equations for momentum losses in a duct
* Rayleigh equations for heating/cooling in a duct
* misc: degree based trigo functions, Newton iterative solve, ODE integration

#### Installation & usage
* `aerokit` folder must be placed in a `PYTHONPATH` structure
* or append `$PWD`to your `PYTHONPATH`

    PYTHONPATH=$PWD:$PYTHONPATH python examples/test.py

### Requirements
* numpy
* examples are plotted using [matplotlib](http://matplotlib.org)

