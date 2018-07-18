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
* `aero` folder must be placed in a `PYTHONPATH` structure
* or append `$PWD`to your `PYTHONPATH`

    PYTHONPATH=$PWD:$PYTHONPATH python examples/test.py

### Requirements
* numpy
* examples are plotted using [matplotlib](http://matplotlib.org)

io
--

The _io_ part of _hades_ adds a set of writers and readers to manipulate different data format without difficulties.

#### Features
### Readers
* Charles X finite volume reader
* Gmsh reader
* Vtk reader
* Icem unstructured mesh reader
### Writers
* Charles X finite volume writer
* Vtk writer

#### Installation & useage
* A new bash script, hades.env, has been added at the root of hades repertory.
    * Adding `source path_to_hades/hades.env` to your `.bashrc` or sourcing it in a given terminal will give you access to all the hades io writers and readers in any python script.
    * A reader can then be loaded via `from hades.io.Readers.XXXX import *` and a writer can be loaded via `from hades.io.Writers.XXXX import *`. For readers, one can choose among `[readVtk, readGmsh, readRestartIC3, readIcemCfdUns]`. For writers, among `[writeVtk, writeRestartIC3]`.
    * Subsequently, it grants access to the `[ReaderVtk, ReaderGmsh, ReaderRestartIC3, ReaderIcemCfdUns]` and `[WriterVtk, WriterRestartIC3]` classes.
* For the icem and ic3 readers and/or writers, compilation of the C extensions is *mandatory*. The icem reader can only be compiled on a machine where icemcfd is already installed.
    * Go to `path_to_hades/hades/io`
    * Modify as necessary either `Makefile.linux` or `Makefile.macosx` with the paths to the different libraries
    * Compile all by doing `make -f Makefile.linux` or `make -f Makefile.macosx`

#### Requirements
* In terms of python modules:
    * numpy
    * collections
    * vtk
    * struct
* Otherwise:
    * compiler C
    * compiler C++
    * libraries and headers IcemCFD

