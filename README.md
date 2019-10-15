[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

# Finite Element Mesh Poisson Solver (FEMPS)

(C) Copyright 2018-2019 UCAR and 2011-2018 John Thuburn, University of Exeter, UK.

This software is licensed under the terms of the Apache Licence Version 2.0 which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

Finite element mesh Poisson solver for cubed sphere and hexagonal icosahedral grids. Originally developed by John Thuburn, University of Exeter. Code here is refactored for use with the JEDI data assimilation system.

**To build:**

Dependencies:
- Fortran compiler
- Git
- Git-lfs (https://git-lfs.github.com/)
- netCDF
- ecbuild (https://github.com/ecmwf/ecbuild)

Build steps:
1. `git clone https://github.com/.../femps`
2. `cd femps`
3. `mkdir build`
4. `cd build`
5. `ecbuild ../` [Optionally specify `--build=debug` for debug flags]
6. `make`
7. `ctest -V`
