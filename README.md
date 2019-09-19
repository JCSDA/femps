# Finite Element Mesh Poisson Solver (FEMPS)

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
