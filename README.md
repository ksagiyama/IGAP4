# IGAP4
An Isogeometric Analysis Program for partial differential equations with 4th-order spatial derivatives.

## Overview
IGAP4 is written in C and it requires PETSc [https://www.mcs.anl.gov/petsc/] for linear/nonlinear solvers. 
IGAP4 is designed to solve (initial) boundary value problems numerically that involve fourth-order spatial derivatives in strong form and second-order in weak form; examples include the (transient) gradient elasticity and the Cahn-Hilliard equation. 

## Installation on Linux

1) Install PETSc 3.7+ [https://www.mcs.anl.gov/petsc/]. Set environmental variables PETSC_DIR and PETSC_ARCH as required by PETSc.

2) Set environmental variable, IGAP4_DIR, on the command line as:
> export IGAP4_DIR=/path/to/your/IGAP4/dir

3) Compile the source and make shared objects (this will take a few minutes):
> cd ${IGAP4_DIR}; make; make install

## Tests

> cd ${IGAP4_DIR}/example/ex1; make; mpiexec -n 8 ./main

## License
BSD 2-Clause License
