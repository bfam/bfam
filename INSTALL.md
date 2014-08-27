# bfam

bfam is a set of tools to develop coupled discontinuous Galerkin
and multi-block summation-by-parts methods to solve multi-physics
problems.

## Install Guide

This document contains information on a how to install bfam

### the easy way!

Assuming that everything works as it ought:
    
    git clone git@ada.uc.nps.edu:bfam
    cd bfam
    mkdir build
    cd build
    cmake ../
    make
    make test

### Requirements
git 1.7 or newer
cmake 1.8 or newer
c and Fortran compiler (successfully built with gcc 4.7, 4.8; Intel 13)
an mpi distribution
optional: an opencl distribution

### Common problems
The problem that we have run into most often is cmake not finding the serial
and/or mpi compilers. If you run into this problem, the following definitions
can be added onto cmake with the '-D' option:
    
    CMAKE_C_COMPILER
    CMAKE_Fortran_COMPILER
    MPI_C_COMPILER
    MPI_Fortran_COMPILER
    MPIEXEC

### Using autotools to build bfam

There is initial support for an autotools build.  It should work like this:

    git clone some/repo/bfam.git
    cd bfam
    autotools/bootstrap
    mkdir build
    cd build
    ../configure CFLAGS="-std=c99 -Wall"
    make

Calling bootstrap is generally not required every time.
Typing make (or make install or make dist ...) should take care of everything.
