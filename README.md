# bfam

bfam is a set of tools to develop coupled discontinuous Galerkin
and multi-block summation-by-parts methods to solve multi-physics
problems.

## Included Software

### `cmake/findopencl`
CMake utilities to find OpenCL libs.  The code was
added to the repository using:

    git subtree add --prefix=cmake/findopencl --squash \
        git://gitorious.org/findopencl/findopencl.git master

and can be updated to the latest upstream version using:

    git subtree merge --prefix=cmake/findopencl --squash \
        git://gitorious.org/findopencl/findopencl.git master

### gopt
We use a slightly modified version Tom Vajzovic's Gopt/8.1 library which
is released in the public domain.  The original files were obtained from
[here](http://www.purposeful.co.uk/software/gopt/).
