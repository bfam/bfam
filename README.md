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
We use a slightly modified version of Tom Vajzovic's Gopt/8.1 library which
is released in the public domain.  The original files were obtained from
[here](http://www.purposeful.co.uk/software/gopt/).

### crit-bit
We use a slightly modified version of Adam Langley's implementation of the
binary crit-bit from Dan Bernstein's qhasm.  The original file were
obtained from [here](https://github.com/agl/critbit.git).

### opencl-helper
The OpenCL helper functions in `src/bfam_opencl.{h,c}` are modified
versions of the ones found in Andreas Kloeckner's cl-helper routines
obtained from the git repo:

    git@github.com:hpc12/lec1-demo

See the copyright statement below for license information of the included
files.

    Copyright (c) 2010, 2012 Andreas Kloeckner

    Permission is hereby granted, free of charge, to any person obtaining
    a copy of this software and associated documentation files (the
    "Software"), to deal in the Software without restriction, including
    without limitation the rights to use, copy, modify, merge, publish,
    distribute, sublicense, and/or sell copies of the Software, and to
    permit persons to whom the Software is furnished to do so, subject to
    the following conditions:

    The above copyright notice and this permission notice shall be
    included in all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
    IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
    CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
    TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
    SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
