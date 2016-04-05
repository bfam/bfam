# Included Software

### crit-bit
We use a slightly modified version of Adam Langley's implementation of the
binary crit-bit from Dan Bernstein's qhasm.  The original file were
obtained from [here](https://github.com/agl/critbit.git).

### Lua
We include a version of Lua from LuaDist by Peter Drahoš and Peter Kapec
redistribution and use of the included files are allowed according to the terms
of the MIT license. Lua was written by Roberto Ierusalimschy, Waldemar Celes,
and Luiz Henrique de Figueiredo.

See the copyright statement below for license information of the included
Lua files.

    Copyright (C) 1994-2012 Lua.org, PUC-Rio.

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
    THE SOFTWARE.

### Catching Exceptions and Printing Stack Traces
The signal handler code is a modified version of the one presented in:

    http://spin.atomicobject.com/2013/01/13/exceptions-stack-traces-c/

by Job Vranish.  The code can be found here:

    https://gist.github.com/4441299

### zlib
The compression library zlib is included in bfam as a third-party library.
See the copyright and license information below for this library.

    (C) 1995-2012 Jean-loup Gailly and Mark Adler

    This software is provided 'as-is', without any express or implied
    warranty.  In no event will the authors be held liable for any damages
    arising from the use of this software.

    Permission is granted to anyone to use this software for any purpose,
    including commercial applications, and to alter it and redistribute it
    freely, subject to the following restrictions:

    1. The origin of this software must not be misrepresented; you must not
       claim that you wrote the original software. If you use this software
       in a product, an acknowledgment in the product documentation would be
       appreciated but is not required.
    2. Altered source versions must be plainly marked as such, and must not be
       misrepresented as being the original software.
    3. This notice may not be removed or altered from any source distribution.

    Jean-loup Gailly        Mark Adler
    jloup@gzip.org          madler@alumni.caltech.edu

### p4est
The octree library p4est is included in bfam as a third-party library.
Basic information about the library from its `README` file is:

    p4est is a C library to manage a collection (a forest) of multiple
    connected adaptive quadtrees or octrees in parallel.

    p4est is written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac
    and released under the GNU General Public Licence version 2 (or later),
    Copyright (C) 2010 The University of Texas System.

    The official web page for source code and documentation is www.p4est.org.
    Please send bug reports and ideas for contribution to info@p4est.org.

See the `INSTALL` file in `third_party/p4est-0.3.4.1.143-e573.tar.gz` for
more information.

### libocca

We include a version of libocca from:

    https://github.com/libocca/occa

with the command:

    git subtree pull --prefix=third_party/OCCA2 --squash git@github.com:libocca/occa.git master

libocca has  the license:

    The MIT License (MIT)

    Copyright (c) 2014 David Medina, Lucas Wilcox, and Tim Warburton

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
