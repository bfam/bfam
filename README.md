# bfam

bfam is a set of tools to develop coupled, multi-physics discontinuous Galerkin
solvers on adapted Cartesian meshes. bfam tools are not a standalone solver, but
designed to be integrated with an application.

# Third Party Dependencies
bfam depends strongly on [p4est][1] for mesh adaptivity.

bfam has weaker (optional) dependence on [lua][2] and [occa][3].

Other included software which is included in bfam can be seen in the [TPL][4]
document

[1]: http://www.p4est.org.
[2]: https://www.lua.org/
[3]: https://github.com/libocca/occa
[4]: https://github.com/bfam/bfam/blob/master/TPL.md
