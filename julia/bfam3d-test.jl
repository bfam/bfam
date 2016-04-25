include("bfam3d.jl")

using bfam3d
bfam = bfam3d
using MPI

MPI.Init()

bfam.init()

conn = bfam.connectivity_brick(2, 3, 4, false, true, false)

comm = MPI.COMM_WORLD
mesh_min_elements = 0
mesh_min_level = 0
domain = bfam.domain_pxest_new_ext_f(comm, conn, mesh_min_elements,
                                              mesh_min_level, 1)
bfam.domain_balance(domain)

MPI.Finalized()
