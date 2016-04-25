include("bfam3d.jl")

using bfam3d
bfam = bfam3d
using MPI

MPI.Init()

bfam.init()

conn = bfam.connectivity_brick(2, 3, 4, false, true, false)


MPI.Finalized()
