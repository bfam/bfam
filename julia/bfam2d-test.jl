include("bfam2d.jl")

using bfam2d
bfam = bfam2d
using MPI

MPI.Init()

bfam.init()

conn = bfam.connectivity_brick(2, 3, false, true)


MPI.Finalized()


