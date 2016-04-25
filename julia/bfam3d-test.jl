include("bfam3d.jl")

using bfam3d
bfam = bfam3d
using MPI

MPI.Init()

bfam.init()


MPI.Finalized()


