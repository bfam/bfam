include("bfam2d.jl")

using bfam2d
bfam = bfam2d
using MPI

MPI.Init()

bfam.init()


MPI.Finalized()


