module bfam2d
const volumedim = 2
const libbfam = Libdl.dlopen("../libbfam2d.so")

include("bfam-base.jl")
end # module bfam2d
