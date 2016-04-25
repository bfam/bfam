module bfam3d
const volumedim = 3
const libbfam = Libdl.dlopen("../libbfam3d.so")

include("bfam-base.jl")
end # module bfam3d
