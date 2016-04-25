using MPI

const BFAM_LC_ALL     =  0
const BFAM_LC_ROOT    =  1
const BFAM_LL_DEFAULT = -1
const BFAM_LL_ALWAYS  =  0
const BFAM_LL_TRACE   =  1
const BFAM_LL_DEBUG   =  2
const BFAM_LL_VERBOSE =  3
const BFAM_LL_INFO    =  4
const BFAM_LL_WARNING =  5
const BFAM_LL_ERROR   =  6
const BFAM_LL_SILENT  =  7

const _bfam_functions = Dict{Symbol, ASCIIString}(
    :bfam_init_helper_f => "bfam_init_helper_f",
    )

const _p4est_functions = Dict{Symbol, ASCIIString}(
    :pxest_connectivity_new_brick => "p4est_connectivity_new_brick",
    )

const _p8est_functions = Dict{Symbol, ASCIIString}(
    :pxest_connectivity_new_brick => "p8est_connectivity_new_brick",
    )

function __init__()
  for (jname, fname) in _bfam_functions
    eval(:(const $jname = Libdl.dlsym(libbfam, $fname)))
  end
  if volumedim == 2
    for (jname, fname) in _p4est_functions
      eval(:(const $jname = Libdl.dlsym(libbfam, $fname)))
    end
  else volumedim == 3
    for (jname, fname) in _p8est_functions
      eval(:(const $jname = Libdl.dlsym(libbfam, $fname)))
    end
  end
end


# initialization functions
function init()
  init(STDOUT)
end
function init(stream::IO)
  @assert MPI.Initialized()
  verbosity = BFAM_LL_DEFAULT
  comm = MPI.COMM_WORLD
  init(stream, comm, verbosity)
end
function init(stream::IO, comm::MPI.Comm, verbosity::Int)
  @assert MPI.Initialized()
  loglevel = max(BFAM_LL_INFO - verbosity, BFAM_LL_ALWAYS);
  ccall(bfam_init_helper_f, Void,
        (Ptr{Void}, Cint, Cint),
        stream.handle, comm.val, loglevel)
end

# pxest connectivity
function connectivity_brick(nx::Int, ny::Int, px::Bool, py::Bool)
  @assert volumedim == 2
  conn = ccall(pxest_connectivity_new_brick,
               Ptr{Void}, (Cint, Cint, Cint, Cint), nx, ny, px, py)
  return conn
end
function connectivity_brick(nx::Int , ny::Int , nz::Int,
                            px::Bool, py::Bool, pz::Bool)
  @assert volumedim == 3
  conn = ccall(pxest_connectivity_new_brick,
               Ptr{Void}, (Cint, Cint, Cint, Cint, Cint, Cint),
               nx, ny, nz, px, py, pz)

end
