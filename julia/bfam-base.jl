using MPI

type Domain
    domain_ptr::Ptr{Void}
    pxest_ptr::Ptr{Void}
    Domain(domain_ptr, pxest_ptr) = new(domain_ptr, pxest_ptr)
end

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
    :bfam_domain_pxest_new_ext_f => "bfam_domain_pxest_new_ext_f",
    :bfam_domain_get_pxest => "bfam_domain_get_pxest",
    :bfam_domain_pxest_balance => "bfam_domain_pxest_balance",
    )

const _p4est_functions = Dict{Symbol, ASCIIString}(
    :pxest_connectivity_new_brick => "p4est_connectivity_new_brick",
    :pxest_vtk_write_all => "p4est_vtk_write_all",
    )

const _p8est_functions = Dict{Symbol, ASCIIString}(
    :pxest_connectivity_new_brick => "p8est_connectivity_new_brick",
    :pxest_vtk_write_all => "p8est_vtk_write_all",
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

# domain_new
function domain_pxest_new_ext_f(comm, pxest_conn, min_quadrants,
                                min_level, fill_uniform)
  domain_ptr = ccall(bfam_domain_pxest_new_ext_f,
                Ptr{Void},  # bfam_domain_pxest_t *
                (Cint,      # int f_domcomm,
                 Ptr{Void}, # p4est_connectivity_t *conn,
                 Cint,      # int min_quadrants,
                 Cint,      # int min_level,
                 Cint),     # int fill_uniform);
                comm.val, pxest_conn, min_quadrants, min_level, fill_uniform)
  pxest_ptr = ccall(bfam_domain_get_pxest,
                Ptr{Void},    # pxest_t *
                (Ptr{Void},), # bfam_domain_pxest_t *domain
                domain_ptr)
  return Domain(domain_ptr, pxest_ptr)
end

function domain_balance(domain::Domain)
  ccall(bfam_domain_pxest_balance,
                Void,
                (Ptr{Void},),
                domain.domain_ptr)
end

function domain_pxest_write(domain::Domain, output)
  ccall(pxest_vtk_write_all,
        Void, (Ptr{Void}, # p4est_t * p4est,
               Ptr{Void}, # p4est_geometry_t * geom,
               Cdouble,   # double         scale,
               Cint,      # int write_tree,
               Cint,      # int write_level,
               Cint,      # int write_rank,
               Cint,      # int wrap_rank,
               Cint,      # int num_scalars,
               Cint,      # int num_vectors,
               Cstring),  # const char * filename,
          domain.pxest_ptr, C_NULL, 1, 1, 1, 1, 0, 0, 0, output)
end
