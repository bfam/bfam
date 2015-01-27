cos = math.cos
sin = math.sin
output_prefix = "convex_cylinder"
data_directory = "data"
elem_order = 4
max_level = 600
min_level = 0

mesh_file = "/home/jekozdon/codes/bfam/examples/beard/SurfaceWaves/convex_cylinder.inp"

default_boundary_tag = "free surface"

-- ABQ->p4est:
-- 1 -> 4
-- 2 -> 5
-- 3 -> 2
-- 4 -> 1
-- 5 -> 3
-- 6 -> 0

function element_size(
  x0,y0,z0,x1,y1,z1,
  x2,y2,z2,x3,y3,z3)


  h = (x0-x1)^2+(y0-y1)^2+(z0-z1)^2
  hmin = h
  hmax = h

  h = (x0-x2)^2+(y0-y2)^2+(z0-z2)^2
  hmax = math.max(hmax,h)
  hmin = math.min(hmin,h)

  h = (x1-x3)^2+(y1-y3)^2+(z1-z3)^2
  hmax = math.max(hmax,h)
  hmin = math.min(hmin,h)

  h = (x2-x3)^2+(y2-y3)^2+(z2-z3)^2
  hmax = math.max(hmax,h)
  hmin = math.min(hmin,h)

  return math.sqrt(hmin), math.sqrt(hmax)

end

function refinement_function(
  x0,y0,z0,x1,y1,z1,
  x2,y2,z2,x3,y3,z3,
  level, treeid)

  if(level < min_level) then
    return 1
  end

  return 0
end

function element_order(
  x0,y0,z0,x1,y1,z1,
  x2,y2,z2,x3,y3,z3,
  level, treeid)

  return elem_order
end

-- material properties
cs = 1
cp = 2
rho = 1
mu  = rho*cs^2
lam = rho*cp^2-2*mu

-- field conditions
S11 = 0
S22 = 0
S33 = 0
S12 = 0
S13 = 0
S23 = 0
v1  = 0
v2  = 0
v3  = 0

-- time stepper to use
lsrk_method  = "KC54"

tend   = 1
tout   = 1
tdisp  = 0.01
nerr   = 0

T = tout

function time_step_parameters(dt)
  dt = 0.5*dt
  N  = math.ceil(T / dt)
  dt = T / N

  nsteps     = tend / dt
  ndisp      = tdisp / dt
  noutput    = tout  / dt
  nstations  = -1
  nfoutput   = -1

  return dt,nsteps, ndisp, noutput, nfoutput, nstations
end

bc_free = {
  type = "boundary",
  tag  = "free surface",
}

bc_nonreflect = {
  type = "boundary",
  tag  = "non-reflecting",
}

bc_rigid = {
  type = "boundary",
  tag  = "rigid wall",
}

glue_info = {
  bc_nonreflect,
  bc_free,
  bc_rigid,
}
