pi   = math.pi
sqrt = math.sqrt
abs  = math.abs

output_prefix = "convex_cylinder"
data_directory = "data"
elem_order = 4
max_level = 600
min_level = 5

mesh_file = "/home/jekozdon/codes/bfam/examples/beard/SurfaceWaves/convex_cylinder.inp"

default_boundary_tag = "free surface"

c1 = sqrt(pi/5)/2
c2 = 1/sqrt(2)

function connectivity_vertices(xin, yin, zin)

  z = 0

  if abs(xin) == 1 then
    x = (xin/abs(xin))*c1
    y = (yin/abs(yin))*c1
  elseif abs(xin) == 2 then
    x = (xin/abs(xin))*c2
    y = (yin/abs(yin))*c2
  end

  return x, y, z
end

function transform_nodes(x, y, z)
  local x1 = x
  local y1 = y
  local  s = 0
  local  r = sqrt(x^2+y^2)
  local x2 = 0
  local y2 = 0

  if     abs(x) <= abs(y) and y >= c1 then
    -- NORTH
    s = (y-c1)/(c2-c1)

    x1 = (c1/y)*x
    y1 = c1

    x2 = x/r
    y2 = y/r
  elseif abs(x) <= abs(y) and y <= -c1 then
    -- SOUTH
    s = (abs(y)-c1)/(c2-c1)

    x1 = -(c1/y)*x
    y1 = -c1

    x2 = x/r
    y2 = y/r
  elseif abs(y) <= abs(x) and x >=  c1 then
    -- WEST
    s = (x-c1)/(c2-c1)

    x1 = c1
    y1 = (c1/x)*y

    x2 = x/r
    y2 = y/r
  elseif abs(y) <= abs(x) and x <= -c1 then
    -- EAST
    s = (abs(x)-c1)/(c2-c1)

    x1 = -c1
    y1 = -(c1/x)*y

    x2 = x/r
    y2 = y/r
  end

  x = (1-s)*x1 + s*x2
  y = (1-s)*y1 + s*y2

  return x, y, z
end

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

  return sqrt(hmin), sqrt(hmax)

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
