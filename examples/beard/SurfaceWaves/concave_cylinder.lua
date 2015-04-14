pi   = math.pi
sqrt = math.sqrt
abs  = math.abs
cos = math.cos
sin = math.sin

output_prefix = "concave_cylinder"
data_directory = "data"
elem_order = 12
max_level = 0
min_level = 0
static_refinement=3

mesh_file = "/home/jekozdon/codes/bfam/examples/beard/SurfaceWaves/concave_cylinder.inp"
glueid_treeid_faceid = mesh_file

default_boundary_tag = "free surface"

c1 = sqrt(pi/5)/2
c2 = 1/sqrt(2)
c3 = 3

function connectivity_vertices(xin, yin, zin)

  z = 0

  if abs(xin) == 1 then
    x = (xin/abs(xin))*c1
  elseif abs(xin) == 2 then
    x = (xin/abs(xin))*c2
  elseif abs(xin) == 3 then
    x = (xin/abs(xin))*c3
  end

  if abs(yin) == 1 then
    y = (yin/abs(yin))*c1
  elseif abs(yin) == 2 then
    y = (yin/abs(yin))*c2
  elseif abs(yin) == 3 then
    y = (yin/abs(yin))*c3
  end

  return x, y, z
end

function transform_nodes(x, y, z)
  local x1 = x
  local y1 = y
  local x2 = x
  local y2 = y
  local  s = 0

  if     abs(x) <= c2 and y >= abs(x) then
    -- NORTH
    s  = (y-c2)/(c3-c2)
    y1 = sqrt(1-x^2)
    y2 = c3
  elseif abs(x) <= c2 and y <= -abs(x) then
    -- SOUTH
    s  = (-y-c2)/(c3-c2)
    y1 = -sqrt(1-x^2)
    y2 = -c3
  elseif abs(y) <= c2 and x <= -abs(y) then
    -- WEST
    s  = (-x-c2)/(c3-c2)
    x1 = -sqrt(1-y^2)
    x2 = -c3
  elseif abs(y) <= c2 and x >= abs(y) then
    -- EAST
    s  = (x-c2)/(c3-c2)
    x1 = sqrt(1-y^2)
    x2 = c3
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
mu  = 1
lam = 1
rho = 1
cp  = sqrt((lam+2*mu)/rho)
cs  = sqrt(mu/rho)

-- time stepper to use
lsrk_method  = "KC54"

tend   = 2*pi
tend   = 2
tout   = tend
tdisp  = tend
-- nerr   = 0
function nerr(dt)
  return terr/dt
end
terr  = tend

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

  nsteps = 1
  noutput = 1

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

bc_exact = {
  type = "boundary",
  tag  = "user defined",
}

glue_info = {
  bc_nonreflect,
  bc_free,
  bc_rigid,
  bc_exact,
}

-- field conditions
function v1(x,y,z,t)
  return 0
end

function v2(x,y,z,t)
  return 0
end

function S11(x,y,z,t)
  return 0
end

function S22(x,y,z,t)
  return 0
end

function S12(x,y,z,t)
  return 0
end

function v3(x,y,z,t)
  return 0
end

function S33(x,y,z,t)
  return 0
end

function S23(x,y,z,t)
  return 0
end

function S13(x,y,z,t)
  return 0
end
