abs  = math.abs

-- refinement parameters
min_level = 0
max_level = 6
output_prefix  = "curve"
data_directory = "data"

-- connectivity info
connectivity = "brick"
brick =
{
  nx = 3,
  ny = 2,
  nz = 2,
  periodic_x = 0,
  periodic_y = 0,
  periodic_z = 0,
  bc0 = 1,
  bc1 = 1,
  bc2 = 1,
  bc3 = 1,
  bc4 = 1,
  bc5 = 1,
}

-- set up the domain
L = 30

Cx = brick.nx/2
Cy = brick.ny/2
Cz = brick.nz

function connectivity_vertices(xin, yin, zin)
  x = L*(xin-Cx)
  y = L*(yin-Cy)
  z = L*(zin-Cz)

  return x,y,z
end

function transform_nodes(x, y)

  if abs(x) < L/2 then
    -- y = y + ((x/10)^2 - (L/10)^2) * (1-abs(y)/L)
    y = y + 1*((x/(L/2))^2 - 1) * (1-abs(y)/L)
  end
  -- if abs(z) < L/2 then
  --   z = z + ((y/10)^2 - (L/10)^2)
  -- end

  return x,y,z
end

function refinement_function(
  x0,y0,z0,x1,y1,z1,
  x2,y2,z2,x3,y3,z3,
  -- x4,y4,z4,x5,y5,z5,
  -- x6,y6,z6,x7,y7,z7,
  level, treeid)

  if level < max_level then
    return 1
  end

  return 0
end

function element_order(
  x0,y0,z0,x1,y1,z1,
  x2,y2,z2,x3,y3,z3,
  x4,y4,z4,x5,y5,z5,
  x6,y6,z6,x7,y7,z7,
  level, treeid)

  N = 6

  return N
end

-- material properties
cs = 3.464
cp = 6
rho = 2.670
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

tend   = 13
tout   = 1
tfout  = 0.1
tdisp  = 0.1
tstations  = 0.01
nerr   = 0

function time_step_parameters(dt)
  dt        = 0.5*dt

  nstations = math.ceil(tstations / dt)
  dt        = tstations / nstations

  ndisp      = tdisp / dt
  nsteps     = tend / dt
  noutput    = tout / dt
  nfoutput   = tfout / dt
  nstations  = tstations / dt

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


nuc_width = 0.1
function S12_fault(x,y,z,t)
  S12_0 = 29.38
  r = math.abs(x)-1.5
  v = 1 - 1/(1+math.exp(-2*r/nuc_width))

  return S12_0 + v*11.6
end
fault = {
  type   = "friction",
  tag    = "slip weakening",
  fs     =    0.677,
  fd     =    0.373,
  Dc     =    0.4,
  S12_0  =   "S12_fault",
  S11_0  =   -60,
  S22_0  =   -60,
}

glue_info = {
  bc_nonreflect,
  bc_free,
  bc_rigid,
  fault,
}

glueid_treeid_faceid = {
  4, 1, 3,
  4, 4, 2,
}
