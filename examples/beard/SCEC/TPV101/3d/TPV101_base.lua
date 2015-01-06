tout
min_level = 0
max_level = 0
output_prefix = "TPV101"
data_directory = "data"
-- connectivity info
connectivity = "brick"
brick =
{
  nx = 2,
  ny = 2,
  nz = 2,
  periodic_x = 0,
  periodic_y = 0,
  periodic_z = 0,
}
default_boundary_tag = "non-reflecting"

-- set up the domain
Lx = 50
Ly = 50
Lz = 50

Cx = brick.nx/2
Cy = brick.ny/2
Cz = brick.nz/2

function connectivity_vertices(x, y, z)
  xout = Lx*(x-Cx)
  yout = Ly*(y-Cy)
  zout = Lz*(z-Cz)
  return xout,yout,zout
end

function refinement_function(
  x0,y0,z0,x1,y1,z1,
  x2,y2,z2,x3,y3,z3,
  x4,y4,z4,x5,y5,z5,
  x6,y6,z6,x7,y7,z7,
  level, treeid)

  xa0 = math.abs(x0)
  xa1 = math.abs(x1)
  xa2 = math.abs(x2)
  xa3 = math.abs(x3)
  xmin = math.min( xa0,xa1)
  xmin = math.min(xmin,xa2)
  xmin = math.min(xmin,xa2)

  ya0 = math.abs(y0)
  ya1 = math.abs(y1)
  ya2 = math.abs(y2)
  ya3 = math.abs(y3)
  ymin = math.min( ya0,ya1)
  ymin = math.min(ymin,ya2)
  ymin = math.min(ymin,ya3)

  zmin = math.min(z0,z1,z2,z3)

  if level < min_level then
    return 1
  elseif level >= max_level then
    return 0
  elseif xmin <= 15 and zmin >= -15 and ymin <= 0.1 then
    return 1
  else
    return 0
  end

  -- if level < min_level then
  --   return 1
  -- elseif level >= max_level then
  --   return 0
  -- else
  --   return 0
  -- end
end

function element_order(
  x0,y0,z0,x1,y1,z1,
  x2,y2,z2,x3,y3,z3,
  x4,y4,z4,x5,y5,z5,
  x6,y6,z6,x7,y7,z7,
  level, treeid)

  -- N = treeid%3+1
  N = 5

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

tend   = 12
tout   = 1
tfout  = 0.1
tdisp  = tout
tstations  = 0.01
nerr   = 0

function time_step_parameters(dt)
  dt      = 0.5*dt
  nfoutput = math.ceil(tfout / dt)
  dt       = tfout / nfoutput

  noutput    = tout  / dt
  ndisp      = tdisp / dt
  nsteps     = tend / dt
  nstations  = tstations / dt

  return dt,nsteps, ndisp, noutput, nfoutput, nstations
end

fault_stations = {
  "faultst-090dp075",    -9.0,  7.5, 0.0,   0,0, 1.0,   0.1,
  "faultst-120dp030",   -12.0,  3.0, 0.0,   0,0, 1.0,   0.1,
  "faultst-120dp120",   -12.0, 12.0, 0.0,   0,0, 1.0,   0.1,
  "faultst000dp030" ,     0.0,  3.0, 0.0,   0,0, 1.0,   0.1,
  "faultst000dp075" ,     0.0,  7.5, 0.0,   0,0, 1.0,   0.1,
  "faultst000dp120" ,     0.0, 12.0, 0.0,   0,0, 1.0,   0.1,
  "faultst090dp075" ,     9.0,  7.5, 0.0,   0,0, 1.0,   0.1,
  "faultst120dp030" ,    12.0,  3.0, 0.0,   0,0, 1.0,   0.1,
  "faultst120dp120" ,    12.0, 12.0, 0.0,   0,0, 1.0,   0.1,
}

-- faults
fault = {
  type   = "friction",
  tag    = "ageing law",
  f0     = 0.6,
  v0     = 1e-6,
  a      = "a_function",
  b      = 0.012,
  L      = 0.02,
  S33_0  = -120.0,
  S23_0  =   75.0,
  psi0   = "psi0_function"
}

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
  fault,
}

-- friction stuff
glueid_treeid_faceid = {
  3, 0, 5,
  3, 1, 5,
  3, 2, 5,
  3, 3, 5,
  3, 4, 4,
  3, 5, 4,
  3, 6, 4,
  3, 7, 4,
}
