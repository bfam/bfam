tanh = math.tanh
sinh = math.sinh
abs  = math.abs
sqrt = math.sqrt
min  = math.min
ln   = math.log
exp  = math.exp
ceil = math.ceil

-- refinement parameters
min_level = 4
max_level = 4
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
  yout = Ly*(y-Cy)+7.5
  zout = Lz*(z-Cz)
  return xout,yout,zout
end

function refinement_function(
  x0,y0,z0,x1,y1,z1,
  x2,y2,z2,x3,y3,z3,
  x4,y4,z4,x5,y5,z5,
  x6,y6,z6,x7,y7,z7,
  level, treeid)

  xa0 = abs(x0)
  xa1 = abs(x1)
  xa2 = abs(x2)
  xa3 = abs(x3)
  xmin = min( xa0,xa1)
  xmin = min(xmin,xa2)
  xmin = min(xmin,xa2)

  ya0 = abs(y0)
  ya1 = abs(y1)
  ya2 = abs(y2)
  ya3 = abs(y3)
  ymin = min( ya0,ya1)
  ymin = min(ymin,ya2)
  ymin = min(ymin,ya3)

  zmin = min(z0,z1,z2,z3)

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
  N = 3

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
  nfoutput = ceil(tfout / dt)
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
function smooth_boxcar(r,w,W)
  r = abs(r)
  if r < W then
    return 1
  elseif r > W+w then
    return 0
  else
    return 0.5*(1+tanh(w/(r-W-w)+w/(r-W)))
  end
end

w = 3
W = 15
x0 = 0
y0 = 7.5
function a_function(x,y,z,t)
  Bx = smooth_boxcar(x-x0,w,    W)
  By = smooth_boxcar(y-y0,w,0.5*W)
  return 0.008*(2 - Bx*By)
end

function psi_function(x,y,z,t)
  f0   = 0.6
  V0   = 1e-6
  Vini = 1e-12
  a    = a_function(x,y,z,t)
  b    = 0.012
  L    = 0.02
  N    = 120.0
  T    =   75.0
  psi  = (L/V0)*exp((a*ln(2*sinh(T/(a*N))) -f0 - a*ln(Vini/V0))/b)
  return psi
end

R = 3.0
function S23_function(x,y,z,t)
  S23 = 75
  r = sqrt((x-x0)^2+(y-y0)^2)
  if r < R then
    S23 = S23 + 25*exp(r^2/(r^2-R^2))
  end
  return S23
end

fault = {
  type   = "friction",
  tag    = "ageing law",
  f0     = 0.6,
  V0     = 1e-6,
  a      = "a_function",
  b      = 0.012,
  L      = 0.02,
  S33_0  = -120.0,
  S23_0  = "S23_function",
  psi    = "psi_function"
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
  bc_rigid,
  fault,
}

-- friction stuff
glueid_treeid_faceid = {
  4, 0, 5,
  4, 1, 5,
  4, 2, 5,
  4, 3, 5,
  4, 4, 4,
  4, 5, 4,
  4, 6, 4,
  4, 7, 4,
}
