max = math.max
min = math.min
abs = math.abs
sqrt = math.sqrt
-- refinement parameters
max_level = 300
height_split = 5
output_prefix = "TPV32_base"
data_directory = "data"
elem_order = 4
h1_targ   = elem_order*0.1
r_targ    = 5
D_targ    = 20
hmax_targ = elem_order*2
-- hmax_targ = h1_targ

-- connectivity info
connectivity = "brick"
BUF = 2
brick =
{
  nx = 2+2*BUF,
  ny = 1+1*BUF,
  nz = 2+2*BUF,
  periodic_x = 0,
  periodic_y = 0,
  periodic_z = 0,
  bc0 = 1,
  bc1 = 1,
  bc2 = 2,
  bc3 = 1,
  bc4 = 1,
  bc5 = 1,
}

-- set up the domain
Lx = 15
Ly = 15
Lz = 15

Cx = brick.nx/2
Cy = 0
Cz = brick.nz/2

function connectivity_vertices(x, y, z)
  xout = Lx*(x-Cx)
  yout = Ly*(y-Cy)
  zout = Lz*(z-Cz)
  return xout,yout,zout
end



--REFINEMENT FUNCTION
function fault_distance(x,y,z)
  xf = max(0,abs(x)-15)
  zf = abs(z)
  yf = max(0,abs(y)-15)
  r  = xf^2 + yf^2 + zf^2
  return sqrt(r)
end

function element_size(
  x0,y0,z0,x1,y1,z1,
  x2,y2,z2,x3,y3,z3,
  x4,y4,z4,x5,y5,z5,
  x6,y6,z6,x7,y7,z7)


  h = (x0-x1)^2+(y0-y1)^2+(z0-z1)^2
  hmin = h
  hmax = h

  h = (x0-x2)^2+(y0-y2)^2+(z0-z2)^2
  hmax = max(hmax,h)
  hmin = min(hmin,h)

  h = (x0-x4)^2+(y0-y4)^2+(z0-z4)^2
  hmax = max(hmax,h)
  hmin = min(hmin,h)

  h = (x1-x3)^2+(y1-y3)^2+(z1-z3)^2
  hmax = max(hmax,h)
  hmin = min(hmin,h)

  h = (x1-x5)^2+(y1-y5)^2+(z1-z5)^2
  hmax = max(hmax,h)
  hmin = min(hmin,h)

  h = (x2-x3)^2+(y2-y3)^2+(z2-z3)^2
  hmax = max(hmax,h)
  hmin = min(hmin,h)

  h = (x2-x6)^2+(y2-y6)^2+(z2-z6)^2
  hmax = max(hmax,h)
  hmin = min(hmin,h)

  h = (x3-x7)^2+(y3-y7)^2+(z3-z7)^2
  hmax = max(hmax,h)
  hmin = min(hmin,h)

  h = (x4-x5)^2+(y4-y5)^2+(z4-z5)^2
  hmax = max(hmax,h)
  hmin = min(hmin,h)

  h = (x4-x6)^2+(y4-y6)^2+(z4-z6)^2
  hmax = max(hmax,h)
  hmin = min(hmin,h)

  h = (x5-x7)^2+(y5-y7)^2+(z5-z7)^2
  hmax = max(hmax,h)
  hmin = min(hmin,h)

  h = (x6-x7)^2+(y6-y7)^2+(z6-z7)^2
  hmax = max(hmax,h)
  hmin = min(hmin,h)

  return sqrt(hmin), sqrt(hmax)

end

function refinement_function(
  x0,y0,z0,x1,y1,z1,
  x2,y2,z2,x3,y3,z3,
  x4,y4,z4,x5,y5,z5,
  x6,y6,z6,x7,y7,z7,
  level, treeid)

  if(level >= max_level) then
    return 0
  end

  r = fault_distance(x0,y0,z0)
  r = min(r,fault_distance(x1,y1,z1))
  r = min(r,fault_distance(x2,y2,z2))
  r = min(r,fault_distance(x3,y3,z3))
  r = min(r,fault_distance(x4,y4,z4))
  r = min(r,fault_distance(x5,y5,z5))
  r = min(r,fault_distance(x6,y6,z6))
  r = min(r,fault_distance(x7,y7,z7))

  hmin, hmax = element_size( x0,y0,z0,x1,y1,z1,
                             x2,y2,z2,x3,y3,z3,
                             x4,y4,z4,x5,y5,z5,
                             x6,y6,z6,x7,y7,z7)

  D = (min(D_targ,max(r_targ,r))-r_targ)/(D_targ-r_targ)
  if hmax > h1_targ*(1-D) + hmax_targ*D then
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

  return elem_order
end

-- material properties
material = {
 D   = {    0, 0.5  , 1.0  , 1.6  , 2.4  , 3.6  , 5.0  , 9.0  , 11.0  , 15.0  ,},
 cp  = {2.200, 3.000, 3.600, 4.400, 4.800, 5.250, 5.500, 5.750,  6.100,  6.300,},
 cs  = {1.050, 1.400, 1.950, 2.500, 2.800, 3.100, 3.250, 3.450,  3.600,  3.700,},
 rho = {2.200, 2.450, 2.550, 2.600, 2.600, 2.620, 2.650, 2.720,  2.750,  2.900,},
}

function transform_nodes(x, y, z)
  Hy = Ly / 2^height_split
  km = 0
  yb = 0
  for key,D in pairs(material.D) do
    kp = math.floor(D/Hy+0.5)
    if Hy*km < y and y <= Hy*kp then
      dy = (y-km*Hy)/(kp*Hy-km*Hy)
      y = yb + dy*(D-yb)
      break
    end
    km = kp
    yb = D
  end
  return x,y,z
end

function cp(x,y,z,t)
  k = 1
  while material.D[k] do
    if material.D[k] > y then
      dy = (material.D[k-1]-y)/(material.D[k-1]-material.D[k])
      V  = material.cp[k-1] - (material.cp[k-1]-material.cp[k])*dy
      break
    end
    V = material.cp[k]
    k = k + 1
  end
  return V
end
function cs(x,y,z,t)
  k = 1
  while material.D[k] do
    if material.D[k] > y then
      dy = (material.D[k-1]-y)/(material.D[k-1]-material.D[k])
      V  = material.cs[k-1] - (material.cs[k-1]-material.cs[k])*dy
      break
    end
    V = material.cs[k]
    k = k + 1
  end
  return V
end
function rho(x,y,z,t)
  k = 1
  while material.D[k] do
    if material.D[k] > y then
      dy = (material.D[k-1]-y)/(material.D[k-1]-material.D[k])
      V  = material.rho[k-1] - (material.rho[k-1]-material.rho[k])*dy
      break
    end
    V = material.rho[k]
    k = k + 1
  end
  return V
end
function mu(x,y,z,t)
  return rho(x,y,z,t)*cs(x,y,z,t)^2
end
function lam(x,y,z,t)
  return rho(x,y,z,t)*cp(x,y,z,t)^2-2*mu(x,y,z,t)
end

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

tend       = 15
-- tout       = 1
-- tfout      = 0.01
tout       = -1
tfout      = -1
tdisp      = 0.01
tstations  = 0.01
nerr       = 0

function time_step_parameters(dt)
  dt      = 0.5*dt

  T       = tstations
  n       = math.ceil(T / dt)
  dt      = T / n

  noutput    = tout      / dt
  ndisp      = tdisp     / dt
  nsteps     = tend      / dt
  nstations  = tstations / dt
  nfoutput   = tfout     / dt
  -- nsteps     = 1
  -- noutput    = 1
  -- nfoutput   = 1

  return dt,nsteps, ndisp, noutput, nfoutput, nstations
end

fault_stations = {
   "faultst000dp000", 00.0, 00.0, 0.0, 0.0, 0.0, 1.0, 0.1,
   "faultst000dp002", 00.0, 00.2, 0.0, 0.0, 0.0, 1.0, 0.1,
   "faultst000dp005", 00.0, 00.5, 0.0, 0.0, 0.0, 1.0, 0.1,
   "faultst000dp010", 00.0, 01.0, 0.0, 0.0, 0.0, 1.0, 0.1,
   "faultst000dp024", 00.0, 02.4, 0.0, 0.0, 0.0, 1.0, 0.1,
   "faultst000dp030", 00.0, 03.0, 0.0, 0.0, 0.0, 1.0, 0.1,
   "faultst000dp050", 00.0, 05.0, 0.0, 0.0, 0.0, 1.0, 0.1,
   "faultst000dp075", 00.0, 07.5, 0.0, 0.0, 0.0, 1.0, 0.1,
   "faultst000dp100", 00.0, 10.0, 0.0, 0.0, 0.0, 1.0, 0.1,
   "faultst000dp120", 00.0, 12.0, 0.0, 0.0, 0.0, 1.0, 0.1,
   "faultst060dp000", 06.0, 00.0, 0.0, 0.0, 0.0, 1.0, 0.1,
   "faultst060dp002", 06.0, 00.2, 0.0, 0.0, 0.0, 1.0, 0.1,
   "faultst060dp005", 06.0, 00.5, 0.0, 0.0, 0.0, 1.0, 0.1,
   "faultst060dp010", 06.0, 01.0, 0.0, 0.0, 0.0, 1.0, 0.1,
   "faultst060dp024", 06.0, 02.4, 0.0, 0.0, 0.0, 1.0, 0.1,
   "faultst060dp030", 06.0, 03.0, 0.0, 0.0, 0.0, 1.0, 0.1,
   "faultst060dp050", 06.0, 05.0, 0.0, 0.0, 0.0, 1.0, 0.1,
   "faultst060dp075", 06.0, 07.5, 0.0, 0.0, 0.0, 1.0, 0.1,
   "faultst060dp100", 06.0, 10.0, 0.0, 0.0, 0.0, 1.0, 0.1,
   "faultst060dp120", 06.0, 12.0, 0.0, 0.0, 0.0, 1.0, 0.1,
   "faultst120dp000", 12.0, 00.0, 0.0, 0.0, 0.0, 1.0, 0.1,
   "faultst120dp002", 12.0, 00.2, 0.0, 0.0, 0.0, 1.0, 0.1,
   "faultst120dp005", 12.0, 00.5, 0.0, 0.0, 0.0, 1.0, 0.1,
   "faultst120dp010", 12.0, 01.0, 0.0, 0.0, 0.0, 1.0, 0.1,
   "faultst120dp024", 12.0, 02.4, 0.0, 0.0, 0.0, 1.0, 0.1,
   "faultst120dp030", 12.0, 03.0, 0.0, 0.0, 0.0, 1.0, 0.1,
   "faultst120dp050", 12.0, 05.0, 0.0, 0.0, 0.0, 1.0, 0.1,
   "faultst120dp075", 12.0, 07.5, 0.0, 0.0, 0.0, 1.0, 0.1,
   "faultst120dp100", 12.0, 10.0, 0.0, 0.0, 0.0, 1.0, 0.1,
   "faultst120dp120", 12.0, 12.0, 0.0, 0.0, 0.0, 1.0, 0.1,
}

volume_stations = {
  "body-030st000dp000",  0.0, 0.0,  -3.0,
  "body-030st000dp005",  0.0, 0.5,  -3.0,
  "body-030st000dp024",  0.0, 2.4,  -3.0,
  "body-030st120dp000", 12.0, 0.0,  -3.0,
  "body-030st120dp005", 12.0, 0.5,  -3.0,
  "body-030st120dp024", 12.0, 2.4,  -3.0,
  "body-090st000dp000",  0.0, 0.0,  -9.0,
  "body-090st000dp005",  0.0, 0.5,  -9.0,
  "body-090st000dp024",  0.0, 2.4,  -9.0,
  "body-090st150dp000", 15.0, 0.0,  -9.0,
  "body-090st150dp005", 15.0, 0.5,  -9.0,
  "body-090st150dp024", 15.0, 2.4,  -9.0,
  "body-150st000dp000",  0.0, 0.0, -15.0,
  "body-150st000dp005",  0.0, 0.5, -15.0,
  "body-150st000dp024",  0.0, 2.4, -15.0,
  "body-150st150dp000", 15.0, 0.0, -15.0,
  "body-150st150dp005", 15.0, 0.5, -15.0,
  "body-150st150dp024", 15.0, 2.4, -15.0,
}

-- faults
mu0  = 2.670*3.464^2
function S11_0_function(x,y,z,t)
  return -60 * mu(x,y,z,t) / mu0
end
function S33_0_function(x,y,z,t)
  return -60 * mu(x,y,z,t) / mu0
end
function S13_0_function(x,y,z,t)
  r = (min(2,max(1.4,sqrt(x^2 + (y-7.5)^2)))-1.4)/0.6
  return (30 + 2.475*(1+math.cos(math.pi*r)))*mu(x,y,z,t) / mu0
end
function c0_function(x,y,z,t)
 return 0.000425*max(0,2.4-y)
end

fault = {
  type   = "friction",
  tag    = "slip weakening",
  fs     = 0.58,
  fd     = 0.45,
  Dc     = 0.18,
  S11_0  = "S11_0_function",
  S12_0  = 0,
  S13_0  = "S13_0_function",
  S22_0  = 0,
  S23_0  = 0,
  S33_0  = "S33_0_function",
  c0     = "c0_function",
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
  4, (Cx-1) + (0)*brick.nx + (Cz-1)*brick.nx*brick.ny, 5,
  4, (Cx-1) + (0)*brick.nx + (Cz  )*brick.nx*brick.ny, 4,
  4, (Cx  ) + (0)*brick.nx + (Cz-1)*brick.nx*brick.ny, 5,
  4, (Cx  ) + (0)*brick.nx + (Cz  )*brick.nx*brick.ny, 4,
}
