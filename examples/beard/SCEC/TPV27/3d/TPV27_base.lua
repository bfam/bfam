max = math.max
min = math.min
abs = math.abs
sqrt = math.sqrt
-- refinement parameters
max_level = 300
output_prefix = "TPV27_base"
data_directory = "data"
elem_order = 4
h1_targ   = elem_order*0.1
r_targ    = 5
D_targ    = 20
hmax_targ = elem_order*2
-- hmax_targ = h1_targ

-- connectivity info
connectivity = "brick"
BUF = 1
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
Lx = 20
Ly = 20
Lz = 20

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
  xf = max(0,abs(x)-20)
  zf = abs(z)
  yf = max(0,abs(y)-20)
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


  return elem_order, "plastic"
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

-- plasticity parameters
plastic = {
  tag = "Duvaut-Lions",
  c  = 1.36,   -- plastic cohesion
  Tr = 0.03,   -- viscoplastic relaxation time
  nu = 0.1934, -- bulk friction
  phi = math.atan(0.1934),  -- angle of friction: atan(nu)
  pf     = "pf_0_function",
  S11_0  = "S11_0_function",
  S12_0  = 0,
  S13_0  = "S13_0_function",
  S22_0  = "S22_0_function",
  S23_0  = 0,
  S33_0  = "S33_0_function",
}

-- time stepper to use
lsrk_method  = "KC54"

tend       = 13
-- tout       = 1
-- tfout      = 0.01
tout       = -1
tfout      = -1
tdisp      = 0.01
tstations  = 0.01
nerr       = 0
output_file_fault = 1

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
  -- nsteps = 1

  return dt,nsteps, ndisp, noutput, nfoutput, nstations
end

fault_stations = {
  "faultst-050dp000", -5.0,  0.0, 0.0, 0.0, 0.0, 1.0, 0.1,
  "faultst150dp000" , 15.0,  0.0, 0.0, 0.0, 0.0, 1.0, 0.1,
  "faultst-050dp050", -5.0,  5.0, 0.0, 0.0, 0.0, 1.0, 0.1,
  "faultst150dp050" , 15.0,  5.0, 0.0, 0.0, 0.0, 1.0, 0.1,
  "faultst-150dp100",-15.0, 10.0, 0.0, 0.0, 0.0, 1.0, 0.1,
  "faultst-050dp100", -5.0, 10.0, 0.0, 0.0, 0.0, 1.0, 0.1,
  "faultst000dp100" ,  0.0, 10.0, 0.0, 0.0, 0.0, 1.0, 0.1,
  "faultst050dp100" ,  5.0, 10.0, 0.0, 0.0, 0.0, 1.0, 0.1,
  "faultst100dp100" , 10.0, 10.0, 0.0, 0.0, 0.0, 1.0, 0.1,
  "faultst150dp100" , 15.0, 10.0, 0.0, 0.0, 0.0, 1.0, 0.1,
  "faultst-050dp150", -5.0, 15.0, 0.0, 0.0, 0.0, 1.0, 0.1,
  "faultst150dp150" , 15.0, 15.0, 0.0, 0.0, 0.0, 1.0, 0.1,
}

-- faults
b11 =  0.926793
b33 =  1.073206
b13 = -0.169029
function Omega_function(x,y,z,t)
  if y < 15 then
    return 1
  elseif y < 20 then
    return (20.0-y)/5.0
  else
    return 0
  end
  -- return min(1.0, max(0.0, (20.0-y)/5.0))
end
function pf_0_function(x,y,z,t)
  return 9.8*y -- in MPa
end
function S22_0_function(x,y,z,t)
  return -2.670*9.8*y -- in MPa
end
function S11_0_function(x,y,z,t)
  W     = Omega_function(x,y,z,t)
  pf    = pf_0_function(x,y,z,t)
  S22   = S22_0_function(x,y,z,t)
  return W*(b11*(S22+pf)-pf)+(1-W)*S22
end
function S33_0_function(x,y,z,t)
  W     = Omega_function(x,y,z,t)
  pf    = pf_0_function(x,y,z,t)
  S22   = S22_0_function(x,y,z,t)
  return W*(b33*(S22 + pf) - pf) + (1 - W)*S22
end
function S13_0_function(x,y,z,t)
  W     = Omega_function(x,y,z,t)
  pf    = pf_0_function(x,y,z,t)
  S22   = S22_0_function(x,y,z,t)
  return W*b13*(S22+pf)
end
function c0_function(x,y,z,t)
  return 0.4 + 0.72*max(0.0,5-y)
end
function Tforce_function(x,y,z,t)
  rcrit = 4
  r = sqrt((x+5)^2+(y-10)^2)
  -- r = sqrt(x^2+(y-Ly/2)^2)
  Vr = 0.7*cs
  if r < rcrit then
    return r/Vr + (0.081*rcrit/Vr)*(1/(1-(r/rcrit)^2) - 1)
  else
    return 1e9
  end
end

fault = {
  type   = "friction",
  tag    = "slip weakening",
  fs     = 0.18,
  fd     = 0.12,
  Dc     = 0.3,
  pf_0   = "pf_0_function",
  S11_0  = "S11_0_function",
  S12_0  = 0,
  S13_0  = "S13_0_function",
  S22_0  = "S22_0_function",
  S23_0  = 0,
  S33_0  = "S33_0_function",
  c0     = "c0_function",
  Tforce = "Tforce_function",
  Tforce_0 = 0.5,
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
