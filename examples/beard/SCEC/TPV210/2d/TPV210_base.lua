-- refinement parameters
min_level = 0
max_level = 1
output_prefix = "TPV210_base"
data_directory = "data"
elem_order = 4

-- connectivity info
connectivity = "brick"
brick =
{
  nx = 20,
  ny = 20,
  nz =  0,
  bc0 = 1,
  bc1 = 1,
  bc2 = 1,
  bc3 = 2,
}

-- set up the domain
L = 1.5

Cx =  brick.nx/2
Cy =  brick.ny
Cz =  brick.nz/2

fault_angle = math.pi/6
cos = math.cos
sin = math.sin
tan = math.tan
function connectivity_vertices(x, y, z)
  y    = L*(y-Cy)*sin(fault_angle)
  zout = L*(z-Cz)

  xout = L*(x-Cx) - y / tan(fault_angle)
  yout = y

  return xout,yout,zout
end

function refinement_function(
  x0,y0,z0,x1,y1,z1,
  x2,y2,z2,x3,y3,z3,
  level, treeid)

  xa = (x0+x1+x2+x3)/4
  ya = (y0+y1+y2+y3)/4

  yf = -15*sin(fault_angle)
  xf = - ya / tan(fault_angle)

  x = math.abs(xa-xf)
  y = math.max(0,yf-ya)

  r = math.sqrt((x^2+y^2)/2)

  -- l = (max_level+1) + (min_level-(max_level+1))*math.sqrt(r)


  if level < max_level and r < L then
    return 1
  else
    return 0
  end
end

function element_order(
  x0,y0,z0,x1,y1,z1,
  x2,y2,z2,x3,y3,z3,
  level, treeid)

  return elem_order
end

-- material properties
cs = 3.3
cp = 5.716
rho = 2.700
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
  dt      = 0.5*dt
  nfoutput = math.ceil(tfout / dt)
  dt       = tfout / nfoutput

  noutput    = tout  / dt
  ndisp      = tdisp / dt
  nsteps     = tend / dt
  nstations  = tstations / dt

  return dt,nsteps, ndisp, noutput, nfoutput, nstations
end


volume_stations = {
  "p12_p3", 12.0, 3.0,
  "m12_p3",-12.0, 3.0,
}

fault_stations = {
  "m4.5_0",  -4.5, 0.0, 0.0, 1.0, L,
  "m7.5_0",  -7.5, 0.0, 0.0, 1.0, L,
  "m12_0" , -12.0, 0.0, 0.0, 1.0, L,
  "0_0"   ,   0.0, 0.0, 0.0, 1.0, L,
  "p4.5_0",   4.5, 0.0, 0.0, 1.0, L,
  "p7.5_0",   7.5, 0.0, 0.0, 1.0, L,
  "p12_0" ,  12.0, 0.0, 0.0, 1.0, L,
}

-- faults
n1 = sin(fault_angle)
n2 = cos(fault_angle)
m1 = n2
m2 = -n1
Snn = -7.378
Snm = 0.55*Snn
Smm = 0

function S11_fault_1(x,y,z,t)
  r = math.sqrt(x^2+y^2)
  S11_0 = (n1*Snn*n1 + n1*Snm*m1 + m1*Snm*n1 + m1*Smm*m1)*r
  return S11_0
end

function S12_fault_1(x,y,z,t)
  r = math.sqrt(x^2+y^2)
  S12_0 = (n1*Snn*n2 + n1*Snm*m2 + m1*Snm*n2 + m1*Smm*m2)*r
  return S12_0
end

function S22_fault_1(x,y,z,t)
  r = math.sqrt(x^2+y^2)
  S22_0 = (n2*Snn*n2 + n2*Snm*m2 + m2*Snm*n2 + m2*Smm*m2)*r
  return S22_0
end

fault_1 = {
  type   = "friction",
  tag    = "slip weakening",
  fs     =    0.760,
  fd     =    0.448,
  Dc     =    0.5,
  S11_0  = "S11_fault_1",
  S12_0  = "S12_fault_1",
  S22_0  = "S22_fault_1",
}

Snm_2 = (0.8)*Snn

function S11_fault_2(x,y,z,t)
  r = math.sqrt(x^2+y^2)
  S11_0 = (n1*Snn*n1 + n1*Snm_2*m1 + m1*Snm_2*n1 + m1*Smm*m1)*r
  return S11_0
end

function S12_fault_2(x,y,z,t)
  r = math.sqrt(x^2+y^2)
  S12_0 = (n1*Snn*n2 + n1*Snm_2*m2 + m1*Snm_2*n2 + m1*Smm*m2)*r
  return S12_0
end

function S22_fault_2(x,y,z,t)
  r = math.sqrt(x^2+y^2)
  S22_0 = (n2*Snn*n2 + n2*Snm_2*m2 + m2*Snm_2*n2 + m2*Smm*m2)*r
  return S22_0
end

fault_2  = {
  type   = "friction",
  tag    = "slip weakening",
  fs     =    0.760,
  fd     =    0.448,
  Dc     =    0.5,
  S11_0  = "S11_fault_2",
  S12_0  = "S12_fault_2",
  S22_0  = "S22_fault_2",
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
  fault_1,
  fault_2,
}

glueid_treeid_faceid = {
  4, (Cx-1) + (Cy-1-0)*brick.nx, 1,
  4, (Cx-1) + (Cy-1-1)*brick.nx, 1,
  4, (Cx-1) + (Cy-1-2)*brick.nx, 1,
  4, (Cx-1) + (Cy-1-3)*brick.nx, 1,
  4, (Cx-1) + (Cy-1-4)*brick.nx, 1,
  4, (Cx-1) + (Cy-1-5)*brick.nx, 1,
  4, (Cx-1) + (Cy-1-6)*brick.nx, 1,
  5, (Cx-1) + (Cy-1-7)*brick.nx, 1,
  5, (Cx-1) + (Cy-1-8)*brick.nx, 1,
  4, (Cx-1) + (Cy-1-9)*brick.nx, 1,
  --
  4, (Cx+0) + (Cy-1-0)*brick.nx, 0,
  4, (Cx+0) + (Cy-1-1)*brick.nx, 0,
  4, (Cx+0) + (Cy-1-2)*brick.nx, 0,
  4, (Cx+0) + (Cy-1-3)*brick.nx, 0,
  4, (Cx+0) + (Cy-1-4)*brick.nx, 0,
  4, (Cx+0) + (Cy-1-5)*brick.nx, 0,
  4, (Cx+0) + (Cy-1-6)*brick.nx, 0,
  5, (Cx+0) + (Cy-1-7)*brick.nx, 0,
  5, (Cx+0) + (Cy-1-8)*brick.nx, 0,
  4, (Cx+0) + (Cy-1-9)*brick.nx, 0,
}
