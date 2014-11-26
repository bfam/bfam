-- refinement parameters
min_level = 0
max_level = 1
output_prefix = "TPV205"
data_directory = "data"
elem_order = 4

-- connectivity info
connectivity = "brick"
brick =
{
  nx = 2*(10+20),
  ny = 2*(10+20),
  nz = 10+10,
  bc0 = 1,
  bc1 = 1,
  bc2 = 1,
  bc3 = 1,
  bc4 = 2,
  bc5 = 1,
}

-- set up the domain
Lx = 1.5
Ly = 1.5
Lz = 1.5

Cx = brick.nx/2
Cy = brick.ny/2
Cz = brick.nz

function connectivity_vertices(x, y, z)
  xout = Lx*(x-Cx)
  yout = Ly*(y-Cy)
  zout = Lz*(z-Cz)
  return xout,yout,zout
end

function refinement_function(
  x0,y0,z0,x1,y1,z1,
  x2,y2,z2,x3,y3,z3,
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

  x = math.max(xmin-15,0) / (Cx*Lx-15)
  y = ymin                / (Cy*Ly)
  r = math.sqrt((x^2+y^2)/2)
  l = (max_level+1) + (min_level-(max_level+1))*math.sqrt(r)


  if level < l then
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
  "m4.5_0",  -4.5, 0.0, 0.0, 1.0, Ly,
  "m7.5_0",  -7.5, 0.0, 0.0, 1.0, Ly,
  "m12_0" , -12.0, 0.0, 0.0, 1.0, Ly,
  "0_0"   ,   0.0, 0.0, 0.0, 1.0, Ly,
  "p4.5_0",   4.5, 0.0, 0.0, 1.0, Ly,
  "p7.5_0",   7.5, 0.0, 0.0, 1.0, Ly,
  "p12_0" ,  12.0, 0.0, 0.0, 1.0, Ly,
}

-- faults
fault_1 = {
  type   = "friction",
  tag    = "slip weakening",
  fs     =    0.677,
  fd     =    0.525,
  Dc     =    0.4,
  S12_0  =   70.0,
  S22_0  = -120.0,
}

fault_2  = {
  type   = "friction",
  tag    = "slip weakening",
  fs     =    0.677,
  fd     =    0.525,
  Dc     =    0.4,
  S12_0  =   78.0,
  S22_0  = -120.0,
}

fault_3  = {
  type   = "friction",
  tag    = "slip weakening",
  fs     =    0.677,
  fd     =    0.525,
  Dc     =    0.4,
  S12_0  =   81.6,
  S22_0  = -120.0,
}

fault_4 = {
  type  = "friction",
  tag   = "slip weakening",
  fs    =    0.677,
  fd    =    0.525,
  Dc    =    0.4,
  S12_0 =   62.0,
  S22_0 = -120.0,
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
  fault_3,
  fault_4,
}

glueid_treeid_faceid = {
  4, (Cx-10) + (Cy-1)*brick.nx, 3,
  4, (Cx- 9) + (Cy-1)*brick.nx, 3,
  4, (Cx- 8) + (Cy-1)*brick.nx, 3,
  4, (Cx- 7) + (Cy-1)*brick.nx, 3,
  5, (Cx- 6) + (Cy-1)*brick.nx, 3,
  5, (Cx- 5) + (Cy-1)*brick.nx, 3,
  4, (Cx- 4) + (Cy-1)*brick.nx, 3,
  4, (Cx- 3) + (Cy-1)*brick.nx, 3,
  4, (Cx- 2) + (Cy-1)*brick.nx, 3,
  6, (Cx- 1) + (Cy-1)*brick.nx, 3,
  6, (Cx+ 0) + (Cy-1)*brick.nx, 3,
  4, (Cx+ 1) + (Cy-1)*brick.nx, 3,
  4, (Cx+ 2) + (Cy-1)*brick.nx, 3,
  4, (Cx+ 3) + (Cy-1)*brick.nx, 3,
  7, (Cx+ 4) + (Cy-1)*brick.nx, 3,
  7, (Cx+ 5) + (Cy-1)*brick.nx, 3,
  4, (Cx+ 6) + (Cy-1)*brick.nx, 3,
  4, (Cx+ 7) + (Cy-1)*brick.nx, 3,
  4, (Cx+ 8) + (Cy-1)*brick.nx, 3,
  4, (Cx+ 9) + (Cy-1)*brick.nx, 3,
  --
  4, (Cx-10) + (Cy+0)*brick.nx, 2,
  4, (Cx- 9) + (Cy+0)*brick.nx, 2,
  4, (Cx- 8) + (Cy+0)*brick.nx, 2,
  4, (Cx- 7) + (Cy+0)*brick.nx, 2,
  5, (Cx- 6) + (Cy+0)*brick.nx, 2,
  5, (Cx- 5) + (Cy+0)*brick.nx, 2,
  4, (Cx- 4) + (Cy+0)*brick.nx, 2,
  4, (Cx- 3) + (Cy+0)*brick.nx, 2,
  4, (Cx- 2) + (Cy+0)*brick.nx, 2,
  6, (Cx- 1) + (Cy+0)*brick.nx, 2,
  6, (Cx+ 0) + (Cy+0)*brick.nx, 2,
  4, (Cx+ 1) + (Cy+0)*brick.nx, 2,
  4, (Cx+ 2) + (Cy+0)*brick.nx, 2,
  4, (Cx+ 3) + (Cy+0)*brick.nx, 2,
  7, (Cx+ 4) + (Cy+0)*brick.nx, 2,
  7, (Cx+ 5) + (Cy+0)*brick.nx, 2,
  4, (Cx+ 6) + (Cy+0)*brick.nx, 2,
  4, (Cx+ 7) + (Cy+0)*brick.nx, 2,
  4, (Cx+ 8) + (Cy+0)*brick.nx, 2,
  4, (Cx+ 9) + (Cy+0)*brick.nx, 2,
}
