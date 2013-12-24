-- refinement parameters
min_level = 0
max_level = 3
output_prefix = "TPV205"
data_directory = "data"

-- connectivity info
connectivity = "brick"
brick =
{
  nx = 96,
  ny = 96,
  nz = 48,
  periodic_x = 0,
  periodic_y = 0,
  periodic_z = 0,
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
tdisp  = 0.01
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

  return dt,1, ndisp, noutput, nfoutput, nstations
end

volume_stations = {
  "p12_p3", 12.0, 3.0,
  "m12_p3",-12.0, 3.0,
}

fault_stations = {
  "m4.5_0",  -4.5, 0.0,
  "m7.5_0",  -7.5, 0.0,
  "m12_0" , -12.0, 0.0,
  "0_0"   ,   0.0, 0.0,
  "p4.5_0",   4.5, 0.0,
  "p7.5_0",   7.5, 0.0,
  "p12_0" ,  12.0, 0.0,
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

glue_info = {
  fault_1,
  fault_2,
  fault_3,
  fault_4,
}

glueid_treeid_faceid = {
1, 3604, 2,
1, 3605, 2,
1, 3648, 2,
1, 3649, 2,
4, 3652, 2,
4, 3653, 2,
1, 3664, 2,
1, 3665, 2,
1, 3668, 2,
3, 3669, 2,
3, 3840, 2,
1, 3841, 2,
1, 3844, 2,
1, 3845, 2,
2, 3856, 2,
2, 3857, 2,
1, 3860, 2,
1, 3861, 2,
1, 3904, 2,
1, 3905, 2,
1, 3262, 3,
1, 3263, 3,
1, 3306, 3,
1, 3307, 3,
4, 3310, 3,
4, 3311, 3,
1, 3322, 3,
1, 3323, 3,
1, 3326, 3,
3, 3327, 3,
3, 3498, 3,
1, 3499, 3,
1, 3502, 3,
1, 3503, 3,
2, 3514, 3,
2, 3515, 3,
1, 3518, 3,
1, 3519, 3,
1, 3562, 3,
1, 3563, 3,
}
