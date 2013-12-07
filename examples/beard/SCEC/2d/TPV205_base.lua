-- refinement parameters
min_level = 3
max_level = 5
output_prefix = "TPV205"
data_directory = "data"

-- connectivity info
connectivity = "brick"
brick =
{
  nx = 40,
  ny = 20,
  nz = 20,
  periodic_x = 1,
  periodic_y = 1,
  periodic_z = 1,
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

  if level < min_level then
    return 1
  elseif level >= max_level then
    return 0
  elseif xmin <= 15 and ymin <= 0.1 then
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
  N = 1

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
tfout  = 0.01
tdisp  = 0.01
nerr   = 0

function time_step_parameters(dt)
  dt      = 0.5*dt
  nfoutput = math.ceil(tfout / dt)
  dt       = tfout / nfoutput

  noutput  = tout  / dt
  ndisp   = tdisp / dt
  nsteps  = tend / dt

  return dt,nsteps, ndisp, noutput, nfoutput
end

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
1, 204, 2,
1, 205, 2,
1, 216, 2,
1, 217, 2,
4, 220, 2,
4, 221, 2,
1, 392, 2,
1, 393, 2,
1, 396, 2,
3, 397, 2,
3, 408, 2,
1, 409, 2,
1, 412, 2,
1, 413, 2,
2, 456, 2,
2, 457, 2,
1, 460, 2,
1, 461, 2,
1, 472, 2,
1, 473, 2,
1, 198, 3,
1, 199, 3,
1, 210, 3,
1, 211, 3,
4, 214, 3,
4, 215, 3,
1, 386, 3,
1, 387, 3,
1, 390, 3,
3, 391, 3,
3, 402, 3,
1, 403, 3,
1, 406, 3,
1, 407, 3,
2, 450, 3,
2, 451, 3,
1, 454, 3,
1, 455, 3,
1, 466, 3,
1, 467, 3,
}
