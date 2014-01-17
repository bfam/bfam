-- refinement parameters
min_level = 0
max_level = 0
output_prefix = "TPV14"
data_directory = "data"
elem_order = 2

-- connectivity info
connectivity = "brick"
brick =
{
  nx = 48,
  ny = 48,
  nz = 48,
  periodic_x = 0,
  periodic_y = 0,
  periodic_z = 0,
}

-- set up the domain
Lx = 3.0
Ly = 3.0
Lz = 3.0

Cx = brick.nx/2
Cy = brick.ny/2
Cz = brick.nz

function connectivity_vertices(x, y, z)
  x = Lx*(x-Cx)
  y = Ly*(y-Cy)
  z = Lz*(z-Cz)

  if x < 0 then
    x = x - 0.25
  end
  if x < -3.25 then
    x = x - 0.25
  end
  if x < -9.5 then
    x = x - 0.25
  end
  if x < -12.75 then
    x = x - 0.25
  end
  if x < -16 then
    x = x - 0.25
  end

  if y < 0 then
    q = math.pi/6
    H = math.max(y,-4*Ly)
    x = x-H*math.cos(q)
    y = y-H*math.sin(q)
  end

  return x,y,z
end

function refinement_function(
  x0,y0,z0,x1,y1,z1,
  x2,y2,z2,x3,y3,z3,
  level, treeid)

  return 0
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

  -- return dt,nsteps, ndisp, noutput, nfoutput, nstations
  return 1,1, 1, 1, 1, 1
end

-- volume_stations = {
--   "p12_p3", 12.0, 3.0,
--   "m12_p3",-12.0, 3.0,
-- }

-- fault_stations = {
--   "m4.5_0",  -4.5, 0.0,
--   "m7.5_0",  -7.5, 0.0,
--   "m12_0" , -12.0, 0.0,
--   "0_0"   ,   0.0, 0.0,
--   "p4.5_0",   4.5, 0.0,
--   "p7.5_0",   7.5, 0.0,
--   "p12_0" ,  12.0, 0.0,
-- }

-- faults
main_fault = {
  type   = "friction",
  tag    = "slip weakening",
  fs     =    0.677,
  fd     =    0.525,
  Dc     =    0.4,
  S12_0  =   70.0,
  S22_0  = -120.0,
}

nucleation_patch  = {
  type   = "friction",
  tag    = "slip weakening",
  fs     =    0.677,
  fd     =    0.525,
  Dc     =    0.4,
  S12_0  =   81.6,
  S22_0  = -120.0,
}

branch_fault  = {
  type   = "friction",
  tag    = "slip weakening",
  fs     =    0.677,
  fd     =    0.525,
  Dc     =    0.4,
  S12_0  =   70.0,
  S22_0  = -120.0,
}

glue_info = {
  main_fault,
  nucleation_patch,
  branch_fault,
}

glueid_treeid_faceid = {
  1, 901, 2,
  1, 912, 2,
  1, 916, 2,
  1, 917, 2,
  1, 960, 2,
  1, 961, 2,
  1, 964, 2,
  1, 965, 2,
  1, 815, 3,
  1, 826, 3,
  1, 830, 3,
  1, 831, 3,
  1, 874, 3,
  1, 875, 3,
  1, 878, 3,
  1, 879, 3,
  2, 913, 2,
  2, 827, 3,
  3, 831, 1,
  3, 829, 1,
  3, 823, 1,
  3, 821, 1,
  3, 874, 0,
  3, 872, 0,
  3, 866, 0,
  3, 864, 0,
}
