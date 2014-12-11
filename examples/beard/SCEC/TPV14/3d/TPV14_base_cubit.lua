cos = math.cos
sin = math.sin
-- refinement parameters
min_level = 0
max_level = 0
output_prefix = "TPV14_base"
data_directory = "data"
elem_order = 1

mesh_file = "/home/jekozdon/codes/bfam/examples/beard/SCEC/TPV14/3d/TPV14.inp"

-- connectivity info
connectivity = "brick"
buffer = 40
brick =
{
  nx = 32+24+2*buffer,
  ny = 24+   2*buffer,
  nz = 30+     buffer,
  periodic_x = 0,
  periodic_y = 0,
  periodic_z = 0,
  bc0 = 1,
  bc1 = 1,
  bc2 = 1,
  bc3 = 1,
  bc4 = 1,
  bc5 = 2,
}

-- set up the domain
L = 0.5

Cx = (brick.nx-(32+24))/2+32
Cy = (brick.ny-( 0+24))/2 +24
Cz = brick.nz

q_branch = math.pi/6

function refinement_function(
  x0,y0,z0,x1,y1,z1,
  x2,y2,z2,x3,y3,z3,
  x4,y4,z4,x5,y5,z5,
  x6,y6,z6,x7,y7,z7,
  level, treeid)

  x = (x0+x1+x2+x3+x4+x5+x6+x7)/8
  y = (y0+y1+y2+y3+y4+y5+y6+y7)/8
  z = (z0+z1+z2+z3+z4+z5+z6+z7)/8

  -- distance to the main fault
  xf = math.max(0,math.abs(x+2)-14)
  yf = math.abs(y)
  zf = math.max(0,math.abs(z)-15)
  r  = math.sqrt(xf^2+yf^2+zf^2)

  if y < 0 then
    if x*cos(q_branch) - y*sin(q_branch) < 12 then
      -- distance to the branch fault
      rb = math.abs(sin(q_branch)*x + cos(q_branch)*y + zf^2)
      r  = math.min(rb,r)
    else
      -- distance to the end of the fault
      xe =  12*cos(q_branch)
      ye = -12*sin(q_branch)
      re = math.sqrt((x-xe)^2 + (y-ye)^2 + zf^2)
      r  = math.min(re,r)
    end
  end

  r = r/10


  l = max_level + (min_level-max_level)*math.sqrt(r)

  if level < l then
    return 1
  else
    return 0
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
  nsteps = 1
  nstations  = tstations / dt

  return dt,nsteps, ndisp, noutput, nfoutput, nstations
end

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

snn    =  -120
snm    =    70
smm    =     0
n1     = sin(q_branch)
n2     = cos(q_branch)
m1     = n2
m2     =-n1

branch_fault  = {
  type   = "friction",
  tag    = "slip weakening",
  fs     =    0.677,
  fd     =    0.525,
  Dc     =    0.4,
  S22_0  = n2*snn*n2 + n2*snm*m2 + m2*snm*n2 + m2*smm*m2,
  S12_0  = n1*snn*n2 + n1*snm*m2 + m1*snm*n2 + m1*smm*m2,
  S11_0  = n1*snn*n1 + n1*snm*m1 + m1*snm*n1 + m1*smm*m1,
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
  main_fault,
  branch_fault,
}

-- ABQ->p4est:
-- 1 -> 4
-- 2 -> 5
-- 3 -> 2
-- 4 -> 1
-- 5 -> 3
-- 6 -> 0

glueid_treeid_faceid = {
      -4,   6, 2,
      -4,  14, 2,
      -4,  22, 2,
      -4,   4, 3,
      -4,  12, 3,
      -4,  20, 3,
      -4,  24, 3,
      -4,  26, 3,
      -4,  30, 3,
      -4,  32, 3,
      -4,  36, 3,
      -4,  38, 3,
      -4,   7, 0,
      -4,  15, 0,
      -4,  23, 0,
      -5,  24, 1,
      -5,  25, 1,
      -5,  30, 1,
      -5,  31, 1,
      -5,  36, 1,
      -5,  37, 1,
}
