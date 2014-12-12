cos = math.cos
sin = math.sin
output_prefix = "TPV15_base"
data_directory = "data"
elem_order = 4
max_level = 600

-- refinement parameters
h1_targ = elem_order*0.1
r1_targ = 5
h2_targ = h1_targ*2

r2_targ = r1_targ*2
h2_targ = h1_targ*2

r3_targ = r2_targ*2
h3_targ = h2_targ*2

nuc_width = h1_targ/elem_order

mesh_file = "/home/jekozdon/codes/bfam/examples/beard/SCEC/TPV15/3d/TPV15.inp"
glueid_treeid_faceid = mesh_file

default_boundary_tag = "non-reflecting"

-- ABQ->p4est:
-- 1 -> 4
-- 2 -> 5
-- 3 -> 2
-- 4 -> 1
-- 5 -> 3
-- 6 -> 0



q_branch = math.pi/6
function fault_distance(x,y,z)
  -- distance to the main fault
  xf = math.max(0,math.abs(x+2)-14)
  yf = math.abs(y)
  zf = math.max(0,math.abs(z)-15)
  r  = xf^2 + yf^2 + zf^2

  if y < 0 then
    if x*cos(q_branch) - y*sin(q_branch) < 12 then
      -- distance to the branch fault
      rb = (sin(q_branch)*x + cos(q_branch)*y)^2 + zf^2
      r  = math.min(rb,r)
    else
      -- distance to the end of the fault
      xe =  12*cos(q_branch)
      ye = -12*sin(q_branch)
      re = (x-xe)^2 + (y-ye)^2 + zf^2
      r  = math.min(re,r)
    end
  end

  return math.sqrt(r)
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
  hmax = math.max(hmax,h)
  hmin = math.min(hmin,h)

  h = (x0-x4)^2+(y0-y4)^2+(z0-z4)^2
  hmax = math.max(hmax,h)
  hmin = math.min(hmin,h)

  h = (x1-x3)^2+(y1-y3)^2+(z1-z3)^2
  hmax = math.max(hmax,h)
  hmin = math.min(hmin,h)

  h = (x1-x5)^2+(y1-y5)^2+(z1-z5)^2
  hmax = math.max(hmax,h)
  hmin = math.min(hmin,h)

  h = (x2-x3)^2+(y2-y3)^2+(z2-z3)^2
  hmax = math.max(hmax,h)
  hmin = math.min(hmin,h)

  h = (x2-x6)^2+(y2-y6)^2+(z2-z6)^2
  hmax = math.max(hmax,h)
  hmin = math.min(hmin,h)

  h = (x3-x7)^2+(y3-y7)^2+(z3-z7)^2
  hmax = math.max(hmax,h)
  hmin = math.min(hmin,h)

  h = (x4-x5)^2+(y4-y5)^2+(z4-z5)^2
  hmax = math.max(hmax,h)
  hmin = math.min(hmin,h)

  h = (x4-x6)^2+(y4-y6)^2+(z4-z6)^2
  hmax = math.max(hmax,h)
  hmin = math.min(hmin,h)

  h = (x5-x7)^2+(y5-y7)^2+(z5-z7)^2
  hmax = math.max(hmax,h)
  hmin = math.min(hmin,h)

  h = (x6-x7)^2+(y6-y7)^2+(z6-z7)^2
  hmax = math.max(hmax,h)
  hmin = math.min(hmin,h)

  return math.sqrt(hmin), math.sqrt(hmax)

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

  -- x = (x0+x1+x2+x3+x4+x5+x6+x7)/8
  -- y = (y0+y1+y2+y3+y4+y5+y6+y7)/8
  -- z = (z0+z1+z2+z3+z4+z5+z6+z7)/8
  -- r = fault_distance(x,y,z)

  r = fault_distance(x0,y0,z0)
  r = math.min(r,fault_distance(x1,y1,z1))
  r = math.min(r,fault_distance(x2,y2,z2))
  r = math.min(r,fault_distance(x3,y3,z3))
  r = math.min(r,fault_distance(x4,y4,z4))
  r = math.min(r,fault_distance(x5,y5,z5))
  r = math.min(r,fault_distance(x6,y6,z6))
  r = math.min(r,fault_distance(x7,y7,z7))

  hmin, hmax = element_size( x0,y0,z0,x1,y1,z1,
                             x2,y2,z2,x3,y3,z3,
                             x4,y4,z4,x5,y5,z5,
                             x6,y6,z6,x7,y7,z7)

  if hmax > math.max(2,r)/2*h1_targ then
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

  return dt,nsteps, ndisp, noutput, nfoutput, nstations
end

-- faults
function S12_main_fault(x,y,z,t)
  S12_0  = -70.0
  r = math.max(math.abs(x+8),math.abs(7.5+z))-1.5;
  v = 1 - 1/(1+math.exp(-2*r/nuc_width))

  return -70 + v*(-81.6+70)
end
main_fault = {
  type   = "friction",
  tag    = "slip weakening",
  fs     =    0.677,
  fd     =    0.525,
  Dc     =    0.4,
  S12_0  =   "S12_main_fault",
  S22_0  = -120.0,
}

snn    =  -120
snm    =   -78
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

rx =  cos(q_branch)
ry = -sin(q_branch)
fault_stations = {
  "faultst-020dp000", -2.0   , 0.0   ,  0.0,  0.0, 1.0, 0.0,  0.1,
  "faultst020dp000" ,  2.0   , 0.0   ,  0.0,  0.0, 1.0, 0.0,  0.1,
  "faultst050dp000" ,  5.0   , 0.0   ,  0.0,  0.0, 1.0, 0.0,  0.1,
  "faultst090dp000" ,  9.0   , 0.0   ,  0.0,  0.0, 1.0, 0.0,  0.1,
  "branchst020dp000",  2.0*rx, 2.0*ry,  0.0,  -ry,  rx, 0.0,  0.1,
  "branchst050dp000",  5.0*rx, 5.0*ry,  0.0,  -ry,  rx, 0.0,  0.1,
  "branchst090dp000",  9.0*rx, 9.0*ry,  0.0,  -ry,  rx, 0.0,  0.1,
  "faultst-020dp075", -2.0   , 0.0   , -7.5,  0.0, 1.0, 0.0,  0.1,
  "faultst020dp075" ,  2.0   , 0.0   , -7.5,  0.0, 1.0, 0.0,  0.1,
  "faultst050dp075" ,  5.0   , 0.0   , -7.5,  0.0, 1.0, 0.0,  0.1,
  "faultst090dp075" ,  9.0   , 0.0   , -7.5,  0.0, 1.0, 0.0,  0.1,
  "branchst020dp075",  2.0*rx, 2.0*ry, -7.5,  -ry,  rx, 0.0,  0.1,
  "branchst050dp075",  5.0*rx, 5.0*ry, -7.5,  -ry,  rx, 0.0,  0.1,
  "branchst090dp075",  9.0*rx, 9.0*ry, -7.5,  -ry,  rx, 0.0,  0.1,
}
