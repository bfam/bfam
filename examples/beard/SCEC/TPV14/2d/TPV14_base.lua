cos = math.cos
sin = math.sin
-- refinement parameters
min_level = 0
max_level = 3
output_prefix = "TPV14_base"
data_directory = "data"
elem_order = 2

-- connectivity info
connectivity = "brick"
buffer = 80
brick =
{
  nx = 32+24+buffer,
  ny = 24+buffer,
  nz = 1,
  periodic_x = 0,
  periodic_y = 0,
  periodic_z = 0,
  bc0 = 1,
  bc1 = 1,
  bc2 = 1,
  bc3 = 1,
}

-- set up the domain
L = 0.5

Cx = (brick.nx-(32+24))/2+32
Cy = (brick.ny-( 0+24))/2 +24
Cz = brick.nz

q_branch = math.pi/6

function connectivity_vertices(x, y, z)
  x = L*(x-Cx)
  y = L*(y-Cy)
  z = L*(z-Cz)

  xout = x
  yout = y
  if y < 0 then
    H = math.max(y,-12)
    alpha = 1-math.max(0,-16-x)/(Cx*L-16)
    xout = xout-H*alpha*cos(q_branch)
    yout = yout-H*sin(q_branch)
  end
  if y > -12 then
    H = math.min(12,12+y)
    alpha = math.max(0,(x-12)/((brick.nx-Cx)*L-12))
    xout = xout + H*alpha*cos(q_branch)
  end

  return xout,yout,z
end

function refinement_function(
  x0,y0,z0,x1,y1,z1,
  x2,y2,z2,x3,y3,z3,
  level, treeid)

  x = (x0+x1+x2+x3)/4
  y = (y0+y1+y2+y3)/4

  -- distance to the main fault
  xf = math.max(0,math.abs(x+2)-14)
  yf = math.abs(y)
  r  = math.sqrt(xf^2+yf^2)

  if y < 0 then
    if x*cos(q_branch) - y*sin(q_branch) < 12 then
      -- distance to the branch fault
      rb = math.abs(sin(q_branch)*x + cos(q_branch)*y)
      r  = math.min(rb,r)
    else
      -- distance to the end of the fault
      xe =  12*cos(q_branch)
      ye = -12*sin(q_branch)
      re = math.sqrt((x-xe)^2 + (y-ye)^2)
      r  = math.min(re,r)
    end
  end

  r = r/4


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
tout   = 0.1
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

rx =  cos(q_branch)
ry = -sin(q_branch)
fault_stations = {
  "faultst-020dp075", -2.0, 0.0, 0.0, 1.0, 0.1,
  "faultst020dp075" ,  2.0, 0.0, 0.0, 1.0, 0.1,
  "faultst050dp075" ,  5.0, 0.0, 0.0, 1.0, 0.1,
  "faultst090dp075" ,  9.0, 0.0, 0.0, 1.0, 0.1,
  "branchst020dp075",  2.0*rx, 2.0*ry, -ry,  rx, 0.1,
  "branchst050dp075",  5.0*rx, 5.0*ry, -ry,  rx, 0.1,
  "branchst090dp075",  9.0*rx, 9.0*ry, -ry,  rx, 0.1,
}

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
  nucleation_patch,
  branch_fault,
}

glueid_treeid_faceid = {
  4, (Cx-32) + (Cy-1)*brick.nx, 3,
  4, (Cx-31) + (Cy-1)*brick.nx, 3,
  4, (Cx-30) + (Cy-1)*brick.nx, 3,
  4, (Cx-29) + (Cy-1)*brick.nx, 3,
  4, (Cx-28) + (Cy-1)*brick.nx, 3,
  4, (Cx-27) + (Cy-1)*brick.nx, 3,
  4, (Cx-26) + (Cy-1)*brick.nx, 3,
  4, (Cx-25) + (Cy-1)*brick.nx, 3,
  4, (Cx-24) + (Cy-1)*brick.nx, 3,
  4, (Cx-23) + (Cy-1)*brick.nx, 3,
  4, (Cx-22) + (Cy-1)*brick.nx, 3,
  4, (Cx-21) + (Cy-1)*brick.nx, 3,
  4, (Cx-20) + (Cy-1)*brick.nx, 3,
  --
  5, (Cx-19) + (Cy-1)*brick.nx, 3,
  5, (Cx-18) + (Cy-1)*brick.nx, 3,
  5, (Cx-17) + (Cy-1)*brick.nx, 3,
  5, (Cx-16) + (Cy-1)*brick.nx, 3,
  5, (Cx-15) + (Cy-1)*brick.nx, 3,
  5, (Cx-14) + (Cy-1)*brick.nx, 3,
  --
  4, (Cx-13) + (Cy-1)*brick.nx, 3,
  4, (Cx-12) + (Cy-1)*brick.nx, 3,
  4, (Cx-11) + (Cy-1)*brick.nx, 3,
  4, (Cx-10) + (Cy-1)*brick.nx, 3,
  4, (Cx- 9) + (Cy-1)*brick.nx, 3,
  4, (Cx- 8) + (Cy-1)*brick.nx, 3,
  4, (Cx- 7) + (Cy-1)*brick.nx, 3,
  4, (Cx- 6) + (Cy-1)*brick.nx, 3,
  4, (Cx- 5) + (Cy-1)*brick.nx, 3,
  4, (Cx- 4) + (Cy-1)*brick.nx, 3,
  4, (Cx- 3) + (Cy-1)*brick.nx, 3,
  4, (Cx- 2) + (Cy-1)*brick.nx, 3,
  4, (Cx- 1) + (Cy-1)*brick.nx, 3,
  4, (Cx+ 0) + (Cy-1)*brick.nx, 3,
  4, (Cx+ 1) + (Cy-1)*brick.nx, 3,
  4, (Cx+ 2) + (Cy-1)*brick.nx, 3,
  4, (Cx+ 3) + (Cy-1)*brick.nx, 3,
  4, (Cx+ 4) + (Cy-1)*brick.nx, 3,
  4, (Cx+ 5) + (Cy-1)*brick.nx, 3,
  4, (Cx+ 6) + (Cy-1)*brick.nx, 3,
  4, (Cx+ 7) + (Cy-1)*brick.nx, 3,
  4, (Cx+ 8) + (Cy-1)*brick.nx, 3,
  4, (Cx+ 9) + (Cy-1)*brick.nx, 3,
  4, (Cx+10) + (Cy-1)*brick.nx, 3,
  4, (Cx+11) + (Cy-1)*brick.nx, 3,
  4, (Cx+12) + (Cy-1)*brick.nx, 3,
  4, (Cx+13) + (Cy-1)*brick.nx, 3,
  4, (Cx+14) + (Cy-1)*brick.nx, 3,
  4, (Cx+15) + (Cy-1)*brick.nx, 3,
  4, (Cx+16) + (Cy-1)*brick.nx, 3,
  4, (Cx+17) + (Cy-1)*brick.nx, 3,
  4, (Cx+18) + (Cy-1)*brick.nx, 3,
  4, (Cx+19) + (Cy-1)*brick.nx, 3,
  4, (Cx+20) + (Cy-1)*brick.nx, 3,
  4, (Cx+21) + (Cy-1)*brick.nx, 3,
  4, (Cx+22) + (Cy-1)*brick.nx, 3,
  4, (Cx+23) + (Cy-1)*brick.nx, 3,
  --
  4, (Cx-32) + (Cy+0)*brick.nx, 2,
  4, (Cx-31) + (Cy+0)*brick.nx, 2,
  4, (Cx-30) + (Cy+0)*brick.nx, 2,
  4, (Cx-29) + (Cy+0)*brick.nx, 2,
  4, (Cx-28) + (Cy+0)*brick.nx, 2,
  4, (Cx-27) + (Cy+0)*brick.nx, 2,
  4, (Cx-26) + (Cy+0)*brick.nx, 2,
  4, (Cx-25) + (Cy+0)*brick.nx, 2,
  4, (Cx-24) + (Cy+0)*brick.nx, 2,
  4, (Cx-23) + (Cy+0)*brick.nx, 2,
  4, (Cx-22) + (Cy+0)*brick.nx, 2,
  4, (Cx-21) + (Cy+0)*brick.nx, 2,
  4, (Cx-20) + (Cy+0)*brick.nx, 2,
  --
  5, (Cx-19) + (Cy+0)*brick.nx, 2,
  5, (Cx-18) + (Cy+0)*brick.nx, 2,
  5, (Cx-17) + (Cy+0)*brick.nx, 2,
  5, (Cx-16) + (Cy+0)*brick.nx, 2,
  5, (Cx-15) + (Cy+0)*brick.nx, 2,
  5, (Cx-14) + (Cy+0)*brick.nx, 2,
  --
  4, (Cx-13) + (Cy+0)*brick.nx, 2,
  4, (Cx-12) + (Cy+0)*brick.nx, 2,
  4, (Cx-11) + (Cy+0)*brick.nx, 2,
  4, (Cx-10) + (Cy+0)*brick.nx, 2,
  4, (Cx- 9) + (Cy+0)*brick.nx, 2,
  4, (Cx- 8) + (Cy+0)*brick.nx, 2,
  4, (Cx- 7) + (Cy+0)*brick.nx, 2,
  4, (Cx- 6) + (Cy+0)*brick.nx, 2,
  4, (Cx- 5) + (Cy+0)*brick.nx, 2,
  4, (Cx- 4) + (Cy+0)*brick.nx, 2,
  4, (Cx- 3) + (Cy+0)*brick.nx, 2,
  4, (Cx- 2) + (Cy+0)*brick.nx, 2,
  4, (Cx- 1) + (Cy+0)*brick.nx, 2,
  4, (Cx+ 0) + (Cy+0)*brick.nx, 2,
  4, (Cx+ 1) + (Cy+0)*brick.nx, 2,
  4, (Cx+ 2) + (Cy+0)*brick.nx, 2,
  4, (Cx+ 3) + (Cy+0)*brick.nx, 2,
  4, (Cx+ 4) + (Cy+0)*brick.nx, 2,
  4, (Cx+ 5) + (Cy+0)*brick.nx, 2,
  4, (Cx+ 6) + (Cy+0)*brick.nx, 2,
  4, (Cx+ 7) + (Cy+0)*brick.nx, 2,
  4, (Cx+ 8) + (Cy+0)*brick.nx, 2,
  4, (Cx+ 9) + (Cy+0)*brick.nx, 2,
  4, (Cx+10) + (Cy+0)*brick.nx, 2,
  4, (Cx+11) + (Cy+0)*brick.nx, 2,
  4, (Cx+12) + (Cy+0)*brick.nx, 2,
  4, (Cx+13) + (Cy+0)*brick.nx, 2,
  4, (Cx+14) + (Cy+0)*brick.nx, 2,
  4, (Cx+15) + (Cy+0)*brick.nx, 2,
  4, (Cx+16) + (Cy+0)*brick.nx, 2,
  4, (Cx+17) + (Cy+0)*brick.nx, 2,
  4, (Cx+18) + (Cy+0)*brick.nx, 2,
  4, (Cx+19) + (Cy+0)*brick.nx, 2,
  4, (Cx+20) + (Cy+0)*brick.nx, 2,
  4, (Cx+21) + (Cy+0)*brick.nx, 2,
  4, (Cx+22) + (Cy+0)*brick.nx, 2,
  4, (Cx+23) + (Cy+0)*brick.nx, 2,
  --
  6, (Cx-1) + (Cy- 1)*brick.nx, 1,
  6, (Cx-1) + (Cy- 2)*brick.nx, 1,
  6, (Cx-1) + (Cy- 3)*brick.nx, 1,
  6, (Cx-1) + (Cy- 4)*brick.nx, 1,
  6, (Cx-1) + (Cy- 5)*brick.nx, 1,
  6, (Cx-1) + (Cy- 6)*brick.nx, 1,
  6, (Cx-1) + (Cy- 7)*brick.nx, 1,
  6, (Cx-1) + (Cy- 8)*brick.nx, 1,
  6, (Cx-1) + (Cy- 9)*brick.nx, 1,
  6, (Cx-1) + (Cy-10)*brick.nx, 1,
  6, (Cx-1) + (Cy-11)*brick.nx, 1,
  6, (Cx-1) + (Cy-12)*brick.nx, 1,
  6, (Cx-1) + (Cy-13)*brick.nx, 1,
  6, (Cx-1) + (Cy-14)*brick.nx, 1,
  6, (Cx-1) + (Cy-15)*brick.nx, 1,
  6, (Cx-1) + (Cy-16)*brick.nx, 1,
  6, (Cx-1) + (Cy-17)*brick.nx, 1,
  6, (Cx-1) + (Cy-18)*brick.nx, 1,
  6, (Cx-1) + (Cy-19)*brick.nx, 1,
  6, (Cx-1) + (Cy-20)*brick.nx, 1,
  6, (Cx-1) + (Cy-21)*brick.nx, 1,
  6, (Cx-1) + (Cy-22)*brick.nx, 1,
  6, (Cx-1) + (Cy-23)*brick.nx, 1,
  6, (Cx-1) + (Cy-24)*brick.nx, 1,
  --
  6, (Cx+0) + (Cy- 1)*brick.nx, 0,
  6, (Cx+0) + (Cy- 2)*brick.nx, 0,
  6, (Cx+0) + (Cy- 3)*brick.nx, 0,
  6, (Cx+0) + (Cy- 4)*brick.nx, 0,
  6, (Cx+0) + (Cy- 5)*brick.nx, 0,
  6, (Cx+0) + (Cy- 6)*brick.nx, 0,
  6, (Cx+0) + (Cy- 7)*brick.nx, 0,
  6, (Cx+0) + (Cy- 8)*brick.nx, 0,
  6, (Cx+0) + (Cy- 9)*brick.nx, 0,
  6, (Cx+0) + (Cy-10)*brick.nx, 0,
  6, (Cx+0) + (Cy-11)*brick.nx, 0,
  6, (Cx+0) + (Cy-12)*brick.nx, 0,
  6, (Cx+0) + (Cy-13)*brick.nx, 0,
  6, (Cx+0) + (Cy-14)*brick.nx, 0,
  6, (Cx+0) + (Cy-15)*brick.nx, 0,
  6, (Cx+0) + (Cy-16)*brick.nx, 0,
  6, (Cx+0) + (Cy-17)*brick.nx, 0,
  6, (Cx+0) + (Cy-18)*brick.nx, 0,
  6, (Cx+0) + (Cy-19)*brick.nx, 0,
  6, (Cx+0) + (Cy-20)*brick.nx, 0,
  6, (Cx+0) + (Cy-21)*brick.nx, 0,
  6, (Cx+0) + (Cy-22)*brick.nx, 0,
  6, (Cx+0) + (Cy-23)*brick.nx, 0,
  6, (Cx+0) + (Cy-24)*brick.nx, 0,
}
