-- default parameters
order = 4
min_level = 3
max_level = 3

-- store random seed
math.randomseed(0)

-- refinement parameters
data_directory = "data"
output_prefix = "solution"

-- connectivity info
connectivity = "brick"
brick =
{
  nx = 4,
  ny = 1,
  nz = 1,
  bc0 = 1,
  bc1 = 1,
  bc2 = 2,
  bc3 = 2,
  bc4 = 2,
  bc5 = 2,
}

-- set up the domain
Lx = 1
Ly = 1
Lz = 1

Cx = Lx*brick.nx/2
Cy = Ly*brick.ny/2
Cz = 0

Sw = 0.5*Lx
Sx = Cx-Sw

function connectivity_vertices(x, y, z)
  -- if x > 0 and x < brick.nx and
  --    y > 0 and y < brick.ny then
  --   x = x + 0.5*(math.random()-0.5)
  --   y = y + 0.5*(math.random()-0.5)
  -- end
  xout = Lx*x-Cx
  yout = Ly*y-Cy
  zout = Lz*z-Cz
  return xout,yout,zout
end

function refinement_function(
  x0,y0,z0,x1,y1,z1,
  x2,y2,z2,x3,y3,z3,
  level, treeid)

  if level < min_level then
    return 1
  else
    return 0
  end
end

function element_order(
  x0,y0,z0,x1,y1,z1,
  x2,y2,z2,x3,y3,z3,
  level, treeid)

  -- V = (x0+x1+x2+x3)/4
  V = math.abs(x0)
  V = math.max(V,x1)
  V = math.max(V,x2)
  V = math.max(V,x3)
  N = 4
  if( math.abs(V) > Cx-0.5*Lx ) then
    return N, "sponge"
  else
    return N, "elastic"
  end
end

-- initial values
rho = 1
lam = 1
mu  = 1
cp = math.sqrt((lam+2*mu)/rho)
S11 = 0
S22 = 0
S33 = 0
S12 = 0
S13 = 0
S23 = 0
v1  = 0
v2  = 0
v3  = 0

function v1(x,y,z,t)
  r2 = x^2+y^2+z^2
  return math.exp(-r2/0.1^2)
end

function a_sponge(x,y,z,t,xc,yc,zc)
  return 100*((math.abs(x)-Sx)/Sw)^2
end


-- time stepper to use
lsrk_method  = "KC54"

tend      = 10*Lx/cp
tout      = tend/100
tdisp     = tend/100
tstations = -1
tfout     = -1
terr      = tend
T         = tdisp

function time_step_parameters(dt)
  dt      = 0.5*dt

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

-- function nerr(dt)
--   return terr/dt
-- end

-- set the glue grids
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
