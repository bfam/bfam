-- refinement parameters
max_level = 3

-- connectivity info
connectivity = "brick"
brick =
{
  nx = 5,
  ny = 5,
  periodic_x = 0,
  periodic_y = 0,
}

function connectivity_vertices(x, y, z)
  base = 12
  angle = math.pi / 6
  cos_angle = math.cos(angle)
  sin_angle = math.sin(angle)
  zout =  0

  xout = -3*base + x*base
  yout = -3*base + y*base

  if x == 5 then
    xout =  base*(x-3+cos_angle)
  elseif x >= 3 and y <= 2 then
    xout =  base*(x-3+cos_angle)
  end

  if y <= 2 then
    yout = -base*(2-y+sin_angle)
  end

  return xout,yout,zout
end

function refinement_function(
  x0,y0,z0,x1,y1,z1,
  x2,y2,z2,x3,y3,z3,
  level, treeid)

  if level >= max_level or level > treeid/4 then
    return 0
  else
    return 1
  end
end

function element_order(
  x0,y0,z0,x1,y1,z1,
  x2,y2,z2,x3,y3,z3,
  level, treeid)

  N = treeid%3+1

  return N
end

-- initial values
rho = 1
lam = 2
mu  = 3
v1  = 4
v2  = 5
v3  = 6
S11 = 7
S22 = 8
S33 = 9
S12 = 10
S13 = 11
S23 = 12


-- time stepper to use
lsrk_method  = "KC54"

function time_step_parameters(dt)
  tout  = 0.1
  tdisp = 0.1
  tend  = 1;

  dt      = 0.5*dt
  noutput = math.ceil(tout / dt)
  dt      = tout / noutput
  ndisp   = tdisp / dt
  nstesp  = tend / dt
  return dt,nstesp, ndisp, noutput
end
