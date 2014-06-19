-- default parameters
N1 = 2
N2 = 2
N3 = 2
min_level = 1
max_level = 2

-- store random seed
math.randomseed(0)

-- refinement parameters
output_prefix = "solution_2d"

-- connectivity info
connectivity = "brick"
brick =
{
  nx = 4,
  ny = 3,
  nz = 2,
  periodic_x = 1,
  periodic_y = 1,
  periodic_z = 1,
}

-- set up the domain
Lx = 25
Ly = 25
Lz = 25
function connectivity_vertices(x, y, z)
  if x > 0 and x < brick.nx and
     y > 0 and y < brick.ny then
    x = x + 0.5*(math.random()-0.5)
    y = y + 0.5*(math.random()-0.5)
  end
  xout = Lx*x
  yout = Ly*y
  zout = 0
  return xout,yout,zout
end

function refinement_function(
  x0,y0,z0,x1,y1,z1,
  x2,y2,z2,x3,y3,z3,
  level, treeid)

  if level < min_level then
    return 1
  elseif level >= max_level or x0+x1-y0-y1 < 0 then
    return 0
  end
  return 0
end

function element_order(
  x0,y0,z0,x1,y1,z1,
  x2,y2,z2,x3,y3,z3,
  level, treeid)

  if treeid%3 == 0 then
    return N1
  elseif treeid%3 == 1 then
    return N2
  end
  return N3
end

-- time stepper to use
lsrk_method  = "KC54"

tend  = 4*Lx
tout  = 2*tend
tdisp = 2*tend
terr  = tend
