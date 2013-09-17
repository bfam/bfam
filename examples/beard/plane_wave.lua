-- refinement parameters
min_refine_level = 2
max_refine_level = 4

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

  if(level < treeid/4) then
    return 1
  else
    return 0
  end
end

-- time stepper to use
lsrk_method  = "KC54"

