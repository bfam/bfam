num_subdomains = 1
connectivity = "brick"
lsrk_method = "KC54"
problem = "slip weakening"
nsteps   = 12000
ndisp    = 1
noutput_body   = 1000
noutput_fault  = 10
nfaults = 4
angle = math.pi / 6
cos_angle = math.cos(angle)
sin_angle = math.sin(angle)
x_end = 12*cos_angle
y_end =-12*sin_angle

function dt_modify ( dt )
  noutput_fault = math.ceil(t_fault/dt)
  dt = t_fault/noutput_fault

  nsteps = math.ceil(t_final/dt)
  ndisp = math.ceil(t_disp/dt)
  noutput_body = math.ceil(t_body/dt)

  return  noutput_fault, noutput_body, ndisp, nsteps, dt
end

-- 2 blocks with fault at y = 0
brick_m = 5
brick_n = 5
brick_a = 0
brick_b = 0

function fault_bricks(f)
  if f == 0 then
    return 11,9
  elseif f == 1 then
    return 14,12
  elseif f == 2 then
    return 15,13
  elseif f == 3 then
    return 13,12
  end
end


function brick_corners(x, y, z )
  base = 12
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

  return zout,yout,xout
end

function fault_refine(x,y,z,level)
  dist = math.huge
  if x < -16 then
    dist = (x+16)^2+y^2
  elseif x > 12 then
    dist = (x-12)^2+y^2
  else
    dist = y^2
  end
  if y < 0 then
    xr = x*cos_angle - y*sin_angle
    yr = x*sin_angle + y*cos_angle
    if xr < 0 then
      dist = math.min(dist,xr^2+yr^2)
    elseif xr > 12 then
      dist = math.min(dist,(xr-12)^2+yr^2)
    else
      dist = math.min(dist,yr^2)
    end
  end
  dist = math.sqrt(dist)
  if(dist < refine_dist*(max_refine_level-level)) then
    return 1
  else
    return 0
  end
end

-- material properties
cs = 3.464
cp = 6
rho = 2.670
mu  = rho*cs^2
lam = rho*cp^2-2*mu

-- friction stuff
fault_fd  =    0.525
fault_Dc  =    0.4

function fault_fs ( t, x, y, z )
  if x > -16 and x < 12 and y > -12*sin_angle then
    return( 0.677 )
  else
    return( 10000 )
  end
end

function fault_Tp1_0 ( t, x, y, z, n1, n2, n3)
  if x > -11 and x < -8 then
    return 81.6*n2
  else
    return 70*n2
  end
end

function fault_Tp2_0 ( t, x, y, z, n1, n2, n3)
  return -70*n1
end

function fault_Tp3_0 ( t, x, y, z, n1, n2, n3)
  return 0
end

function fault_Tn_0 ( t, x, y, z, n1, n2, n3)
  return -120
end
