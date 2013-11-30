num_subdomains = 1
connectivity = "brick"
lsrk_method = "KC54"
problem = "slip weakening"
nsteps   = 12000
ndisp    = 1
noutput_body   = 1000
noutput_fault  = 10
nfaults = 2

function fault_bricks(f)
  if f == 0 then
    return 0,2
  elseif f == 1 then
    return 1,3
  end
end


function fault_refine(x,y,z,level)
  if x < -15 then
    dist = (x+15)^2+y^2
  elseif x > 15 then
    dist = (x-15)^2+y^2
  else
    dist = y^2
  end
  dist = math.sqrt(dist)
  if(dist < refine_dist*(max_refine_level-level)) then
    return 1
  else
    return 0
  end
end

function dt_modify ( dt )
  noutput_fault = math.ceil(t_fault/dt)
  dt = t_fault/noutput_fault

  nsteps = math.ceil(t_final/dt)
  ndisp = math.ceil(t_disp/dt)
  noutput_body = math.ceil(t_body/dt)
  return  noutput_fault, noutput_body, ndisp, nsteps, dt
end

-- 2 blocks with fault at y = 0
brick_m = 2
brick_n = 2
brick_a = 0
brick_b = 0

brick_Lx = 30
brick_Ly = 30

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
  if math.abs(x) <= 15 then
    return( 0.677 )
  else
    return( 10000 )
  end
end

function fault_Tp1_0 ( t, x, y, z, n1, n2, n3)
  if math.abs(x) <= 1.5 then
    return 81.6*n2
  elseif math.abs(x-7.5) <= 1.5 then
    return 62*n2
  elseif math.abs(x+7.5) <= 1.5 then
    return 78*n2
  else
    return 70*n2
  end
end

function fault_Tp2_0 ( t, x, y, z, n1, n2, n3)
  return 0
end

function fault_Tp3_0 ( t, x, y, z, n1, n2, n3)
  return 0
end

function fault_Tn_0 ( t, x, y, z, n1, n2, n3)
  return -120
end
