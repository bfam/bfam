-- simulation data
N = 2
num_subdomains = 2
connectivity = "brick"
lsrk_method = "KC54"
problem = "slip weakening"
refine_level = 3
max_refine_level = 8

-- 2 blocks with fault at y = 0
brick_m = 1
brick_n = 2
brick_a = 0
brick_b = 0

brick_Lx = 30
brick_Ly = brick_Lx

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
