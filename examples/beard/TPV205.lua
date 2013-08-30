-- simulation data
N = 4
N_fault = 2
num_subdomains = 1
connectivity = "brick"
lsrk_method = "KC54"
problem = "slip weakening"
refine_level = 5
max_refine_level = 9
nsteps   = 12000
ndisp    = 1
noutput_body   = 1000
noutput_fault  = 10

function dt_modify ( dt )
  noutput_fault = math.ceil(0.1/dt)
  dt = 0.05/noutput_fault
  nsteps = 12/dt
  noutput_body = 1/dt
  noutput_fault = 0.1/dt
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
