N = 4
num_subdomains = 1
connectivity = "brick"
lsrk_method = "KC54"
problem = "stress free box"
refine_level = 3
max_refine_level = 3

nsteps   = 10000
ndisp    = 10
noutput  = 10
-- noutput  = math.huge

dt_scale = 1

brick_m = 2
brick_n = 2
brick_a = 1
brick_b = 1
brick_random = 1

cs = 3.464
cp = 6
rho = 2.670
mu  = rho*cs^2
lam = rho*cp^2-2*mu

function fault_refine(x,y,z,level)
  return 0
end

function dt_modify ( dt )
  tf = 1
  t_final = tf/cs/math.sqrt(2)
  t_body  = t_final/tf/10
  noutput_body = math.ceil(t_body/dt)
  dt = t_body/noutput_body
  nsteps = math.ceil(t_final/dt)
  ndisp = math.ceil(t_body/dt)
  noutput_fault = 2*noutput_body
  return  noutput_fault, noutput_body, ndisp, nsteps, dt
end

