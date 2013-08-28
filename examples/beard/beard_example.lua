N = 5
num_subdomains = 1
connectivity = "brick"
lsrk_method = "KC54"
problem = "stress free box"
refine_level = 0
max_refine_level = 2

nsteps   = 10000
ndisp    = 10
noutput  = 10
-- noutput  = math.huge

dt_scale = 0.1

brick_m = 2
brick_n = 2
brick_a = 1
brick_b = 1
brick_random = 0

cs = 3.464
cp = 6
rho = 2.670
mu  = rho*cs^2
lam = rho*cp^2-2*mu
