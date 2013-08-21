-- simulation data
N = 5
num_subdomains = 2
connectivity = "brick"
lsrk_method = "KC54"
problem = "slip weakening"
refine_level = 5

-- 2 blocks with fault at y = 0
brick_m = 1
brick_n = 2
brick_a = 0
brick_b = 0

brick_Lx = 3
brick_Ly = brick_Lx

-- material properties
cs = 3.464
cp = 6
rho = 2.670
mu  = rho*cs^2
lam = rho*cp^2-2*mu

-- friction stuff
friction_fs  =    0.67
friction_fd  =    0.525
friction_Dc  =    0.4
friction_S12 =   62.0
friction_S22 = -120

friction_nuc_dS12 = 81.6-friction_S12
friction_nuc_x   =  0;
friction_nuc_y   =  0;
friction_nuc_z   =  0;
friction_nuc_R   =  0.3;

