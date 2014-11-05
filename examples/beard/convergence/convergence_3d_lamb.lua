-- default convergence parameters
min_level = 0
max_level = min_level
static_refinement = 1
N1 = 4
N2 = 4
N3 = 4

-- convenience functions
printf = function(s,...)
  return io.write(s:format(...))
end -- function

-- store random seed
math.randomseed(0)
-- assume using double for lua numbers
eps = 2.220446049250313e-16

-- refinement parameters
output_prefix = "solution"

-- connectivity info
connectivity = "brick"
scale = 2
brick =
{
  nx = 2*scale,
  ny = 1*scale,
  nz = 1*scale,
  periodic_x = 1,
  periodic_y = 0,
  periodic_z = 1,
  p4est_brick = 1,
}

default_boundary_tag = "free surface"

-- set up the domain
function connectivity_vertices(x, y, z)
  if x > 0 and x < brick.nx and
     y > 0 and y < brick.ny and
     z > 0 and z < brick.nz then
    x = x + 0.5*(math.random()-0.5)
    y = y + 0.5*(math.random()-0.5)
    z = z + 0.5*(math.random()-0.5)
  end
  x = x/scale - 1.0
  y = y/scale - 0.5
  z = z/scale - 0.5
  return x,y,z
end

function refinement_function(
  x0,y0,z0,x1,y1,z1,
  x2,y2,z2,x3,y3,z3,
  x4,y4,z4,x5,y5,z5,
  x6,y6,z6,x7,y7,z7,
  level, treeid)

  if level < min_level then
    return 1
  elseif level >= max_level or x0+x1-y0-y1 < 0 then
    return 0
  else
    return 1
  end
end

function element_order(
  x0,y0,z0,x1,y1,z1,
  x2,y2,z2,x3,y3,z3,
  x4,y4,z4,x5,y5,z5,
  x6,y6,z6,x7,y7,z7,
  level, treeid)

  if treeid%3 == 0 then
    return N1
  elseif treeid%3 == 1 then
    return N2
  end
  return N3
end

-- initial values
rho = 1
lam = 2
mu  = 1
wavelength = 1.0
ymax = 0.5 -- assumed that ymin = -ymax

-- derived values
cpsq = (lam + 2.0 * mu) / rho
cssq = mu / rho
cp = math.sqrt(cpsq)
cs = math.sqrt(cssq)
k = 2.0 * math.pi / wavelength

function lamb_f(omega)
  p = math.sqrt(omega * omega / cpsq - k * k)
  q = math.sqrt(omega * omega / cssq - k * k)

  assert(omega * omega / cpsq - k * k >= -eps)
  assert(omega * omega / cssq - k * k >= -eps)

  f = 4.0 * mu * k * k * p * q * math.sin (p * ymax) * math.cos (q * ymax)
    - (-k * k + q * q) * (-lam * k * k - (lam + 2.0 * mu) * p * p)
    * math.cos (p * ymax) * math.sin (q * ymax)

  return f
end

-- Use bisection method to calculate omega
omega_step = k/10.0
omega_max_steps = 1000
omega_bi_steps  = 1000

omega_min = k * cp + 0.01
f_omega_min = lamb_f(omega_min)

-- Find a positive and negative interval for f(omega) */
for index = 1,omega_max_steps do
  omega_max = index * omega_step + omega_min
  f_omega_max = lamb_f(omega_max)

  if f_omega_min*f_omega_max < -eps*eps  then
    printf("Lamb: index_interval = %d\n", index)
    break
  end

  assert(index ~= omega_max_steps, "Interval not found for bisection")
end


printf("Lamb: omega_min      = %25.16e\n", omega_min)
printf("Lamb: omega_max      = %25.16e\n", omega_max)
printf("Lamb: Start Bisection\n");

for index = 1,omega_bi_steps do
  omega = 0.5 * (omega_min + omega_max)
  f_omega = lamb_f(omega)

  if math.abs (f_omega) < 1e4 * eps then
    printf ("Lamb: index = %d\n", index);
    break
  end

  if f_omega * f_omega_min < -eps * eps then
      omega_max = omega
      f_omega_max = f_omega
  elseif f_omega * f_omega_max < -eps * eps then
    omega_min = omega
    f_omega_min = f_omega
  else
    printf("Lamb: Error: f_omega_min = %25.16e\n", f_omega_min)
    printf("Lamb: Error: f_omega_max = %25.16e\n", f_omega_max)
    printf("Lamb: Error: f_omega     = %25.16e\n", f_omega)
    asset(0, "Should not be reached")
  end

  assert(index ~= omega_bi_steps, "Bisection not converged");
end

assert(omega * omega / cpsq - k * k > 0.0)
assert(omega * omega / cssq - k * k > 0.0)

p = math.sqrt(omega * omega / cpsq - k * k)
q = math.sqrt(omega * omega / cssq - k * k)

A = 2.0 * mu * k * q * math.cos (q * ymax)
B = (lam * k * k + (lam + 2.0 * mu) * p * p) * math.cos (p * ymax)


printf("Lamb: omega_min   = %25.16e\n", omega_min);
printf("Lamb: omega_max   = %25.16e\n", omega_max);
printf("Lamb: omega       = %25.16e\n", omega);

printf("Lamb: f_omega_min = %25.16e\n", f_omega_min);
printf("Lamb: f_omega_max = %25.16e\n", f_omega_max);
printf("Lamb: f_omega     = %25.16e\n", f_omega);

printf("Lamb: k    = %25.16e\n", k);
printf("Lamb: mu   = %25.16e\n", mu);
printf("Lamb: lam  = %25.16e\n", lam);
printf("Lamb: rho  = %25.16e\n", rho);

printf("Lamb: A    = %25.16e\n", A);
printf("Lamb: B    = %25.16e\n", B);
printf("Lamb: p    = %25.16e\n", p);
printf("Lamb: q    = %25.16e\n", q);
printf("Lamb: ymax = %25.16e\n", ymax);

function v(x1,x2,x3,t,i)
  sink =      math.sin(k * x1 - omega * t)
  cosk =      math.cos(k * x1 - omega * t)
  Asinp = A * math.sin(p * x2)
  Acosp = A * math.cos(p * x2)
  Bsinq = B * math.sin(q * x2)
  Bcosq = B * math.cos(q * x2)

  if     i == 1 then return omega * ( k * Acosp + q * Bcosq) * cosk
  elseif i == 2 then return omega * (-p * Asinp + k * Bsinq) * sink
  elseif i == 3 then return 0.0
  else
    assert(0, "Never should be reached")
  end
end

function v1(x,y,z,t)
  return v(x,y,z,t,1)
end
function v2(x,y,z,t)
  return v(x,y,z,t,2)
end
function v3(x,y,z,t)
  return v(x,y,z,t,3)
end

function S(x1,x2,x3,t,i,j)
  sink =      math.sin(k * x1 - omega * t)
  cosk =      math.cos(k * x1 - omega * t)
  Asinp = A * math.sin(p * x2)
  Acosp = A * math.cos(p * x2)
  Bsinq = B * math.sin(q * x2)
  Bcosq = B * math.cos(q * x2)

  E11 =     k * (-k * Acosp - q * Bcosq) * cosk
  E22 = (-p * p * Acosp + k * q * Bcosq) * cosk
  E33 = 0.0
  E12 = 0.5 * ((k * p * Asinp + q * q * Bsinq) * sink
                - k * (-p * Asinp + k * Bsinq) * sink)
  E13 = 0.0
  E23 = 0.0

  trE = E11 + E22 + E33

  if     i == 1 and j == 1 then
    return lam*trE + 2*mu * E11
  elseif i == 2 and j == 2 then
    return lam*trE + 2*mu * E22
  elseif i == 3 and j == 3 then
    return lam*trE + 2*mu * E33
  elseif i == 1 and j == 2  then
    return           2*mu * E12
  elseif i == 1 and j == 3  then
    return           2*mu * E13
  elseif i == 2 and j == 3  then
   return            2*mu * E23
  else
    assert(0, "Never should be reached")
  end
end

function S11(x,y,z,t)
  return S(x,y,z,t,1,1)
end
function S22(x,y,z,t)
  return S(x,y,z,t,2,2)
end
function S33(x,y,z,t)
  return S(x,y,z,t,3,3)
end
function S12(x,y,z,t)
  return S(x,y,z,t,1,2)
end
function S13(x,y,z,t)
  return S(x,y,z,t,1,3)
end
function S23(x,y,z,t)
  return S(x,y,z,t,2,3)
end



-- time stepper to use
-- lsrk_method  = "KC54"

tend  = brick.ny/cs
tout  = 2*tend
tdisp = 2*tend
terr  = 0.01*tend

function time_step_parameters(dt)
  dt      = dt_fudge*dt
  nsteps = math.ceil(tend / dt)
  dt      = tend / nsteps
  ndisp   = tdisp / dt
  noutput  = 0
  nfoutput = -1
  nstations = -1
  return dt,nsteps, ndisp, noutput, nfoutput, nstations
end

function nerr(dt)
  return terr/dt
end
