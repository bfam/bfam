pi   = math.pi
sqrt = math.sqrt
abs  = math.abs
cos = math.cos
sin = math.sin

output_prefix = "concave_cylinder"
data_directory = "data"
elem_order = 12
max_level = 0
min_level = 0
static_refinement=3

mesh_file = "/home/jekozdon/codes/bfam/examples/beard/SurfaceWaves/concave_cylinder.inp"

default_boundary_tag = "free surface"

c1 = sqrt(pi/5)/2
c2 = 1/sqrt(2)
c3 = 3

function connectivity_vertices(xin, yin, zin)

  z = 0

  if abs(xin) == 1 then
    x = (xin/abs(xin))*c1
  elseif abs(xin) == 2 then
    x = (xin/abs(xin))*c2
  elseif abs(xin) == 3 then
    x = (xin/abs(xin))*c3
  end

  if abs(yin) == 1 then
    y = (yin/abs(yin))*c1
  elseif abs(yin) == 2 then
    y = (yin/abs(yin))*c2
  elseif abs(yin) == 3 then
    y = (yin/abs(yin))*c3
  end

  return x, y, z
end

function transform_nodes(x, y, z)
  local x1 = x
  local y1 = y
  local x2 = x
  local y2 = y
  local  s = 0

  if     abs(x) <= c2 and y >= abs(x) then
    -- NORTH
    s  = (y-c2)/(c3-c2)
    y1 = sqrt(1-x^2)
    y2 = c3
  elseif abs(x) <= c2 and y <= -abs(x) then
    -- SOUTH
    s  = (-y-c2)/(c3-c2)
    y1 = -sqrt(1-x^2)
    y2 = -c3
  elseif abs(y) <= c2 and x <= -abs(y) then
    -- WEST
    s  = (-x-c2)/(c3-c2)
    x1 = -sqrt(1-y^2)
    x2 = -c3
  elseif abs(y) <= c2 and x >= abs(y) then
    -- EAST
    s  = (x-c2)/(c3-c2)
    x1 = sqrt(1-y^2)
    x2 = c3
  end

  x = (1-s)*x1 + s*x2
  y = (1-s)*y1 + s*y2

  return x, y, z
end

function element_size(
  x0,y0,z0,x1,y1,z1,
  x2,y2,z2,x3,y3,z3)


  h = (x0-x1)^2+(y0-y1)^2+(z0-z1)^2
  hmin = h
  hmax = h

  h = (x0-x2)^2+(y0-y2)^2+(z0-z2)^2
  hmax = math.max(hmax,h)
  hmin = math.min(hmin,h)

  h = (x1-x3)^2+(y1-y3)^2+(z1-z3)^2
  hmax = math.max(hmax,h)
  hmin = math.min(hmin,h)

  h = (x2-x3)^2+(y2-y3)^2+(z2-z3)^2
  hmax = math.max(hmax,h)
  hmin = math.min(hmin,h)

  return sqrt(hmin), sqrt(hmax)

end

function refinement_function(
  x0,y0,z0,x1,y1,z1,
  x2,y2,z2,x3,y3,z3,
  level, treeid)

  if(level < min_level) then
    return 1
  end

  return 0
end

function element_order(
  x0,y0,z0,x1,y1,z1,
  x2,y2,z2,x3,y3,z3,
  level, treeid)

  return elem_order
end

-- material properties

OMEGA_R = {  0.495719213980678,  1.376300729324809,  2.247966737740946,   3.122343319168747,   4.000739349461448,   4.883072039743391,   5.768949287751974}
OMEGA_I = {  0.482489802859255,  0.562001476905459,  0.595314872100875,   0.610590759912388,   0.616417207729038,   0.616677688937145,   0.613435428254969}
B_R     =  {-0.119481693313368, -0.433288418131060, -0.908447005479862,  -1.634838389505964,  -2.741262468068052,  -4.415934969322963,  -6.934676718614770}
B_I     =  {-2.996924150548532, -4.741971316305029, -7.300093250067868, -11.096286394747006, -16.740023608389567, -25.126574554563344, -37.575882485238743}

KI = 1

mu  = 1
lam = 1
rho = 1
cp  = sqrt((lam+2*mu)/rho)
cs  = sqrt(mu/rho)
CW = {
  omega_r = OMEGA_R[KI],
  omega_i = OMEGA_I[KI],
  b_y     = B_R[KI],
  b_i     = B_I[KI],
  A     = 1,
  n       = KI+1,
  Ka = OMEGA[KI][KJ]/cp,
  Kb = OMEGA[KI][KJ]/cs,
}


-- time stepper to use
lsrk_method  = "KC54"

tend   = 2*pi/CW.omega
tend = 2
tout   = tend
tdisp  = tend
-- nerr   = 0
function nerr(dt)
  return terr/dt
end
terr  = tend

T = tout

function time_step_parameters(dt)
  dt = 0.5*dt
  N  = math.ceil(T / dt)
  dt = T / N

  nsteps     = tend / dt
  ndisp      = tdisp / dt
  noutput    = tout  / dt
  nstations  = -1
  nfoutput   = -1

  nsteps = 1
  noutput = 1

  return dt,nsteps, ndisp, noutput, nfoutput, nstations
end

bc_free = {
  type = "boundary",
  tag  = "free surface",
}

bc_nonreflect = {
  type = "boundary",
  tag  = "non-reflecting",
}

bc_rigid = {
  type = "boundary",
  tag  = "rigid wall",
}

glue_info = {
  bc_nonreflect,
  bc_free,
  bc_rigid,
}




-- field conditions
S11 = 0
S22 = 0
S33 = 0
S12 = 0
S13 = 0
S23 = 0
v2  = 0
v3  = 0

function CW_prelims(x,y,z,t)
  local r = sqrt(x^2+y^2)
  if r == 0 then
    return 0,0,0,0
  end

  local S = math.atan2(y,x)

  local A  = CW.A
  local bi = CW.bi
  local n  = CW.n
  local Ka = CW.Ka
  local Kb = CW.Kb



  local q1 =  A*Ka/2*(jn(n-1,Ka*r)-jn(n+1,Ka*r));
  local q2 =-bi*(n/r)*jn(n,Kb*r);
  local w1 = A*(n/r)*jn(n,Ka*r);
  local w2 = -bi*Kb/2*(jn(n-1,Kb*r)-jn(n+1,Kb*r));

  local W = w1+w2;
  local Q = q1+q2;
  return r,S,W,Q
end

function djn(n,r)
  return 0.5*(jn(n-1,r) - jn(n+1,r))
end

function CW_deriv_prelims(x,y,z,t)
  local r,S,W,Q = CW_prelims(x,y,z,t)
  local A     = CW.A
  local bi    = CW.bi
  local n     = CW.n
  local Ka    = CW.Ka
  local Kb    = CW.Kb
  local omega = CW.omega

  local r_x = x/r
  local r_y = y/r
  local S_x = -y/r^2
  local S_y =  x/r^2

  local q1_r =  A*Ka^2/2*(djn(n-1,Ka*r)-djn(n+1,Ka*r));
  local q2_r = bi*(n/r^2)*jn(n,Kb*r) - Kb*bi*(n/r)*djn(n,Kb*r);
  local w1_r =-A*(n/r^2)*jn(n,Ka*r) + Ka*A*(n/r)*djn(n,Ka*r);
  local w2_r = -bi*Kb^2/2*(djn(n-1,Kb*r)-djn(n+1,Kb*r));

  if r== 0 then
    r_x = 0
    r_y = 0
    S_x = 0
    S_y = 0
    q2_r = 0
    w1_r = 0
  end

  local Q_r = q1_r+q2_r;
  local W_r = w1_r+w2_r;
  local ux_S =   -Q*cos(n*S + omega*t)*sin(S) +  n*W*sin(S)*cos(n*S + omega*t) -
                n*Q*sin(n*S + omega*t)*cos(S) +    W*cos(S)*sin(n*S + omega*t);
  local ux_Q = cos(S) * cos(omega*t+n*S);
  local ux_W = sin(S) * sin(omega*t+n*S);
  local ux_r = ux_Q * Q_r + ux_W * W_r;

  local ux_x = ux_r*r_x + ux_S*S_x;
  local ux_y = ux_r*r_y + ux_S*S_y;

  local uy_S =    Q*cos(S)*cos(n*S + omega*t) - n*W*cos(S)*cos(n*S + omega*t) -
               n*Q*sin(S)*sin(n*S + omega*t) +    W*sin(S)*sin(n*S + omega*t);

  local uy_Q = sin(S) * cos(omega*t+n*S);
  local uy_W = - cos(S) * sin(omega*t+n*S);
  local uy_r = uy_Q * Q_r + uy_W * W_r;

  local uy_x = uy_r*r_x + uy_S*S_x;
  local uy_y = uy_r*r_y + uy_S*S_y;

  return ux_x, ux_y, uy_x, uy_y
end


function v1(x,y,z,t)
  local r,S,W,Q = CW_prelims(x,y,z,t)
  local omega   = CW.omega
  local n       = CW.n
  return -omega * cos(S)*Q * sin(omega*t+n*S) +
          omega * sin(S)*W * cos(omega*t+n*S)
end

function v2(x,y,z,t)
  local r,S,W,Q = CW_prelims(x,y,z,t)
  local omega   = CW.omega
  local n       = CW.n
  return -omega * sin(S)*Q * sin(omega*t+n*S) -
          omega * cos(S)*W * cos(omega*t+n*S)
end

function S11(x,y,z,t)
  local ux_x, ux_y, uy_x, uy_y = CW_deriv_prelims(x,y,z,t)
  return lam*(ux_x + uy_y) + 2*mu*ux_x;
end

function S22(x,y,z,t)
  local ux_x, ux_y, uy_x, uy_y = CW_deriv_prelims(x,y,z,t)
  return lam*(ux_x + uy_y) + 2*mu*uy_y;
end

function S12(x,y,z,t)
  local ux_x, ux_y, uy_x, uy_y = CW_deriv_prelims(x,y,z,t)
  return mu*(ux_y+uy_x);
end

function v3(x,y,z,t)
  return 0
end

function S33(x,y,z,t)
  local ux_x, ux_y, uy_x, uy_y = CW_deriv_prelims(x,y,z,t)
  return lam*(ux_x + uy_y)
end

function S23(x,y,z,t)
  return 0
end

function S13(x,y,z,t)
  return 0
end
