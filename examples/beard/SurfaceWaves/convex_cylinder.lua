pi   = math.pi
sqrt = math.sqrt
abs  = math.abs
cos = math.cos
sin = math.sin

output_prefix = "convex_cylinder"
data_directory = "data"
elem_order = 4
max_level = 0
min_level = 0
static_refinement=3

mesh_file = "/home/jekozdon/codes/bfam/examples/beard/SurfaceWaves/convex_cylinder.inp"

default_boundary_tag = "free surface"

c1 = sqrt(pi/5)/2
c2 = 1/sqrt(2)

function connectivity_vertices(xin, yin, zin)

  z = 0

  if abs(xin) == 1 then
    x = (xin/abs(xin))*c1
    y = (yin/abs(yin))*c1
  elseif abs(xin) == 2 then
    x = (xin/abs(xin))*c2
    y = (yin/abs(yin))*c2
  end

  return x, y, z
end

function transform_nodes(x, y, z)
  local x1 = x
  local y1 = y
  local  s = 0
  local  r = sqrt(x^2+y^2)
  local x2 = 0
  local y2 = 0

  if     abs(x) <= abs(y) and y >= c1 then
    -- NORTH
    s = (y-c1)/(c2-c1)

    x1 = (c1/y)*x
    y1 = c1

    x2 = x/r
    y2 = y/r
  elseif abs(x) <= abs(y) and y <= -c1 then
    -- SOUTH
    s = (abs(y)-c1)/(c2-c1)

    x1 = -(c1/y)*x
    y1 = -c1

    x2 = x/r
    y2 = y/r
  elseif abs(y) <= abs(x) and x >=  c1 then
    -- WEST
    s = (x-c1)/(c2-c1)

    x1 = c1
    y1 = (c1/x)*y

    x2 = x/r
    y2 = y/r
  elseif abs(y) <= abs(x) and x <= -c1 then
    -- EAST
    s = (abs(x)-c1)/(c2-c1)

    x1 = -c1
    y1 = -(c1/x)*y

    x2 = x/r
    y2 = y/r
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
MU = {1,0.1,0.01,0.001}

OMEGA = {{6.712844035888430,10.598136345660269,13.113329869474558,16.021054935234528,17.730944523371178},
         {2.174059402269135, 3.475669329435217, 4.373718805991318, 5.375869130532646, 6.404312117277238},
         {0.690252461874728, 1.104283873555503, 1.392820917483943, 1.707902134335568, 2.031628962752277},
         {0.218370110186628, 0.349373578609341, 0.440766027573439, 0.540340250707020, 0.642652150833394}}

BI = {{0.259859421131935    ,0.551757806713995    ,-1.111726166772277    ,0.958432243767886    ,-0.732513991983086    },
      {0.008979356435568    ,0.061797475082992    ,-0.195334916685379    ,0.716217080313942    ,-1.855812650827352    },
      {1.800969171774775e-05,1.608778558352172e-04,-5.996559491070863e-04,0.002791241585515    ,-0.010369491587158    },
      {1.946610327968597e-08,1.791032555369920e-07,-6.806632052120861e-07,3.250907109941674e-06,-1.252501580039413e-05}}

KI = 1
KJ = 1

mu    = MU[KI]

lam = 1
rho = 1
cp  = sqrt((lam+2*mu)/rho)
cs  = sqrt(mu/rho)
CW = {
  omega = OMEGA[KI][KJ],
  bi    = BI[KI][KJ],
  A     = 1,
  n     = 6,
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
