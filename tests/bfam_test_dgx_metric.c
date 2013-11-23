#include <bfam.h>

#define REAL_APPROX_EQ(x, y, K)                                              \
  BFAM_APPROX_EQ((x), (y), (K), BFAM_REAL_ABS, BFAM_REAL_EPS, 10*BFAM_REAL_EPS)

static int
check_field(bfam_subdomain_t *s, const char *name, const bfam_real_t *ex)
{

  bfam_subdomain_dgx_t *sub = (bfam_subdomain_dgx_t *) s;
  int failures = 0;
  bfam_real_t *field =
    bfam_dictionary_get_value_ptr(&sub->base.fields, name);

  for(bfam_locidx_t i = 0; i < sub->K; ++i)
  {
    for(int j=0; j<sub->Np; ++j)
    {
      size_t idx = i*sub->Np + j;

      int fail = !REAL_APPROX_EQ(field[idx], ex[idx], 10);

      if(fail)
        BFAM_INFO("Fail match (%s) %d %25.15"BFAM_REAL_PRIe
            " %25.15"BFAM_REAL_PRIe " %d", name, idx,
            field[idx], ex[idx],
            BFAM_REAL_ABS(field[idx]-ex[idx]) < BFAM_REAL_MIN);

      failures += fail;
    }
  }

  return failures;
}

int
test_0d()
{
  return 0;
}

int
test_1d()
{
  const bfam_long_real_t  Vx[] = {1,-5};
  const bfam_long_real_t  Vy[] = { 0,1};
  const bfam_long_real_t  Vz[] = {-1,20};
  const bfam_long_real_t *Vi[] = {Vx,Vy,Vz};

  const bfam_real_t      x0r0[] = {0.5*(Vx[1]-Vx[0]),
                                    0.5*(Vx[1]-Vx[0]),
                                    0.5*(Vx[1]-Vx[0])};
  const bfam_real_t      x1r0[] = {0.5*(Vy[1]-Vy[0]),
                                    0.5*(Vy[1]-Vy[0]),
                                    0.5*(Vy[1]-Vy[0])};
  const bfam_real_t      x2r0[] = {0.5*(Vz[1]-Vz[0]),
                                    0.5*(Vz[1]-Vz[0]),
                                    0.5*(Vz[1]-Vz[0])};
  const bfam_real_t      *xr[] = {x0r0,x1r0,x2r0};
  const bfam_real_t      Jex[] = {BFAM_REAL_ABS(x0r0[0]),
                                  BFAM_REAL_ABS(x0r0[1]),
                                  BFAM_REAL_ABS(x0r0[2])};
  const bfam_real_t      xex[] = {(bfam_real_t)(Vx[0]),
                                  (bfam_real_t)(0.5*(Vx[0]+Vx[1])),
                                  (bfam_real_t)(Vx[1])};
  const bfam_real_t      yex[] = {(bfam_real_t)(Vy[0]),
                                  (bfam_real_t)(0.5*(Vy[0]+Vy[1])),
                                  (bfam_real_t)(Vy[1])};
  const bfam_real_t      zex[] = {(bfam_real_t)(Vz[0]),
                                  (bfam_real_t)(0.5*(Vz[0]+Vz[1])),
                                  (bfam_real_t)(Vz[1])};
  const bfam_real_t     *vex[] = {xex,yex,zex};

  const bfam_locidx_t EToV[] = {0,1};
  const bfam_locidx_t EToE[] = {0,0};
  const int8_t        EToF[] = {0,1};

  /* first check 1d */
  int failures = 0;
  for(int d = 0;d < 3;d++)
  {
    bfam_subdomain_dgx_t *d1 =
      bfam_subdomain_dgx_new_1(0, "1d", 2, 2, d+1, Vi, 1, EToV, EToE, EToF, 1);
    for(int v = 0;v < d+1;v++)
    {
      char name[BFAM_BUFSIZ];
      snprintf(name,BFAM_BUFSIZ,"_grid_x%d",v);
      failures += check_field((bfam_subdomain_t*)d1, name, vex[v]);
      if(d > 0)
      {
        snprintf(name,BFAM_BUFSIZ,"_grid_x%dr0",v);
        failures += check_field((bfam_subdomain_t*)d1, name, xr[v]);
      }
    }
    if(d == 0)
    {
      failures += check_field((bfam_subdomain_t*)d1, "_grid_Jr0x0", x0r0);
      failures += check_field((bfam_subdomain_t*)d1, "_grid_J",     Jex);
    }
    bfam_subdomain_free((bfam_subdomain_t*)d1);
    bfam_free(d1);
  }
  return 0;
}

int
test_2d()
{
  int failures = 0;
  const bfam_locidx_t   EToV[] = {0,1,2,3};
  const bfam_locidx_t   EToE[] = {0,0,0,0};
  const int8_t          EToF[] = {0,1,2,3};
  const bfam_locidx_t      K   =  1;
  const bfam_long_real_t  Vx[] = {   0, 0.2, -0.1,0.4};
  const bfam_long_real_t  Vy[] = {   0,-0.1,  0.1,0.4};
  const bfam_long_real_t  Vz[] = {-0.1,-0.2,  0.2,0.1};
  const bfam_long_real_t *Vi[] = {Vx,Vy,Vz};
  const int               nV   = 4;
  const int                N   = 2;

  const bfam_real_t      xex[] = {     Vx[0]       ,  0.5*(Vx[0]+Vx[1]            ),     Vx[1],
                                  0.5*(Vx[0]+Vx[2]), 0.25*(Vx[0]+Vx[1]+Vx[2]+Vx[3]),0.5*(Vx[1]+Vx[3]),
                                             Vx[2] ,  0.5*(            Vx[2]+Vx[3]),           Vx[3]};
  const bfam_real_t      yex[] = {     Vy[0]       ,  0.5*(Vy[0]+Vy[1]            ),     Vy[1],
                                  0.5*(Vy[0]+Vy[2]), 0.25*(Vy[0]+Vy[1]+Vy[2]+Vy[3]),0.5*(Vy[1]+Vy[3]),
                                             Vy[2] ,  0.5*(            Vy[2]+Vy[3]),           Vy[3]};
  const bfam_real_t      zex[] = {     Vz[0]       ,  0.5*(Vz[0]+Vz[1]            ),     Vz[1],
                                  0.5*(Vz[0]+Vz[2]), 0.25*(Vz[0]+Vz[1]+Vz[2]+Vz[3]),0.5*(Vz[1]+Vz[3]),
                                             Vz[2] ,  0.5*(            Vz[2]+Vz[3]),           Vz[3]};
  const bfam_real_t     *vex[] = {xex,yex,zex};

  const bfam_real_t       xr[] = {0.50*(Vx[1]      -Vx[0]      ),0.50*(Vx[1]      -Vx[0]      ),0.50*(Vx[1]      -Vx[0]      ),
                                  0.25*(Vx[1]+Vx[3]-Vx[0]-Vx[2]),0.25*(Vx[1]+Vx[3]-Vx[0]-Vx[2]),0.25*(Vx[1]+Vx[3]-Vx[0]-Vx[2]),
                                  0.50*(      Vx[3]      -Vx[2]),0.50*(      Vx[3]      -Vx[2]),0.50*(      Vx[3]      -Vx[2])};
  const bfam_real_t       xs[] = {0.50*(Vx[2]      -Vx[0]      ),0.25*(Vx[2]+Vx[3]-Vx[0]-Vx[1]),0.50*(      Vx[3]      -Vx[1]),
                                  0.50*(Vx[2]      -Vx[0]      ),0.25*(Vx[2]+Vx[3]-Vx[0]-Vx[1]),0.50*(      Vx[3]      -Vx[1]),
                                  0.50*(Vx[2]      -Vx[0]      ),0.25*(Vx[2]+Vx[3]-Vx[0]-Vx[1]),0.50*(      Vx[3]      -Vx[1])};

  const bfam_real_t       yr[] = {0.50*(Vy[1]      -Vy[0]      ),0.50*(Vy[1]      -Vy[0]      ),0.50*(Vy[1]      -Vy[0]      ),
                                  0.25*(Vy[1]+Vy[3]-Vy[0]-Vy[2]),0.25*(Vy[1]+Vy[3]-Vy[0]-Vy[2]),0.25*(Vy[1]+Vy[3]-Vy[0]-Vy[2]),
                                  0.50*(      Vy[3]      -Vy[2]),0.50*(      Vy[3]      -Vy[2]),0.50*(      Vy[3]      -Vy[2])};
  const bfam_real_t       ys[] = {0.50*(Vy[2]      -Vy[0]      ),0.25*(Vy[2]+Vy[3]-Vy[0]-Vy[1]),0.50*(      Vy[3]      -Vy[1]),
                                  0.50*(Vy[2]      -Vy[0]      ),0.25*(Vy[2]+Vy[3]-Vy[0]-Vy[1]),0.50*(      Vy[3]      -Vy[1]),
                                  0.50*(Vy[2]      -Vy[0]      ),0.25*(Vy[2]+Vy[3]-Vy[0]-Vy[1]),0.50*(      Vy[3]      -Vy[1])};

  const bfam_real_t       zr[] = {0.50*(Vz[1]      -Vz[0]      ),0.50*(Vz[1]      -Vz[0]      ),0.50*(Vz[1]      -Vz[0]      ),
                                  0.25*(Vz[1]+Vz[3]-Vz[0]-Vz[2]),0.25*(Vz[1]+Vz[3]-Vz[0]-Vz[2]),0.25*(Vz[1]+Vz[3]-Vz[0]-Vz[2]),
                                  0.50*(      Vz[3]      -Vz[2]),0.50*(      Vz[3]      -Vz[2]),0.50*(      Vz[3]      -Vz[2])};
  const bfam_real_t       zs[] = {0.50*(Vz[2]      -Vz[0]      ),0.25*(Vz[2]+Vz[3]-Vz[0]-Vz[1]),0.50*(      Vz[3]      -Vz[1]),
                                  0.50*(Vz[2]      -Vz[0]      ),0.25*(Vz[2]+Vz[3]-Vz[0]-Vz[1]),0.50*(      Vz[3]      -Vz[1]),
                                  0.50*(Vz[2]      -Vz[0]      ),0.25*(Vz[2]+Vz[3]-Vz[0]-Vz[1]),0.50*(      Vz[3]      -Vz[1])};
  const bfam_real_t    *xrex[] = {xr,xs,yr,ys,zr,zs};

  bfam_real_t J[9],Jrx[9],Jry[9],Jsx[9],Jsy[9];
  for(int n = 0; n<9; n++)
  {
    J[n]   = xr[n]*ys[n]-xs[n]*yr[n];
    Jrx[n] = ys[n];
    Jry[n] =-xs[n];
    Jsx[n] =-yr[n];
    Jsy[n] = xr[n];
  }

  for(int d = 1;d < 3;d++)
  {
    bfam_subdomain_dgx_t *d2 =
      bfam_subdomain_dgx_new_2(0, "2d", N, nV, d+1, Vi, K, EToV, EToE, EToF, 2);
    for(int v = 0;v < d+1;v++)
    {
      char name[BFAM_BUFSIZ];
      snprintf(name,BFAM_BUFSIZ,"_grid_x%d",v);
      failures += check_field((bfam_subdomain_t*)d2, name, vex[v]);
    }
    if(d > 1)
    {
      for(int v = 0;v < d+1;v++)
      {
        for(int i = 0;i < 2;i++)
        {
          char name[BFAM_BUFSIZ];
          snprintf(name,BFAM_BUFSIZ,"_grid_x%dr%d",v,i);
          failures += check_field((bfam_subdomain_t*)d2, name, xrex[i+v*2]);
        }
      }
    }
    else
    {
      failures += check_field((bfam_subdomain_t*)d2, "_grid_J"    , J  );
      failures += check_field((bfam_subdomain_t*)d2, "_grid_Jr0x0", Jrx);
      failures += check_field((bfam_subdomain_t*)d2, "_grid_Jr0x1", Jry);
      failures += check_field((bfam_subdomain_t*)d2, "_grid_Jr1x0", Jsx);
      failures += check_field((bfam_subdomain_t*)d2, "_grid_Jr1x1", Jsy);
    }

    bfam_subdomain_free((bfam_subdomain_t*)d2);
    bfam_free(d2);
  }
  return failures;
}

int
test_3d()
{
  int failures = 0;
  const bfam_locidx_t   EToV[] = {0,1,2,3,4,5,6,7};
  const bfam_locidx_t   EToE[] = {0,0,0,0,0,0,0,0};
  const int8_t          EToF[] = {0,1,2,3,4,5,6,7};

  const bfam_long_real_t  Vx[] = {-1, 1,-2, 3, 1, 4,-4,-1};
  const bfam_long_real_t  Vy[] = {-1,-2, 2, 5,-1, 2, 5, 6};
  const bfam_long_real_t  Vz[] = { 0,-1,-2, 1, 5, 4, 6, 3};
  const bfam_long_real_t  Va[] = { 0,-5,3, 4, 8, -1, 5,10};
  const bfam_long_real_t *Vi[] = {Vx,Vy,Vz,Va};

  const bfam_locidx_t      K   =  1;
  const int               nV   = 8;
  const int                N   = 2;

  const bfam_real_t      xex[] = {     Vx[0]       ,  0.5*(Vx[0]+Vx[1]            ),     Vx[1],
                                  0.5*(Vx[0]+Vx[2]), 0.25*(Vx[0]+Vx[1]+Vx[2]+Vx[3]),0.5*(Vx[1]+Vx[3]),
                                             Vx[2] ,  0.5*(            Vx[2]+Vx[3]),           Vx[3],
                                  0.5*(Vx[0]+Vx[4]), 0.25*(Vx[0]+Vx[1]+Vx[4]+Vx[5]),0.5*(Vx[1]+Vx[5]),
                                  0.25*(+Vx[0]+Vx[4]+Vx[2]+Vx[6]), 0.125*(Vx[0]+Vx[1]+Vx[4]+Vx[5]+Vx[2]+Vx[3]+Vx[6]+Vx[7]),0.25*(Vx[1]+Vx[5]+Vx[3]+Vx[7]),
                                  0.5*(Vx[2]+Vx[6]), 0.25*(Vx[2]+Vx[3]+Vx[6]+Vx[7]),0.5*(Vx[3]+Vx[7]),
                                       Vx[4]       ,  0.5*(Vx[4]+Vx[5]            ),     Vx[5],
                                  0.5*(Vx[4]+Vx[6]), 0.25*(Vx[4]+Vx[5]+Vx[6]+Vx[7]),0.5*(Vx[5]+Vx[7]),
                                             Vx[6] ,  0.5*(            Vx[6]+Vx[7]),           Vx[7]  };
  const bfam_real_t      yex[] = {     Vy[0]       ,  0.5*(Vy[0]+Vy[1]            ),     Vy[1],
                                  0.5*(Vy[0]+Vy[2]), 0.25*(Vy[0]+Vy[1]+Vy[2]+Vy[3]),0.5*(Vy[1]+Vy[3]),
                                             Vy[2] ,  0.5*(            Vy[2]+Vy[3]),           Vy[3],
                                  0.5*(Vy[0]+Vy[4]), 0.25*(Vy[0]+Vy[1]+Vy[4]+Vy[5]),0.5*(Vy[1]+Vy[5]),
                                  0.25*(+Vy[0]+Vy[4]+Vy[2]+Vy[6]), 0.125*(Vy[0]+Vy[1]+Vy[4]+Vy[5]+Vy[2]+Vy[3]+Vy[6]+Vy[7]),0.25*(Vy[1]+Vy[5]+Vy[3]+Vy[7]),
                                  0.5*(Vy[2]+Vy[6]), 0.25*(Vy[2]+Vy[3]+Vy[6]+Vy[7]),0.5*(Vy[3]+Vy[7]),
                                       Vy[4]       ,  0.5*(Vy[4]+Vy[5]            ),     Vy[5],
                                  0.5*(Vy[4]+Vy[6]), 0.25*(Vy[4]+Vy[5]+Vy[6]+Vy[7]),0.5*(Vy[5]+Vy[7]),
                                             Vy[6] ,  0.5*(            Vy[6]+Vy[7]),           Vy[7]  };
  const bfam_real_t      zex[] = {     Vz[0]       ,  0.5*(Vz[0]+Vz[1]            ),     Vz[1],
                                  0.5*(Vz[0]+Vz[2]), 0.25*(Vz[0]+Vz[1]+Vz[2]+Vz[3]),0.5*(Vz[1]+Vz[3]),
                                             Vz[2] ,  0.5*(            Vz[2]+Vz[3]),           Vz[3],
                                  0.5*(Vz[0]+Vz[4]), 0.25*(Vz[0]+Vz[1]+Vz[4]+Vz[5]),0.5*(Vz[1]+Vz[5]),
                                  0.25*(+Vz[0]+Vz[4]+Vz[2]+Vz[6]), 0.125*(Vz[0]+Vz[1]+Vz[4]+Vz[5]+Vz[2]+Vz[3]+Vz[6]+Vz[7]),0.25*(Vz[1]+Vz[5]+Vz[3]+Vz[7]),
                                  0.5*(Vz[2]+Vz[6]), 0.25*(Vz[2]+Vz[3]+Vz[6]+Vz[7]),0.5*(Vz[3]+Vz[7]),
                                       Vz[4]       ,  0.5*(Vz[4]+Vz[5]            ),     Vz[5],
                                  0.5*(Vz[4]+Vz[6]), 0.25*(Vz[4]+Vz[5]+Vz[6]+Vz[7]),0.5*(Vz[5]+Vz[7]),
                                             Vz[6] ,  0.5*(            Vz[6]+Vz[7]),           Vz[7]  };
  const bfam_real_t      aex[] = {     Va[0]       ,  0.5*(Va[0]+Va[1]            ),     Va[1],
                                  0.5*(Va[0]+Va[2]), 0.25*(Va[0]+Va[1]+Va[2]+Va[3]),0.5*(Va[1]+Va[3]),
                                             Va[2] ,  0.5*(            Va[2]+Va[3]),           Va[3],
                                  0.5*(Va[0]+Va[4]), 0.25*(Va[0]+Va[1]+Va[4]+Va[5]),0.5*(Va[1]+Va[5]),
                                  0.25*(+Va[0]+Va[4]+Va[2]+Va[6]), 0.125*(Va[0]+Va[1]+Va[4]+Va[5]+Va[2]+Va[3]+Va[6]+Va[7]),0.25*(Va[1]+Va[5]+Va[3]+Va[7]),
                                  0.5*(Va[2]+Va[6]), 0.25*(Va[2]+Va[3]+Va[6]+Va[7]),0.5*(Va[3]+Va[7]),
                                       Va[4]       ,  0.5*(Va[4]+Va[5]            ),     Va[5],
                                  0.5*(Va[4]+Va[6]), 0.25*(Va[4]+Va[5]+Va[6]+Va[7]),0.5*(Va[5]+Va[7]),
                                             Va[6] ,  0.5*(            Va[6]+Va[7]),           Va[7]  };
  const bfam_real_t     *vex[] = {xex,yex,zex,aex};

  const bfam_real_t       xr[] = {0.50*(xex[ 2]-xex[ 0]),0.50*(xex[ 2]-xex[ 0]),0.50*(xex[ 2]-xex[ 0]),
                                  0.50*(xex[ 5]-xex[ 3]),0.50*(xex[ 5]-xex[ 3]),0.50*(xex[ 5]-xex[ 3]),
                                  0.50*(xex[ 8]-xex[ 6]),0.50*(xex[ 8]-xex[ 6]),0.50*(xex[ 8]-xex[ 6]),
                                  0.50*(xex[11]-xex[ 9]),0.50*(xex[11]-xex[ 9]),0.50*(xex[11]-xex[ 9]),
                                  0.50*(xex[14]-xex[12]),0.50*(xex[14]-xex[12]),0.50*(xex[14]-xex[12]),
                                  0.50*(xex[17]-xex[15]),0.50*(xex[17]-xex[15]),0.50*(xex[17]-xex[15]),
                                  0.50*(xex[20]-xex[18]),0.50*(xex[20]-xex[18]),0.50*(xex[20]-xex[18]),
                                  0.50*(xex[23]-xex[21]),0.50*(xex[23]-xex[21]),0.50*(xex[23]-xex[21]),
                                  0.50*(xex[26]-xex[24]),0.50*(xex[26]-xex[24]),0.50*(xex[26]-xex[24])};
  const bfam_real_t       xs[] = {0.50*(xex[ 6]-xex[ 0]),0.50*(xex[ 7]-xex[ 1]),0.50*(xex[ 8]-xex[ 2]),
                                  0.50*(xex[ 6]-xex[ 0]),0.50*(xex[ 7]-xex[ 1]),0.50*(xex[ 8]-xex[ 2]),
                                  0.50*(xex[ 6]-xex[ 0]),0.50*(xex[ 7]-xex[ 1]),0.50*(xex[ 8]-xex[ 2]),
                                  0.50*(xex[15]-xex[ 9]),0.50*(xex[16]-xex[10]),0.50*(xex[17]-xex[11]),
                                  0.50*(xex[15]-xex[ 9]),0.50*(xex[16]-xex[10]),0.50*(xex[17]-xex[11]),
                                  0.50*(xex[15]-xex[ 9]),0.50*(xex[16]-xex[10]),0.50*(xex[17]-xex[11]),
                                  0.50*(xex[24]-xex[18]),0.50*(xex[25]-xex[19]),0.50*(xex[26]-xex[20]),
                                  0.50*(xex[24]-xex[18]),0.50*(xex[25]-xex[19]),0.50*(xex[26]-xex[20]),
                                  0.50*(xex[24]-xex[18]),0.50*(xex[25]-xex[19]),0.50*(xex[26]-xex[20])};
  const bfam_real_t       xt[] = {0.50*(xex[18]-xex[ 0]),0.50*(xex[19]-xex[ 1]),0.50*(xex[20]-xex[ 2]),
                                  0.50*(xex[21]-xex[ 3]),0.50*(xex[22]-xex[ 4]),0.50*(xex[23]-xex[ 5]),
                                  0.50*(xex[24]-xex[ 6]),0.50*(xex[25]-xex[ 7]),0.50*(xex[26]-xex[ 8]),
                                  0.50*(xex[18]-xex[ 0]),0.50*(xex[19]-xex[ 1]),0.50*(xex[20]-xex[ 2]),
                                  0.50*(xex[21]-xex[ 3]),0.50*(xex[22]-xex[ 4]),0.50*(xex[23]-xex[ 5]),
                                  0.50*(xex[24]-xex[ 6]),0.50*(xex[25]-xex[ 7]),0.50*(xex[26]-xex[ 8]),
                                  0.50*(xex[18]-xex[ 0]),0.50*(xex[19]-xex[ 1]),0.50*(xex[20]-xex[ 2]),
                                  0.50*(xex[21]-xex[ 3]),0.50*(xex[22]-xex[ 4]),0.50*(xex[23]-xex[ 5]),
                                  0.50*(xex[24]-xex[ 6]),0.50*(xex[25]-xex[ 7]),0.50*(xex[26]-xex[ 8])};

  const bfam_real_t       yr[] = {0.50*(yex[ 2]-yex[ 0]),0.50*(yex[ 2]-yex[ 0]),0.50*(yex[ 2]-yex[ 0]),
                                  0.50*(yex[ 5]-yex[ 3]),0.50*(yex[ 5]-yex[ 3]),0.50*(yex[ 5]-yex[ 3]),
                                  0.50*(yex[ 8]-yex[ 6]),0.50*(yex[ 8]-yex[ 6]),0.50*(yex[ 8]-yex[ 6]),
                                  0.50*(yex[11]-yex[ 9]),0.50*(yex[11]-yex[ 9]),0.50*(yex[11]-yex[ 9]),
                                  0.50*(yex[14]-yex[12]),0.50*(yex[14]-yex[12]),0.50*(yex[14]-yex[12]),
                                  0.50*(yex[17]-yex[15]),0.50*(yex[17]-yex[15]),0.50*(yex[17]-yex[15]),
                                  0.50*(yex[20]-yex[18]),0.50*(yex[20]-yex[18]),0.50*(yex[20]-yex[18]),
                                  0.50*(yex[23]-yex[21]),0.50*(yex[23]-yex[21]),0.50*(yex[23]-yex[21]),
                                  0.50*(yex[26]-yex[24]),0.50*(yex[26]-yex[24]),0.50*(yex[26]-yex[24])};
  const bfam_real_t       ys[] = {0.50*(yex[ 6]-yex[ 0]),0.50*(yex[ 7]-yex[ 1]),0.50*(yex[ 8]-yex[ 2]),
                                  0.50*(yex[ 6]-yex[ 0]),0.50*(yex[ 7]-yex[ 1]),0.50*(yex[ 8]-yex[ 2]),
                                  0.50*(yex[ 6]-yex[ 0]),0.50*(yex[ 7]-yex[ 1]),0.50*(yex[ 8]-yex[ 2]),
                                  0.50*(yex[15]-yex[ 9]),0.50*(yex[16]-yex[10]),0.50*(yex[17]-yex[11]),
                                  0.50*(yex[15]-yex[ 9]),0.50*(yex[16]-yex[10]),0.50*(yex[17]-yex[11]),
                                  0.50*(yex[15]-yex[ 9]),0.50*(yex[16]-yex[10]),0.50*(yex[17]-yex[11]),
                                  0.50*(yex[24]-yex[18]),0.50*(yex[25]-yex[19]),0.50*(yex[26]-yex[20]),
                                  0.50*(yex[24]-yex[18]),0.50*(yex[25]-yex[19]),0.50*(yex[26]-yex[20]),
                                  0.50*(yex[24]-yex[18]),0.50*(yex[25]-yex[19]),0.50*(yex[26]-yex[20])};
  const bfam_real_t       yt[] = {0.50*(yex[18]-yex[ 0]),0.50*(yex[19]-yex[ 1]),0.50*(yex[20]-yex[ 2]),
                                  0.50*(yex[21]-yex[ 3]),0.50*(yex[22]-yex[ 4]),0.50*(yex[23]-yex[ 5]),
                                  0.50*(yex[24]-yex[ 6]),0.50*(yex[25]-yex[ 7]),0.50*(yex[26]-yex[ 8]),
                                  0.50*(yex[18]-yex[ 0]),0.50*(yex[19]-yex[ 1]),0.50*(yex[20]-yex[ 2]),
                                  0.50*(yex[21]-yex[ 3]),0.50*(yex[22]-yex[ 4]),0.50*(yex[23]-yex[ 5]),
                                  0.50*(yex[24]-yex[ 6]),0.50*(yex[25]-yex[ 7]),0.50*(yex[26]-yex[ 8]),
                                  0.50*(yex[18]-yex[ 0]),0.50*(yex[19]-yex[ 1]),0.50*(yex[20]-yex[ 2]),
                                  0.50*(yex[21]-yex[ 3]),0.50*(yex[22]-yex[ 4]),0.50*(yex[23]-yex[ 5]),
                                  0.50*(yex[24]-yex[ 6]),0.50*(yex[25]-yex[ 7]),0.50*(yex[26]-yex[ 8])};

  const bfam_real_t       zr[] = {0.50*(zex[ 2]-zex[ 0]),0.50*(zex[ 2]-zex[ 0]),0.50*(zex[ 2]-zex[ 0]),
                                  0.50*(zex[ 5]-zex[ 3]),0.50*(zex[ 5]-zex[ 3]),0.50*(zex[ 5]-zex[ 3]),
                                  0.50*(zex[ 8]-zex[ 6]),0.50*(zex[ 8]-zex[ 6]),0.50*(zex[ 8]-zex[ 6]),
                                  0.50*(zex[11]-zex[ 9]),0.50*(zex[11]-zex[ 9]),0.50*(zex[11]-zex[ 9]),
                                  0.50*(zex[14]-zex[12]),0.50*(zex[14]-zex[12]),0.50*(zex[14]-zex[12]),
                                  0.50*(zex[17]-zex[15]),0.50*(zex[17]-zex[15]),0.50*(zex[17]-zex[15]),
                                  0.50*(zex[20]-zex[18]),0.50*(zex[20]-zex[18]),0.50*(zex[20]-zex[18]),
                                  0.50*(zex[23]-zex[21]),0.50*(zex[23]-zex[21]),0.50*(zex[23]-zex[21]),
                                  0.50*(zex[26]-zex[24]),0.50*(zex[26]-zex[24]),0.50*(zex[26]-zex[24])};
  const bfam_real_t       zs[] = {0.50*(zex[ 6]-zex[ 0]),0.50*(zex[ 7]-zex[ 1]),0.50*(zex[ 8]-zex[ 2]),
                                  0.50*(zex[ 6]-zex[ 0]),0.50*(zex[ 7]-zex[ 1]),0.50*(zex[ 8]-zex[ 2]),
                                  0.50*(zex[ 6]-zex[ 0]),0.50*(zex[ 7]-zex[ 1]),0.50*(zex[ 8]-zex[ 2]),
                                  0.50*(zex[15]-zex[ 9]),0.50*(zex[16]-zex[10]),0.50*(zex[17]-zex[11]),
                                  0.50*(zex[15]-zex[ 9]),0.50*(zex[16]-zex[10]),0.50*(zex[17]-zex[11]),
                                  0.50*(zex[15]-zex[ 9]),0.50*(zex[16]-zex[10]),0.50*(zex[17]-zex[11]),
                                  0.50*(zex[24]-zex[18]),0.50*(zex[25]-zex[19]),0.50*(zex[26]-zex[20]),
                                  0.50*(zex[24]-zex[18]),0.50*(zex[25]-zex[19]),0.50*(zex[26]-zex[20]),
                                  0.50*(zex[24]-zex[18]),0.50*(zex[25]-zex[19]),0.50*(zex[26]-zex[20])};
  const bfam_real_t       zt[] = {0.50*(zex[18]-zex[ 0]),0.50*(zex[19]-zex[ 1]),0.50*(zex[20]-zex[ 2]),
                                  0.50*(zex[21]-zex[ 3]),0.50*(zex[22]-zex[ 4]),0.50*(zex[23]-zex[ 5]),
                                  0.50*(zex[24]-zex[ 6]),0.50*(zex[25]-zex[ 7]),0.50*(zex[26]-zex[ 8]),
                                  0.50*(zex[18]-zex[ 0]),0.50*(zex[19]-zex[ 1]),0.50*(zex[20]-zex[ 2]),
                                  0.50*(zex[21]-zex[ 3]),0.50*(zex[22]-zex[ 4]),0.50*(zex[23]-zex[ 5]),
                                  0.50*(zex[24]-zex[ 6]),0.50*(zex[25]-zex[ 7]),0.50*(zex[26]-zex[ 8]),
                                  0.50*(zex[18]-zex[ 0]),0.50*(zex[19]-zex[ 1]),0.50*(zex[20]-zex[ 2]),
                                  0.50*(zex[21]-zex[ 3]),0.50*(zex[22]-zex[ 4]),0.50*(zex[23]-zex[ 5]),
                                  0.50*(zex[24]-zex[ 6]),0.50*(zex[25]-zex[ 7]),0.50*(zex[26]-zex[ 8])};

  const bfam_real_t       ar[] = {0.50*(aex[ 2]-aex[ 0]),0.50*(aex[ 2]-aex[ 0]),0.50*(aex[ 2]-aex[ 0]),
                                  0.50*(aex[ 5]-aex[ 3]),0.50*(aex[ 5]-aex[ 3]),0.50*(aex[ 5]-aex[ 3]),
                                  0.50*(aex[ 8]-aex[ 6]),0.50*(aex[ 8]-aex[ 6]),0.50*(aex[ 8]-aex[ 6]),
                                  0.50*(aex[11]-aex[ 9]),0.50*(aex[11]-aex[ 9]),0.50*(aex[11]-aex[ 9]),
                                  0.50*(aex[14]-aex[12]),0.50*(aex[14]-aex[12]),0.50*(aex[14]-aex[12]),
                                  0.50*(aex[17]-aex[15]),0.50*(aex[17]-aex[15]),0.50*(aex[17]-aex[15]),
                                  0.50*(aex[20]-aex[18]),0.50*(aex[20]-aex[18]),0.50*(aex[20]-aex[18]),
                                  0.50*(aex[23]-aex[21]),0.50*(aex[23]-aex[21]),0.50*(aex[23]-aex[21]),
                                  0.50*(aex[26]-aex[24]),0.50*(aex[26]-aex[24]),0.50*(aex[26]-aex[24])};
  const bfam_real_t       as[] = {0.50*(aex[ 6]-aex[ 0]),0.50*(aex[ 7]-aex[ 1]),0.50*(aex[ 8]-aex[ 2]),
                                  0.50*(aex[ 6]-aex[ 0]),0.50*(aex[ 7]-aex[ 1]),0.50*(aex[ 8]-aex[ 2]),
                                  0.50*(aex[ 6]-aex[ 0]),0.50*(aex[ 7]-aex[ 1]),0.50*(aex[ 8]-aex[ 2]),
                                  0.50*(aex[15]-aex[ 9]),0.50*(aex[16]-aex[10]),0.50*(aex[17]-aex[11]),
                                  0.50*(aex[15]-aex[ 9]),0.50*(aex[16]-aex[10]),0.50*(aex[17]-aex[11]),
                                  0.50*(aex[15]-aex[ 9]),0.50*(aex[16]-aex[10]),0.50*(aex[17]-aex[11]),
                                  0.50*(aex[24]-aex[18]),0.50*(aex[25]-aex[19]),0.50*(aex[26]-aex[20]),
                                  0.50*(aex[24]-aex[18]),0.50*(aex[25]-aex[19]),0.50*(aex[26]-aex[20]),
                                  0.50*(aex[24]-aex[18]),0.50*(aex[25]-aex[19]),0.50*(aex[26]-aex[20])};
  const bfam_real_t       at[] = {0.50*(aex[18]-aex[ 0]),0.50*(aex[19]-aex[ 1]),0.50*(aex[20]-aex[ 2]),
                                  0.50*(aex[21]-aex[ 3]),0.50*(aex[22]-aex[ 4]),0.50*(aex[23]-aex[ 5]),
                                  0.50*(aex[24]-aex[ 6]),0.50*(aex[25]-aex[ 7]),0.50*(aex[26]-aex[ 8]),
                                  0.50*(aex[18]-aex[ 0]),0.50*(aex[19]-aex[ 1]),0.50*(aex[20]-aex[ 2]),
                                  0.50*(aex[21]-aex[ 3]),0.50*(aex[22]-aex[ 4]),0.50*(aex[23]-aex[ 5]),
                                  0.50*(aex[24]-aex[ 6]),0.50*(aex[25]-aex[ 7]),0.50*(aex[26]-aex[ 8]),
                                  0.50*(aex[18]-aex[ 0]),0.50*(aex[19]-aex[ 1]),0.50*(aex[20]-aex[ 2]),
                                  0.50*(aex[21]-aex[ 3]),0.50*(aex[22]-aex[ 4]),0.50*(aex[23]-aex[ 5]),
                                  0.50*(aex[24]-aex[ 6]),0.50*(aex[25]-aex[ 7]),0.50*(aex[26]-aex[ 8])};

  const bfam_real_t    *xrex[] = {xr,xs,xt,yr,ys,yt,zr,zs,zt,ar,as,at};

  bfam_real_t J[27],
              Jrx[27],Jry[27],Jrz[27],
              Jsx[27],Jsy[27],Jsz[27],
              Jtx[27],Jty[27],Jtz[27];
  for(int n = 0; n<27; n++)
  {
    J[n]   = xr[n]*(ys[n]*zt[n]-yt[n]*zs[n])
            -xs[n]*(yr[n]*zt[n]-yt[n]*zr[n])
            +xt[n]*(yr[n]*zs[n]-ys[n]*zr[n]);
    Jrx[n] =  ys[n]*zt[n]-yt[n]*zs[n];
    Jry[n] =  zs[n]*xt[n]-zt[n]*xs[n];
    Jrz[n] =  xs[n]*yt[n]-xt[n]*ys[n];
    Jsx[n] =  yt[n]*zr[n]-yr[n]*zt[n];
    Jsy[n] =  zt[n]*xr[n]-zr[n]*xt[n];
    Jsz[n] =  xt[n]*yr[n]-xr[n]*yt[n];
    Jtx[n] =  yr[n]*zs[n]-ys[n]*zr[n];
    Jty[n] =  zr[n]*xs[n]-zs[n]*xr[n];
    Jtz[n] =  xr[n]*ys[n]-xs[n]*yr[n];
  }

  for(int d = 2;d < 4; d++)
  {
    bfam_subdomain_dgx_t *d3 =
      bfam_subdomain_dgx_new_3(0, "3d", N, nV, d+1, Vi, K, EToV, EToE, EToF, 3);
    for(int v = 0;v < d+1;v++)
    {
      char name[BFAM_BUFSIZ];
      snprintf(name,BFAM_BUFSIZ,"_grid_x%d",v);
      failures += check_field((bfam_subdomain_t*)d3, name, vex[v]);
    }
    if(d > 2)
    {
      for(int v = 0;v < d+1;v++)
      {
        for(int i = 0;i < 3;i++)
        {
          char name[BFAM_BUFSIZ];
          snprintf(name,BFAM_BUFSIZ,"_grid_x%dr%d",v,i);
          failures += check_field((bfam_subdomain_t*)d3, name, xrex[i+v*3]);
        }
      }
    }
    else
    {
      failures += check_field((bfam_subdomain_t*)d3, "_grid_J"    , J  );
      failures += check_field((bfam_subdomain_t*)d3, "_grid_Jr0x0", Jrx);
      failures += check_field((bfam_subdomain_t*)d3, "_grid_Jr0x1", Jry);
      failures += check_field((bfam_subdomain_t*)d3, "_grid_Jr0x2", Jrz);
      failures += check_field((bfam_subdomain_t*)d3, "_grid_Jr1x0", Jsx);
      failures += check_field((bfam_subdomain_t*)d3, "_grid_Jr1x1", Jsy);
      failures += check_field((bfam_subdomain_t*)d3, "_grid_Jr1x2", Jsz);
      failures += check_field((bfam_subdomain_t*)d3, "_grid_Jr2x0", Jtx);
      failures += check_field((bfam_subdomain_t*)d3, "_grid_Jr2x1", Jty);
      failures += check_field((bfam_subdomain_t*)d3, "_grid_Jr2x2", Jtz);
    }

    bfam_subdomain_free((bfam_subdomain_t*)d3);
    bfam_free(d3);
  }

  return failures;
}

int
main (int argc, char *argv[])
{
  int rank;
  BFAM_MPI_CHECK(MPI_Init(&argc,&argv));
  BFAM_MPI_CHECK(MPI_Comm_rank(MPI_COMM_WORLD, &rank));

  int check = 0;

  check += test_1d();
  check += test_2d();
  check += test_3d();

  BFAM_MPI_CHECK(MPI_Finalize());
  return check;
}

