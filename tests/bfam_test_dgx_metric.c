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
  return 0;
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

