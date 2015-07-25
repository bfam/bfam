#include <bfam.h>

int test_0d()
{
  bfam_subdomain_dgx_t *d0 = bfam_subdomain_dgx_new_0(
      1, -1, "dim0", 0, 1, 1, NULL, 1, NULL, NULL, NULL, NULL, NULL, NULL, 0);

  d0->base.free((bfam_subdomain_t *)d0);
  bfam_free(d0);
  return 0;
}

int test_1d()
{
  const bfam_long_real_t Vx[] = {-1.2, 0, 3.4};
  const bfam_long_real_t Vy[] = {-1, 0, -1};
  const bfam_long_real_t Vz[] = {1, 2, 3};
  const bfam_long_real_t *Vi[] = {Vx, Vy, Vz};

  const bfam_locidx_t EToV[] = {0, 1, 1, 2};
  const bfam_locidx_t EToE[] = {0, 1, 0, 1};
  const int8_t EToF[] = {0, 0, 1, 1};

  for (int d = 0; d < 3; d++)
  {
    bfam_domain_t domain;
    bfam_domain_init(&domain, MPI_COMM_WORLD);

    bfam_subdomain_dgx_t *d1 = bfam_subdomain_dgx_new_1(
        0, -1, "1d", 8, 3, d + 1, Vi, 2, NULL, EToV, EToE, EToF, NULL, NULL, 1);

    bfam_domain_add_subdomain(&domain, (bfam_subdomain_t *)d1);

    const char *volume[] = {NULL};
    const char *ps[] = {"_grid_x0", NULL};
    const char *ve[] = {"g", NULL};
    const char *vc[] = {"_grid_x0", "_grid_x0", "_grid_x0", NULL};

    char name[BFAM_BUFSIZ];
    snprintf(name, BFAM_BUFSIZ, "d1_%d", d);
    bfam_vtk_write_file(&domain, BFAM_DOMAIN_AND, volume, NULL, name, 0, ps, ve,
                        vc, 1, 0, 0);

    bfam_domain_free(&domain);
  }
  return 0;
}

int test_2d()
{
  const bfam_locidx_t EToV[] = {0, 1, 3, 4, 1, 2, 4, 5, 4, 5, 3, 6};
  const bfam_locidx_t EToE[] = {0, 1, 0, 2, 0, 1, 1, 2, 0, 2, 1, 2};
  const int8_t EToF[] = {0, 0, 2, 0 + 4, 1, 1, 2, 2, 3 + 4, 1, 3, 3};
  const bfam_locidx_t K = 3;

  const bfam_long_real_t Vx[] = {0, 0.25, 1, 0, 0.25, 0.5, 0};
  const bfam_long_real_t Vy[] = {0, 0, 0, 0.5, 0.25, 0.5, 1};
  const bfam_long_real_t Vz[] = {0.1, -0.1, 0.1, 0.2, -0.2, 0.3, -0.1};
  const bfam_long_real_t *Vi[] = {Vx, Vy, Vz};
  const int nV = 7;

  const int N = 8;

  for (int d = 1; d < 3; d++)
  {

    bfam_domain_t domain;
    bfam_domain_init(&domain, MPI_COMM_WORLD);

    bfam_subdomain_dgx_t *d2 =
        bfam_subdomain_dgx_new_2(0, -1, "2d", N, nV, d + 1, Vi, K, NULL, EToV,
                                 EToE, EToF, NULL, NULL, 2);

    bfam_domain_add_subdomain(&domain, (bfam_subdomain_t *)d2);

    const char *volume[] = {NULL};
    const char *ps[] = {"_grid_x0", NULL};
    const char *ve[] = {"g", NULL};
    const char *vc[] = {"_grid_x0", "_grid_x1", "_grid_x1", NULL};

    char name[BFAM_BUFSIZ];
    snprintf(name, BFAM_BUFSIZ, "d2_%d", d);
    bfam_vtk_write_file(&domain, BFAM_DOMAIN_AND, volume, NULL, name, 0, ps, ve,
                        vc, 1, 0, 0);

    bfam_domain_free(&domain);
  }

  return 0;
}

int test_3d()
{
  /* connectivity from p8est_connectivity_new_rotcubes */
  const bfam_locidx_t K = 6;
  const bfam_locidx_t EToV[6 * 8] = {
      0,  17, 3,  4,  15, 11, 13, 14, 7,  2,  6,  17, 9,  12, 8,  11,
      2,  12, 5,  10, 17, 11, 4,  14, 19, 13, 18, 14, 16, 15, 1,  11,
      14, 11, 21, 25, 18, 1,  22, 23, 21, 20, 25, 24, 14, 10, 11, 12,
  };
  const bfam_locidx_t EToE[6 * 6] = {
      0, 2, 0, 0, 0, 3, 1, 2, 1, 1, 1, 1, 2, 5, 1, 2, 2, 0,
      3, 0, 3, 4, 3, 3, 4, 4, 3, 4, 5, 4, 4, 5, 5, 5, 5, 2,
  };
  const int8_t EToF[6 * 6] = {
      0, 5,  2, 3, 4, 13, 0, 2, 2, 3, 4,  5, 0,  23, 1, 3, 4, 1,
      0, 17, 2, 8, 4, 5,  0, 1, 9, 3, 12, 5, 16, 1,  2, 3, 4, 19,
  };

  const bfam_long_real_t Vx[26] = {0, 1, 2, 0, 1, 2, 1, 2,   1, 2, 2, 1,   2,
                                   0, 1, 0, 0, 1, 1, 0, 2.5, 2, 2, 2, 2.5, 2};
  const bfam_long_real_t Vy[26] = {0,  0, 0,   1,   1,   1,  -1, -1, -1,
                                   -1, 1, 0,   0,   1,   1,  0,  0,  0,
                                   1,  1, 1.5, 1.5, 1.5, .5, .5, .5};
  const bfam_long_real_t Vz[26] = {0, 2, 0, 0, 0, 0, 0, 0, 1, 1,   1,   1, 1,
                                   1, 1, 1, 2, 0, 2, 2, 2, 2, 2.5, 2.5, 2, 2};
  const bfam_long_real_t *Vi[] = {Vx, Vy, Vz};
  const int nV = 26;

  const int N = 4;

  bfam_domain_t domain;
  bfam_domain_init(&domain, MPI_COMM_WORLD);

  bfam_subdomain_dgx_t *d3 = bfam_subdomain_dgx_new_3(
      0, -1, "3d", N, nV, 3, Vi, K, NULL, EToV, EToE, EToF, NULL, NULL, 3);

  bfam_domain_add_subdomain(&domain, (bfam_subdomain_t *)d3);

  const char *volume[] = {NULL};
  const char *ps[] = {"_grid_x0", NULL};
  const char *ve[] = {"g", NULL};
  const char *vc[] = {"_grid_x0", "_grid_x1", "_grid_x2", NULL};
  bfam_vtk_write_file(&domain, BFAM_DOMAIN_AND, volume, NULL, "d3", 0, ps, ve,
                      vc, 1, 0, 0);

  bfam_domain_free(&domain);

  return 0;
}

int main(int argc, char *argv[])
{
  int rank;
  BFAM_MPI_CHECK(MPI_Init(&argc, &argv));
  BFAM_MPI_CHECK(MPI_Comm_rank(MPI_COMM_WORLD, &rank));

  int check = 0;

  check += test_0d();
  check += test_1d();
  check += test_2d();
  check += test_3d();

  BFAM_MPI_CHECK(MPI_Finalize());
  return check;
}
