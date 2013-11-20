#include <bfam.h>

int
main (int argc, char *argv[])
{
  const bfam_long_real_t  Vx[] = {0,1,2,0,1,2,0,1,2,0,1,2};
  const bfam_long_real_t  Vy[] = {0,0,0,1,1,1,0,0,0,1,1,1};
  const bfam_long_real_t  Vz[] = {0,0,0,0,0,0,1,1,1,1,1,1};
  const bfam_long_real_t *Vi[] = {Vx,Vy,Vz};

  bfam_subdomain_dgx_t *d0 =
    bfam_subdomain_dgx_new_0(1, "dim0", 0, 1, 1, NULL, 1, NULL, NULL, NULL, 0);

  bfam_subdomain_dgx_t *d1_1d =
    bfam_subdomain_dgx_new_1(2, "dim1", 8,  3, 1, Vi, 1, NULL, NULL, NULL, 1);

  bfam_subdomain_dgx_t *d1_2d =
    bfam_subdomain_dgx_new_1(2, "dim1", 8,  3, 2, Vi, 1, NULL, NULL, NULL, 1);

  bfam_subdomain_dgx_t *d1_3d =
    bfam_subdomain_dgx_new_1(2, "dim1", 8,  3, 3, Vi, 1, NULL, NULL, NULL, 1);

  bfam_subdomain_dgx_t *d2_2d =
    bfam_subdomain_dgx_new_2(3, "dim2", 6,  6, 2, Vi, 1, NULL, NULL, NULL, 2);

  bfam_subdomain_dgx_t *d2_3d =
    bfam_subdomain_dgx_new_2(3, "dim2", 6,  6, 3, Vi, 1, NULL, NULL, NULL, 2);

  bfam_subdomain_dgx_t *d3 =
    bfam_subdomain_dgx_new_3(5, "dim3", 4, 12, 3, Vi, 1, NULL, NULL, NULL, 3);

  d0->base.free((bfam_subdomain_t*)d0); bfam_free(d0);

  d1_1d->base.free((bfam_subdomain_t*)d1_1d); bfam_free(d1_1d);
  d1_2d->base.free((bfam_subdomain_t*)d1_2d); bfam_free(d1_2d);
  d1_3d->base.free((bfam_subdomain_t*)d1_3d); bfam_free(d1_3d);

  d2_2d->base.free((bfam_subdomain_t*)d2_2d); bfam_free(d2_2d);
  d2_3d->base.free((bfam_subdomain_t*)d2_3d); bfam_free(d2_3d);

  d3->base.free((bfam_subdomain_t*)d3); bfam_free(d3);

  return 0;
}

