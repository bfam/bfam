#include <bfam.h>

int
main (int argc, char *argv[])
{

  bfam_subdomain_dgx_t *d0 =
    bfam_subdomain_dgx_new_0(1, "dim0", 1, 1, 1, NULL, 1, NULL, NULL, NULL, 0);
  bfam_subdomain_dgx_t *d1 =
    bfam_subdomain_dgx_new_1(2, "dim1", 1, 1, 1, NULL, 1, NULL, NULL, NULL, 1);
  bfam_subdomain_dgx_t *d2 =
    bfam_subdomain_dgx_new_2(3, "dim2", 1, 1, 1, NULL, 1, NULL, NULL, NULL, 2);
  bfam_subdomain_dgx_t *d3 =
    bfam_subdomain_dgx_new_3(5, "dim3", 1, 1, 1, NULL, 1, NULL, NULL, NULL, 3);

  d0->base.free((bfam_subdomain_t*)d0);
  bfam_free(d0);

  d1->base.free((bfam_subdomain_t*)d1);
  bfam_free(d1);

  d2->base.free((bfam_subdomain_t*)d2);
  bfam_free(d2);

  d3->base.free((bfam_subdomain_t*)d3);
  bfam_free(d3);

  return 0;
}

