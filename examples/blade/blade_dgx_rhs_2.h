#ifndef BLADE_DGX_RHS_2_H
#define BLADE_DGX_RHS_2_H

#include <bfam.h>

#define X(order)                                                     \
  void blade_dgx_energy_2_##order(int N,                             \
      bfam_real_t* energy_sq, bfam_subdomain_dgx_t *sub,             \
      const char *field_prefix);
BFAM_LIST_OF_DGX_NORDERS
#undef X

void blade_dgx_energy_2_(int N, bfam_real_t *energy_sq,
    bfam_subdomain_dgx_t *sub, const char *field_prefix);

#endif
