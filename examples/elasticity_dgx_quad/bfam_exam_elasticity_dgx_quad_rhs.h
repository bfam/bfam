#ifndef BFAM_EXAM_ELASTICITY_DGX_QUAD_RHS_H
#define BFAM_EXAM_ELASTICITY_DGX_QUAD_RHS_H

#include <bfam.h>

#define X(order) \
  void bfam_elasticity_dgx_quad_print_order_##order(int N);      \
  void bfam_elasticity_dgx_quad_intra_rhs_##order(int N,         \
      bfam_subdomain_dgx_quad_t *sub, const char *rate_prefix,   \
      const char *field_prefix);
BFAM_LIST_OF_DGX_QUAD_NORDERS
#undef X

void bfam_elasticity_dgx_quad_print_order_(int N);
void bfam_elasticity_dgx_quad_intra_rhs_(int N, bfam_subdomain_dgx_quad_t *sub,
    const char *rate_prefix, const char *field_prefix);

#endif
