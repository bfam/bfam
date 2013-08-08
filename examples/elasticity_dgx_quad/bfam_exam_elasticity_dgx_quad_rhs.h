#ifndef BFAM_EXAM_ELASTICITY_DGX_QUAD_RHS_H
#define BFAM_EXAM_ELASTICITY_DGX_QUAD_RHS_H

#include <bfam.h>

#define X(order) \
  void bfam_elasticity_dgx_quad_print_order_##order(int N);        \
  void bfam_elasticity_dgx_quad_intra_rhs_elastic_##order(int N,   \
      bfam_subdomain_dgx_quad_t *sub, const char *rate_prefix,     \
      const char *field_prefix, const bfam_long_real_t t);         \
  void bfam_elasticity_dgx_quad_scale_rates_elastic_##order(int N, \
      bfam_subdomain_dgx_quad_t *sub, const char *rate_prefix,     \
      const bfam_long_real_t a);                                   \
  void bfam_elasticity_dgx_quad_add_rates_elastic_##order(int N,     \
      bfam_subdomain_dgx_quad_t *sub, const char *field_prefix_lhs,\
      const char *field_prefix_rhs, const char *rate_prefix,       \
      const bfam_long_real_t a);
BFAM_LIST_OF_DGX_QUAD_NORDERS
#undef X

void bfam_elasticity_dgx_quad_print_order_(int N);
void bfam_elasticity_dgx_quad_intra_rhs_elastic_(int N,
    bfam_subdomain_dgx_quad_t *sub, const char *rate_prefix,
    const char *field_prefix, const bfam_long_real_t t);
void bfam_elasticity_dgx_quad_scale_rates_elastic_(int N,
    bfam_subdomain_dgx_quad_t *sub, const char *rate_prefix,
    const bfam_long_real_t a);
void bfam_elasticity_dgx_quad_add_rates_elastic_(int N,
    bfam_subdomain_dgx_quad_t *sub, const char *field_prefix_lhs,
    const char *field_prefix_rhs, const char *rate_prefix,
    const bfam_long_real_t a);

#endif
