#ifndef BFAM_EXAM_ELASTICITY_DGX_QUAD_RHS_H
#define BFAM_EXAM_ELASTICITY_DGX_QUAD_RHS_H

#include <bfam.h>

#define X(order) \
  void bfam_elasticity_dgx_quad_print_order_##order(int N);

BFAM_LIST_OF_DGX_QUAD_NORDERS
#undef X

void bfam_elasticity_dgx_quad_print_order_(int N);

#endif
