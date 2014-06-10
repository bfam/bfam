#ifndef BEARD_DGX_RHS_2_H
#define BEARD_DGX_RHS_2_H

#include <bfam.h>

#define X(order)                                                     \
  void beard_dgx_intra_rhs_elastic_2_##order(int N,                  \
      bfam_subdomain_dgx_t *sub, const char *rate_prefix,            \
      const char *field_prefix, const bfam_long_real_t t);           \
  void beard_dgx_scale_rates_elastic_2_##order(int N,                \
      bfam_subdomain_dgx_t *sub, const char *rate_prefix,            \
      const bfam_long_real_t a);                                     \
  void beard_dgx_add_rates_elastic_2_##order(int N,                  \
      bfam_subdomain_dgx_t *sub, const char *field_prefix_lhs,       \
      const char *field_prefix_rhs, const char *rate_prefix,         \
      const bfam_long_real_t a);                                     \
  void beard_dgx_scale_rates_slip_weakening_2_##order(int N,         \
      bfam_subdomain_dgx_t *sub, const char *rate_prefix,            \
      const bfam_long_real_t a);                                     \
  void beard_dgx_add_rates_slip_weakening_2_##order(int N,           \
      bfam_subdomain_dgx_t *sub, const char *field_prefix_lhs,       \
      const char *field_prefix_rhs, const char *rate_prefix,         \
      const bfam_long_real_t a);                                     \
  void beard_dgx_add_rates_glue_m_##order(int N,                     \
      bfam_subdomain_dgx_t *sub, const char *field_prefix_lhs,       \
      const char *field_prefix_rhs, const char *rate_prefix,         \
      const bfam_long_real_t a,                                      \
      const char** scals, const char** vects, const char** tens);    \
  void beard_dgx_inter_rhs_boundary_2_##order(int N,                 \
      bfam_subdomain_dgx_t *sub, const char *rate_prefix,            \
      const char *field_prefix, const bfam_long_real_t t,            \
      const bfam_real_t R);                                          \
  void beard_dgx_inter_rhs_interface_2_##order(int N,                \
      bfam_subdomain_dgx_t *sub, const char *rate_prefix,            \
      const char *field_prefix, const bfam_long_real_t t);           \
  void beard_dgx_inter_rhs_slip_weakening_interface_2_##order(int N, \
      bfam_subdomain_dgx_t *sub, const char *rate_prefix,            \
      const char *minus_rate_prefix, const char *field_prefix,       \
      const bfam_long_real_t t);                                     \
  void beard_dgx_energy_2_##order(int N,                             \
      bfam_real_t* energy_sq, bfam_subdomain_dgx_t *sub,             \
      const char *field_prefix);
BFAM_LIST_OF_DGX_NORDERS
#undef X

void beard_dgx_intra_rhs_elastic_2_(int N,
    bfam_subdomain_dgx_t *sub, const char *rate_prefix,
    const char *field_prefix, const bfam_long_real_t t);
void beard_dgx_scale_rates_elastic_2_(int N,
    bfam_subdomain_dgx_t *sub, const char *rate_prefix,
    const bfam_long_real_t a);
void beard_dgx_add_rates_elastic_2_(int N,
    bfam_subdomain_dgx_t *sub, const char *field_prefix_lhs,
    const char *field_prefix_rhs, const char *rate_prefix,
    const bfam_long_real_t a);
void beard_dgx_scale_rates_slip_weakening_2_(int N,
    bfam_subdomain_dgx_t *sub, const char *rate_prefix,
    const bfam_long_real_t a);
void beard_dgx_add_rates_slip_weakening_2_(int N,
    bfam_subdomain_dgx_t *sub, const char *field_prefix_lhs,
    const char *field_prefix_rhs, const char *rate_prefix,
    const bfam_long_real_t a);
void beard_dgx_add_rates_glue_m_2(int N,
    bfam_subdomain_dgx_t *sub, const char *field_prefix_lhs,
    const char *field_prefix_rhs, const char *rate_prefix,
    const bfam_long_real_t a,
    const char** scals, const char** vects, const char** tens);
void beard_dgx_inter_rhs_boundary_2_(int N,
    bfam_subdomain_dgx_t *sub, const char *rate_prefix,
    const char *field_prefix, const bfam_long_real_t t, const bfam_real_t R);
void beard_dgx_inter_rhs_interface_2_(int N,
    bfam_subdomain_dgx_t *sub, const char *rate_prefix,
    const char *field_prefix, const bfam_long_real_t t);
void beard_dgx_inter_rhs_slip_weakening_interface_2_(int N,
    bfam_subdomain_dgx_t *sub, const char *rate_prefix,
    const char *minus_rate_prefix, const char *field_prefix,
    const bfam_long_real_t t);
void beard_dgx_energy_2_(int N, bfam_real_t *energy_sq,
    bfam_subdomain_dgx_t *sub, const char *field_prefix);

#endif
