#ifndef BEARD_DGX_RHS_3_H
#define BEARD_DGX_RHS_3_H

#include <bfam.h>

#define X(order)                                                               \
  void beard_dgx_intra_rhs_elastic_3_##order(                                  \
      int N, bfam_subdomain_dgx_t *sub, const char *rate_prefix,               \
      const char *field_prefix, const bfam_long_real_t t);                     \
  void beard_dgx_intra_rhs_sponge_3_##order(                                   \
      int N, bfam_subdomain_dgx_t *sub, const char *rate_prefix,               \
      const char *field_prefix, const bfam_long_real_t t);                     \
  void beard_dgx_rupture_time_3_##order(                                       \
      int N, bfam_subdomain_dgx_t *sub, const char *field_prefix,              \
      const bfam_long_real_t t, const bfam_real_t Vrup);                       \
  void beard_dgx_duvaut_lions_return_map_3_##order(                            \
      int N, bfam_subdomain_dgx_t *sub, const char *field_prefix,              \
      const bfam_long_real_t t, const bfam_long_real_t dt);                    \
  void beard_dgx_scale_rates_elastic_3_##order(                                \
      int N, bfam_subdomain_dgx_t *sub, const char *rate_prefix,               \
      const bfam_long_real_t a);                                               \
  void beard_dgx_scale_rates_interface_3_##order(                              \
      int N, bfam_subdomain_dgx_t *sub, const char *rate_prefix,               \
      const bfam_long_real_t a, const char **f_names);                         \
  void beard_dgx_scale_rates_3_##order(                                        \
      int N, bfam_subdomain_dgx_t *sub, const char *rate_prefix,               \
      const bfam_long_real_t a, bfam_dictionary_t *fields, const char **scals, \
      const char **vects, const char **tens, const bfam_locidx_t num_pts);     \
  void beard_dgx_add_rates_elastic_3_##order(                                  \
      int N, bfam_subdomain_dgx_t *sub, const char *field_prefix_lhs,          \
      const char *field_prefix_rhs, const char *rate_prefix,                   \
      const bfam_long_real_t a);                                               \
  void beard_dgx_add_rates_interface_3_##order(                                \
      int N, bfam_subdomain_dgx_t *sub, const char *field_prefix_lhs,          \
      const char *field_prefix_rhs, const char *rate_prefix,                   \
      const bfam_long_real_t a, const char **f_names);                         \
  void beard_dgx_add_rates_glue_p_3_##order(                                   \
      int N, bfam_subdomain_dgx_t *sub, const char *field_prefix_lhs,          \
      const char *field_prefix_rhs, const char *rate_prefix,                   \
      const bfam_long_real_t a, const char **scals, const char **vects,        \
      const char **tens);                                                      \
  void beard_dgx_add_rates_3_##order(                                          \
      int N, bfam_subdomain_dgx_t *sub, const char *field_prefix_lhs,          \
      const char *field_prefix_rhs, const char *rate_prefix,                   \
      const bfam_long_real_t a, bfam_dictionary_t *fields, const char **scals, \
      const char **vects, const char **tens, const bfam_locidx_t num_pts);     \
  void beard_dgx_inter_rhs_boundary_3_##order(                                 \
      int N, bfam_subdomain_dgx_t *sub, const char *rate_prefix,               \
      const char *field_prefix, const bfam_long_real_t t,                      \
      const bfam_real_t R);                                                    \
  void beard_dgx_inter_rhs_interface_3_##order(                                \
      int N, bfam_subdomain_dgx_t *sub, const char *rate_prefix,               \
      const char *field_prefix, const bfam_long_real_t t);                     \
  void beard_dgx_inter_rhs_slip_weakening_interface_3_##order(                 \
      int N, bfam_subdomain_dgx_t *sub, const char *rate_prefix,               \
      const char *minus_rate_prefix, const char *field_prefix,                 \
      const bfam_long_real_t t);                                               \
  void beard_dgx_inter_rhs_ageing_law_interface_3_##order(                     \
      int N, bfam_subdomain_dgx_t *sub, const char *rate_prefix,               \
      const char *minus_rate_prefix, const char *field_prefix,                 \
      const bfam_long_real_t t);                                               \
  void beard_dgx_energy_3_##order(int N, bfam_real_t *energy_sq,               \
                                  bfam_subdomain_dgx_t *sub,                   \
                                  const char *field_prefix);
BFAM_LIST_OF_DGX_NORDERS
#undef X

void beard_dgx_intra_rhs_elastic_3_(int N, bfam_subdomain_dgx_t *sub,
                                    const char *rate_prefix,
                                    const char *field_prefix,
                                    const bfam_long_real_t t);
void beard_dgx_intra_rhs_sponge_3_(int N, bfam_subdomain_dgx_t *sub,
                                   const char *rate_prefix,
                                   const char *field_prefix,
                                   const bfam_long_real_t t);
void beard_dgx_rupture_time_3_(int N, bfam_subdomain_dgx_t *sub,
                               const char *field_prefix,
                               const bfam_long_real_t t,
                               const bfam_real_t Vrup);
void beard_dgx_duvaut_lions_return_map_3_(int N, bfam_subdomain_dgx_t *sub,
                                          const char *field_prefix,
                                          const bfam_long_real_t t,
                                          const bfam_long_real_t dt);
void beard_dgx_scale_rates_elastic_3_(int N, bfam_subdomain_dgx_t *sub,
                                      const char *rate_prefix,
                                      const bfam_long_real_t a);
void beard_dgx_scale_rates_interface_3_(int N, bfam_subdomain_dgx_t *sub,
                                        const char *rate_prefix,
                                        const bfam_long_real_t a,
                                        const char **f_names);
void beard_dgx_scale_rates_3_(int N, bfam_subdomain_dgx_t *sub,
                              const char *rate_prefix, const bfam_long_real_t a,
                              bfam_dictionary_t *fields, const char **scals,
                              const char **vects, const char **tens,
                              const bfam_locidx_t num_pts);
void beard_dgx_add_rates_elastic_3_(int N, bfam_subdomain_dgx_t *sub,
                                    const char *field_prefix_lhs,
                                    const char *field_prefix_rhs,
                                    const char *rate_prefix,
                                    const bfam_long_real_t a);
void beard_dgx_add_rates_interface_3_(int N, bfam_subdomain_dgx_t *sub,
                                      const char *field_prefix_lhs,
                                      const char *field_prefix_rhs,
                                      const char *rate_prefix,
                                      const bfam_long_real_t a,
                                      const char **f_names);
void beard_dgx_add_rates_glue_p_3_(int N, bfam_subdomain_dgx_t *sub,
                                   const char *field_prefix_lhs,
                                   const char *field_prefix_rhs,
                                   const char *rate_prefix,
                                   const bfam_long_real_t a, const char **scals,
                                   const char **vects, const char **tens);
void beard_dgx_add_rates_3_(int N, bfam_subdomain_dgx_t *sub,
                            const char *field_prefix_lhs,
                            const char *field_prefix_rhs,
                            const char *rate_prefix, const bfam_long_real_t a,
                            bfam_dictionary_t *fields, const char **scals,
                            const char **vects, const char **tens,
                            const bfam_locidx_t num_pts);
void beard_dgx_inter_rhs_boundary_3_(int N, bfam_subdomain_dgx_t *sub,
                                     const char *rate_prefix,
                                     const char *field_prefix,
                                     const bfam_long_real_t t,
                                     const bfam_real_t R);
void beard_dgx_inter_rhs_interface_3_(int N, bfam_subdomain_dgx_t *sub,
                                      const char *rate_prefix,
                                      const char *field_prefix,
                                      const bfam_long_real_t t);
void beard_dgx_inter_rhs_slip_weakening_interface_3_(
    int N, bfam_subdomain_dgx_t *sub, const char *rate_prefix,
    const char *minus_rate_prefix, const char *field_prefix,
    const bfam_long_real_t t);
void beard_dgx_inter_rhs_ageing_law_interface_3_(int N,
                                                 bfam_subdomain_dgx_t *sub,
                                                 const char *rate_prefix,
                                                 const char *minus_rate_prefix,
                                                 const char *field_prefix,
                                                 const bfam_long_real_t t);
void beard_dgx_energy_3_(int N, bfam_real_t *energy_sq,
                         bfam_subdomain_dgx_t *sub, const char *field_prefix);

#endif
