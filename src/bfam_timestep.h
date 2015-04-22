#ifndef BFAM_TIMESTEP_H
#define BFAM_TIMESTEP_H

#include <bfam_base.h>
#include <bfam_domain.h>

/**
 * structure comtaining the necessary features of a time step routine
 */
struct bfam_ts;

typedef void (*bfam_ts_step_t)(struct bfam_ts *ts, bfam_long_real_t dt,
                               void *user_data);

typedef void (*step_extended_t)(struct bfam_ts *a_ts, bfam_long_real_t dt,
                                const char *rate_prefix,
                                const char *field_prefix_lhs,
                                const char *field_prefix_rhs, void *user_data);

/* scale rates function */
typedef void (*scale_rates_t)(bfam_subdomain_t *thisSubdomain,
                              const char *rate_prefix, const bfam_long_real_t a,
                              void *user_data);

/* compute rhs that does not require communication */
typedef void (*intra_rhs_t)(bfam_subdomain_t *thisSubdomain,
                            const char *rate_prefix,
                            const char *minus_rate_prefix,
                            const char *field_prefix, const bfam_long_real_t t,
                            void *user_data);

/* compute rhs that does require communication */
typedef void (*inter_rhs_t)(bfam_subdomain_t *thisSubdomain,
                            const char *rate_prefix,
                            const char *minus_rate_prefix,
                            const char *field_prefix, const bfam_long_real_t t,
                            void *user_data);

/* add the rates to the fields: q_lhs := q_rhs + a*dq */
/* NOTE: should handle case of in place addition */
typedef void (*add_rates_t)(bfam_subdomain_t *thisSubdomain,
                            const char *field_prefix_lhs,
                            const char *field_prefix_rhs,
                            const char *rate_prefix, const bfam_long_real_t a,
                            void *user_data);

typedef void (*aux_rates_t)(bfam_subdomain_t *thisSubdomain, const char *prefix,
                            void *user_data);

typedef struct bfam_ts
{
  bfam_domain_t *domain; /**< my domain */

  bfam_long_real_t tmp; /* field forces alignment of pointer. Better soln?*/

  /* do a step of size dt with time stepper ts*/
  bfam_ts_step_t step;
} bfam_ts_t;

/** initialize a time step routine
 *
 * \warning It is the callers responsibility to ensure that
 *          \a dom is freed after the time stepper is
 *
 * \param [in,out]  ts       pointer to time stepper to initialize
 * \param [in]      dom      pointer to the domain
 */
void bfam_ts_init(bfam_ts_t *ts, bfam_domain_t *dom);

/** free a time step routine
 *
 * \param [in,out]  ts       pointer to time stepper to free
 */
void bfam_ts_free(bfam_ts_t *ts);

#endif
