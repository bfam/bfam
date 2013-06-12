#ifndef BFAM_TIMESTEP_LSRK_H
#define BFAM_TIMESTEP_LSRK_H

#include <bfam_base.h>
#include <bfam_domain.h>
#include <bfam_timestep.h>
#include <bfam_communicator.h>

/**
 * structure comtaining the necessary features of a low memory RK scheme
 *
 * Each stages the form
 * t  := t0       + c*dt
 * dq := RHS(q,t) + a*dq
 * q  := q        + b*dt*dq
 */
typedef struct bfam_ts_lsrk
{
  bfam_ts_t base;            /**< parent timestepper */
  bfam_long_real_t* A;       /**< low memory RK A: rate scale */
  bfam_long_real_t* B;       /**< low memory RK B: update scale */
  bfam_long_real_t* C;       /**< low memory RK C: time scale*/
  int nStages;               /**< number of stages */
  bfam_long_real_t  t;       /**< domain time */
  bfam_communicator_t *comm; /**< communicator I handle */
  bfam_dictionary_t elems;   /**< dictionary of subdomains I step */
} bfam_ts_lsrk_t;

typedef enum bfam_ts_lsrk_method
{
  BFAM_TS_LSRK_KC54,
  BFAM_TS_LSRK_FE,
  BFAM_TS_LSRK_HEUN,
  BFAM_TS_LSRK_W33
} bfam_ts_lsrk_method_t;

/** create a low storage RK scheme
 *
 * \warning It is the callers responsibility to ensure that
 *          \a dom is freed after this LSRK is
 *
 * \param [in]  dom      pointer to the domain
 * \param [in]  comm     pointer to the communicator I use
 * \param [in]  method   Low storage RK tyoe we are using
 *
 * \return the newly created low storage RK time stepper
 */
bfam_ts_lsrk_t*
bfam_ts_lsrk_new(bfam_domain_t* dom,
    bfam_communicator_t *comm, bfam_ts_lsrk_method_t method);

/** initialize a low storage RK scheme
 *
 * \warning It is the callers responsibility to ensure that
 *          \a dom is freed after this LSRK is
 *
 * \param [in,out]  ts       pointer to time stepper to initialize
 * \param [in]      dom      pointer to the domain
 * \param [in]  comm     pointer to the communicator I use
 * \param [in]      method   Low storage RK tyoe we are using
 */
void
bfam_ts_lsrk_init(bfam_ts_lsrk_t* ts, bfam_domain_t* dom,
    bfam_communicator_t *comm, bfam_ts_lsrk_method_t method);

/** free a low storage RK scheme
 *
 * \param [in,out]  ts       pointer to time stepper to free
 */
void
bfam_ts_lsrk_free(bfam_ts_lsrk_t* ts);

/** set the time of the scheme
 *
 * \param [in,out]  ts       pointer to time stepper to set
 * \param [in]      time     time to set
 */
void
bfam_ts_lsrk_set_time(bfam_ts_lsrk_t* ts,bfam_long_real_t time);


/** get the time of the scheme
 *
 * \param [in]  ts       pointer to lsrk to get time
 */
bfam_long_real_t
bfam_ts_lsrk_get_time(bfam_ts_lsrk_t* ts);


/** add subdomains to this time stepper
 *
 * \param [in,out] ts            pointer to time stepper
 * \param [in]     matchType     type of match, \c BFAM_DOMAIN_OR will
 *                               match subdomains with any of the tags
 *                               and \c BFAM_DOMAIN_AND will match subdomains
 *                               with all of the tags.
 * \param [in]     tags          \c NULL terminated array of the tags to match
 * \param [in]     fields        \c NULL terminated array of fields names
 * \param [in]     scale_rates   function handle to use to scale rates
 * \param [in]     intra_rhs     function handle to use to function to use to
 *                               updae rhs function for non-communication work
 * \param [in]     inter_rhs     function handle to use to function to use to
 *                               updae rhs function for communication work
 * \param [in]     add_rates     function pointer to add rates back into
 *                               solution
 */
void
bfam_ts_lsrk_add_subdomains(bfam_ts_lsrk_t *ts, bfam_domain_match_t match,
                            const char *tags[], const char *fields[],
  void (*scale_rates) (bfam_subdomain_t *thisSubdomain, void *rates,
      const bfam_real_t a),
  void (*intra_rhs) (bfam_subdomain_t *thisSubdomain, void *rates,
      const void *fields, const bfam_real_t t),
  void (*inter_rhs) (bfam_subdomain_t *thisSubdomain, void *rates,
      const void *fields, const bfam_real_t t),
  void (*add_rates) (bfam_subdomain_t *thisSubdomain, void *fields,
      const void *rates, const bfam_real_t a));

#endif
