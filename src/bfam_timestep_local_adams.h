#ifndef BFAM_TIMESTEP_LOCAL_ADAMS_H
#define BFAM_TIMESTEP_LOCAL_ADAMS_H

#include <bfam_base.h>
#include <bfam_domain.h>
#include <bfam_timestep.h>
#include <bfam_timestep_lsrk.h>
#include <bfam_communicator.h>

/**
 * structure comtaining the necessary features of an explicit, local Adams
 * method
 *
 * Each step is of the form
 * t_{n+1}  := t_n + dt
 * q_{n+1}  := q_{n} + dt \sum_{k=0}^{m} a_{k} dq_{n-k}
 * dq_{n+1} := RHS(q_{n+1},t_{n+1})
 *
 */
typedef struct bfam_ts_local_adams
{
  bfam_ts_t base;            /**< parent timestepper */
  int nStages;               /**< number of steps */
  int currentStage;          /**< current stage counter */
  int numSteps;              /**< number of steps completed */
  bfam_long_real_t  t;       /**< domain time */
  bfam_communicator_t *comm; /**< communicator I handle */
  bfam_dictionary_t elems;   /**< dictionary of subdomains I step */

  /* LSRK method for initialization */
  bfam_ts_lsrk_t* lsrk;

  /* scale rates function */
  void (*scale_rates) (bfam_subdomain_t *thisSubdomain,
      const char *rate_prefix, const bfam_long_real_t a);

  /* compute rhs that does not require communication */
  void (*intra_rhs) (bfam_subdomain_t *thisSubdomain,
      const char *rate_prefix, const char *field_prefix,
      const bfam_long_real_t t);

  /* compute rhs that does require communication */
  void (*inter_rhs) (bfam_subdomain_t *thisSubdomain,
      const char *rate_prefix, const char *field_prefix,
      const bfam_long_real_t t);

  /* add the rates to the fields: q_lhs := q_rhs + a*dq */
  /* NOTE: should handle case of in place addition */
  void (*add_rates) (bfam_subdomain_t *thisSubdomain,
      const char *field_prefix_lhs, const char *field_prefix_rhs,
      const char *rate_prefix, const bfam_long_real_t a);
} bfam_ts_local_adams_t;

typedef enum bfam_ts_local_adams_method
{
  BFAM_TS_LOCAL_ADAMS_1,
  BFAM_TS_LOCAL_ADAMS_2,
  BFAM_TS_LOCAL_ADAMS_3,
  BFAM_TS_LOCAL_ADAMS_4,
  BFAM_TS_LOCAL_ADAMS_NOOP,
} bfam_ts_local_adams_method_t;

/** fill the level tag for a given time level
 *
 * \param [out]  tag        string to fill with tag
 * \param [in]   buf_siz    size of the buffer to fill
 * \param [in]   level      level number to fill
 */
void
bfam_ts_local_adams_fill_level_tag(char* tag, size_t buf_siz, int level);

#endif
