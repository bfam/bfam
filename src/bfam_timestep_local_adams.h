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
  bfam_ts_t base;                   /**< parent timestepper */
  int nStages;                      /**< number of steps */
  bfam_locidx_t *currentStageArray; /**< array of current stage counter */
  int numSteps;                     /**< number of steps completed */
  bfam_locidx_t numLevels;          /**< number of levels to time step */
  bfam_long_real_t  t;              /**< domain time */
  bfam_communicator_t **comm_array; /**< NULL terminated array of communicators
                                         for the local level communicator I
                                         handle
                                      */
  bfam_dictionary_t elems;    /**< dictionary of subdomains I step */

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

/** create an explicit Local Adams scheme
 *
 * \warning It is the callers responsibility to ensure that
 *          \a dom is freed after this Local Adams is
 *
 * \param [in]  dom              pointer to the domain
 * \param [in]  method           Adams tyoe we are using
 * \param [in]  max_level        number of levels to time step
 * \param [in]  subdom_match     match type for subdomains
 * \param [in]  subdom_tags      tags for the subdomains to time step
 * \param [in]  comm_match       match type for communication
 * \param [in]  comm_tags        tags for the communication required for RHS
 * \param [in]  mpicomm          MPI_Comm to use for communication
 * \param [in]  mpitag           tag to use for MPI communcation
 * \param [in]  comm_data        user data passed to the communicator new
 * \param [in]  aux_rates        create rate field with given prefix
 * \param [in]  glue_rates       create glue rate fields with given prefix
 * \param [in]  scale_rates      scale rates function
 * \param [in]  intra_rhs        function handle to intra RHS routine
 * \param [in]  inter_rhs        function handle to inter RHS routine
 * \param [in]  add_rates        function handle to add rates routine
 * \param [in]  RK_init          boolean which if true signifies using LSRK to
 *                               init Local Adams method
 *
 * \return the newly created low storage RK time stepper
 */
bfam_ts_local_adams_t*
bfam_ts_local_adams_new(bfam_domain_t* dom, bfam_ts_local_adams_method_t method,
    bfam_locidx_t max_level, bfam_domain_match_t subdom_match, const char**
    subdom_tags, bfam_domain_match_t comm_match, const char** comm_tags,
    MPI_Comm mpicomm, int mpitag, void * comm_data,
    void (*aux_rates)  (bfam_subdomain_t *thisSubdomain, const char *prefix),
    void (*glue_rates) (bfam_subdomain_t *thisSubdomain, const char *prefix),
    void (*scale_rates) (bfam_subdomain_t *thisSubdomain,
      const char *rate_prefix, const bfam_long_real_t a),
    void (*intra_rhs) (bfam_subdomain_t *thisSubdomain,
      const char *rate_prefix, const char *field_prefix,
      const bfam_long_real_t t),
    void (*inter_rhs) (bfam_subdomain_t *thisSubdomain,
      const char *rate_prefix, const char *field_prefix,
      const bfam_long_real_t t),
    void (*add_rates) (bfam_subdomain_t *thisSubdomain,
      const char *field_prefix_lhs, const char *field_prefix_rhs,
      const char *rate_prefix, const bfam_long_real_t a),
    const int RK_init);

/** initialize an Local Adams scheme
 *
 * \warning It is the callers responsibility to ensure that
 *          \a dom is freed after this Local Adams is
 *
 * \param [in,out]  ts           pointer to time stepper to initialize
 * \param [in]  dom              pointer to the domain
 * \param [in]  method           Low storage RK tyoe we are using
 * \param [in]  max_level        number of levels to time step
 * \param [in]  subdom_match     match type for subdomains
 * \param [in]  subdom_tags      tags for the subdomains to time step
 * \param [in]  comm_match       match type for communication
 * \param [in]  comm_tags        tags for the communication required for RHS
 * \param [in]  mpicomm          MPI_Comm to use for communication
 * \param [in]  mpitag           tag to use for MPI communcation
 * \param [in]  comm_data        user data passed to the communicator new
 * \param [in]  aux_rates        create rate field with given prefix
 * \param [in]  glue_rates       create glue rate fields with given prefix
 * \param [in]  scale_rates      scale rates function
 * \param [in]  intra_rhs        function handle to intra RHS routine
 * \param [in]  inter_rhs        function handle to inter RHS routine
 * \param [in]  add_rates        function handle to add rates routine
 * \param [in]  RK_init          boolean which if true signifies using LSRK to
 *                               init Local Adams method
 */
void
bfam_ts_local_adams_init(bfam_ts_local_adams_t* ts,
    bfam_domain_t* dom, bfam_ts_local_adams_method_t method, bfam_locidx_t
    max_level, bfam_domain_match_t subdom_match, const char** subdom_tags,
    bfam_domain_match_t comm_match, const char** comm_tags, MPI_Comm mpicomm,
    int mpitag, void * comm_data,
    void (*aux_rates) (bfam_subdomain_t *thisSubdomain, const char *prefix),
    void (*glue_rates) (bfam_subdomain_t *thisSubdomain, const char *prefix),
    void (*scale_rates) (bfam_subdomain_t *thisSubdomain,
      const char *rate_prefix, const bfam_long_real_t a),
    void (*intra_rhs) (bfam_subdomain_t *thisSubdomain,
      const char *rate_prefix, const char *field_prefix,
      const bfam_long_real_t t),
    void (*inter_rhs) (bfam_subdomain_t *thisSubdomain,
      const char *rate_prefix, const char *field_prefix,
      const bfam_long_real_t t),
    void (*add_rates) (bfam_subdomain_t *thisSubdomain,
      const char *field_prefix_lhs, const char *field_prefix_rhs,
      const char *rate_prefix, const bfam_long_real_t a),
    const int RK_init);

/** free an Local Adams scheme
 *
 * \param [in,out]  ts       pointer to time stepper to free
 */
void
bfam_ts_local_adams_free(bfam_ts_local_adams_t* ts);

/** set the time of the scheme
 *
 * \param [in,out]  ts       pointer to time stepper to set
 * \param [in]      time     time to set
 */
void
bfam_ts_local_adams_set_time(bfam_ts_local_adams_t* ts, bfam_long_real_t time);


/** get the time of the scheme
 *
 * \param [in]  ts       pointer to local_adams to get time
 */
bfam_long_real_t
bfam_ts_local_adams_get_time(bfam_ts_local_adams_t* ts);

#endif
