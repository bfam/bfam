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
  bfam_long_real_t *A;       /**< low memory RK A: rate scale */
  bfam_long_real_t *B;       /**< low memory RK B: update scale */
  bfam_long_real_t *C;       /**< low memory RK C: time scale*/
  int nStages;               /**< number of stages */
  bfam_long_real_t t;        /**< domain time */
  bfam_communicator_t *comm; /**< communicator I handle */
  bfam_dictionary_t elems;   /**< dictionary of subdomains I step */

  update_soln_t update_soln;
  scale_rates_t scale_rates;
  intra_rhs_t intra_rhs;
  inter_rhs_t inter_rhs;
  add_rates_t add_rates;
  step_extended_t step_extended;
} bfam_ts_lsrk_t;

typedef enum bfam_ts_lsrk_method
{
  BFAM_TS_LSRK_KC54,
  BFAM_TS_LSRK_FE,
  BFAM_TS_LSRK_HEUN,
  BFAM_TS_LSRK_W33,
  BFAM_TS_LSRK_NOOP,
} bfam_ts_lsrk_method_t;

/** create a low storage RK scheme with update kernel
 *
 * \warning It is the callers responsibility to ensure that
 *          \a dom is freed after this LSRK is
 *
 * \param [in]  dom              pointer to the domain
 * \param [in]  method           Low storage RK tyoe we are using
 * \param [in]  subdom_match     match type for subdomains
 * \param [in]  subdom_tags      tags for the subdomains to time step
 * \param [in]  comm_match       match type for communication
 * \param [in]  comm_tags        tags for the communication required for RHS
 * \param [in]  mpicomm          MPI_Comm to use for communication
 * \param [in]  mpitag           tag to use for MPI communcation
 * \param [in]  comm_late_send   late send tag to pass to communicator
 * \param [in]  comm_data        user data passed to the communicator new
 * \param [in]  update_sol       function handle to update
 * \param [in]  aux_rates        create rate field with given prefix
 * \param [in]  intra_rhs        function handle to intra RHS routine
 * \param [in]  inter_rhs        function handle to inter RHS routine
 * \param [in]  add_rates        function handle to add rates routine
 * \param [in]  user_data        pointer to data to be passed to aux_rates
 *
 * \return the newly created low storage RK time stepper
 */
bfam_ts_lsrk_t *
bfam_ts_lsrk_new_up(bfam_domain_t *dom, bfam_ts_lsrk_method_t method,
                    bfam_domain_match_t subdom_match, const char **subdom_tags,
                    bfam_domain_match_t comm_match, const char **comm_tags,
                    MPI_Comm mpicomm, int mpitag, int comm_late_send,
                    void *comm_data, update_soln_t update_soln,
                    aux_rates_t aux_rates, intra_rhs_t intra_rhs,
                    inter_rhs_t inter_rhs, void *user_data);

/** create a low storage RK scheme
 *
 * \warning It is the callers responsibility to ensure that
 *          \a dom is freed after this LSRK is
 *
 * \param [in]  dom              pointer to the domain
 * \param [in]  method           Low storage RK tyoe we are using
 * \param [in]  subdom_match     match type for subdomains
 * \param [in]  subdom_tags      tags for the subdomains to time step
 * \param [in]  comm_match       match type for communication
 * \param [in]  comm_tags        tags for the communication required for RHS
 * \param [in]  mpicomm          MPI_Comm to use for communication
 * \param [in]  mpitag           tag to use for MPI communcation
 * \param [in]  comm_late_send   late send tag to pass to communicator
 * \param [in]  comm_data        user data passed to the communicator new
 * \param [in]  aux_rates        create rate field with given prefix
 * \param [in]  scale_rates      function handle to scale_rates function
 * \param [in]  intra_rhs        function handle to intra RHS routine
 * \param [in]  inter_rhs        function handle to inter RHS routine
 * \param [in]  add_rates        function handle to add rates routine
 * \param [in]  user_data        pointer to data to be passed to aux_rates
 *
 * \return the newly created low storage RK time stepper
 */
bfam_ts_lsrk_t *
bfam_ts_lsrk_new(bfam_domain_t *dom, bfam_ts_lsrk_method_t method,
                 bfam_domain_match_t subdom_match, const char **subdom_tags,
                 bfam_domain_match_t comm_match, const char **comm_tags,
                 MPI_Comm mpicomm, int mpitag, int comm_late_send,
                 void *comm_data, aux_rates_t aux_rates,
                 scale_rates_t scale_rates, intra_rhs_t intra_rhs,
                 inter_rhs_t inter_rhs, add_rates_t add_rates, void *user_data);

/** initialize a low storage RK scheme
 *
 * \warning It is the callers responsibility to ensure that
 *          \a dom is freed after this LSRK is
 *
 * \param [in,out]  ts           pointer to time stepper to initialize
 * \param [in]  dom              pointer to the domain
 * \param [in]  method           Low storage RK tyoe we are using
 * \param [in]  subdom_match     match type for subdomains
 * \param [in]  subdom_tags      tags for the subdomains to time step
 * \param [in]  comm_match       match type for communication
 * \param [in]  comm_tags        tags for the communication required for RHS
 * \param [in]  mpicomm          MPI_Comm to use for communication
 * \param [in]  mpitag           tag to use for MPI communcation
 * \param [in]  comm_late_send   late send tag to pass to communicator
 * \param [in]  comm_data        user data passed to the communicator new
 * \param [in]  aux_rates        create rate field with given prefix
 * \param [in]  scale_rates      function handle to scale_rates function
 * \param [in]  intra_rhs        function handle to intra RHS routine
 * \param [in]  inter_rhs        function handle to inter RHS routine
 * \param [in]  add_rates        function handle to add rates routine
 * \param [in]  user_data        pointer to data to be passed to aux_rates
 */
void bfam_ts_lsrk_init(bfam_ts_lsrk_t *ts, bfam_domain_t *dom,
                       bfam_ts_lsrk_method_t method,
                       bfam_domain_match_t subdom_match,
                       const char **subdom_tags, bfam_domain_match_t comm_match,
                       const char **comm_tags, MPI_Comm mpicomm, int mpitag,
                       int comm_late_send, void *comm_data,
                       aux_rates_t aux_rates, scale_rates_t scale_rates,
                       intra_rhs_t intra_rhs, inter_rhs_t inter_rhs,
                       add_rates_t add_rates, void *user_data);

/** create a low storage RK scheme
 *
 * \warning It is the callers responsibility to ensure that
 *          \a dom is freed after this LSRK is
 *
 * \param [in]  dom              pointer to the domain
 * \param [in]  method           Low storage RK tyoe we are using
 * \param [in]  subdom_match     match type for subdomains
 * \param [in]  subdom_tags      tags for the subdomains to time step
 * \param [in]  comm_match       match type for communication
 * \param [in]  comm_tags        tags for the communication required for RHS
 * \param [in]  mpicomm          MPI_Comm to use for communication
 * \param [in]  mpitag           tag to use for MPI communcation
 * \param [in]  comm_late_send   late send tag to pass to communicator
 * \param [in]  comm_data        user data passed to the communicator new
 * \param [in]  aux_rates        create rate field with given prefix
 * \param [in]  scale_rates      function handle to scale_rates function
 * \param [in]  intra_rhs        function handle to intra RHS routine
 * \param [in]  inter_rhs        function handle to inter RHS routine
 * \param [in]  add_rates        function handle to add rates routine
 * \param [in]  make_rates       boolean for whether LSRK should make rates
 * \param [in]  user_data        pointer to data to be passed to aux_rates
 *
 * \return the newly created low storage RK time stepper
 */
bfam_ts_lsrk_t *bfam_ts_lsrk_new_extended(
    bfam_domain_t *dom, bfam_ts_lsrk_method_t method,
    bfam_domain_match_t subdom_match, const char **subdom_tags,
    bfam_domain_match_t comm_match, const char **comm_tags, MPI_Comm mpicomm,
    int mpitag, int comm_late_send, void *comm_data, aux_rates_t aux_rates,
    scale_rates_t scale_rates, intra_rhs_t intra_rhs, inter_rhs_t inter_rhs,
    add_rates_t add_rates, int make_rates, void *user_data);

/** initialize a low storage RK scheme
 *
 * \warning It is the callers responsibility to ensure that
 *          \a dom is freed after this LSRK is
 *
 * \param [in,out]  ts           pointer to time stepper to initialize
 * \param [in]  dom              pointer to the domain
 * \param [in]  method           Low storage RK tyoe we are using
 * \param [in]  subdom_match     match type for subdomains
 * \param [in]  subdom_tags      tags for the subdomains to time step
 * \param [in]  comm_match       match type for communication
 * \param [in]  comm_tags        tags for the communication required for RHS
 * \param [in]  mpicomm          MPI_Comm to use for communication
 * \param [in]  mpitag           tag to use for MPI communcation
 * \param [in]  comm_late_send   late send tag to pass to communicator
 * \param [in]  comm_data        user data passed to the communicator new
 * \param [in]  update_soln      add rates to the soln and scale the rates
 * \param [in]  aux_rates        create rate field with given prefix
 * \param [in]  scale_rates      function handle to scale_rates function
 * \param [in]  intra_rhs        function handle to intra RHS routine
 * \param [in]  inter_rhs        function handle to inter RHS routine
 * \param [in]  add_rates        function handle to add rates routine
 * \param [in]  make_rates       boolean for whether LSRK should make rates
 * \param [in]  user_data        pointer to data to be passed to aux_rates
 */
void bfam_ts_lsrk_init_extended(
    bfam_ts_lsrk_t *ts, bfam_domain_t *dom, bfam_ts_lsrk_method_t method,
    bfam_domain_match_t subdom_match, const char **subdom_tags,
    bfam_domain_match_t comm_match, const char **comm_tags, MPI_Comm mpicomm,
    int mpitag, int comm_late_send, void *comm_data, update_soln_t update_soln,
    aux_rates_t aux_rates, scale_rates_t scale_rates, intra_rhs_t intra_rhs,
    inter_rhs_t inter_rhs, add_rates_t add_rates, int make_rates,
    void *user_data);

/** free a low storage RK scheme
 *
 * \param [in,out]  ts       pointer to time stepper to free
 */
void bfam_ts_lsrk_free(bfam_ts_lsrk_t *ts);

/** set the time of the scheme
 *
 * \param [in,out]  ts       pointer to time stepper to set
 * \param [in]      time     time to set
 */
void bfam_ts_lsrk_set_time(bfam_ts_lsrk_t *ts, bfam_long_real_t time);

/** get the time of the scheme
 *
 * \param [in]  ts       pointer to lsrk to get time
 */
bfam_long_real_t bfam_ts_lsrk_get_time(bfam_ts_lsrk_t *ts);

#endif
