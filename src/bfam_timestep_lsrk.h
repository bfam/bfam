#ifndef BFAM_TIMESTEP_LSRK_H
#define BFAM_TIMESTEP_LSRK_H

#include <bfam_base.h>
#include <bfam_domain.h>
#include <bfam_timestep.h>

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
  bfam_ts_t* p_ts;                     /**< parent timestepper */
  bfam_real_t* A;                      /**< low memory RK A: rate scale */
  bfam_real_t* B;                      /**< low memory RK B: update scale */
  bfam_real_t* C;                      /**< low memory RK C: time scale*/
  int nStages;                         /**< number of stages */
  bfam_real_t  t;                      /**< domain time */
  bfam_real_t  dt;                     /**< domain dt   */
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
 * \param [in]  method   Low storage RK tyoe we are using
 *
 * \return the newly created low storage RK time stepper
 */
bfam_ts_lsrk_t*
bfam_ts_lsrk_new(bfam_domain_t* dom, bfam_ts_lsrk_method_t method);

/** initialize a low storage RK scheme
 *
 * \warning It is the callers responsibility to ensure that
 *          \a dom is freed after this LSRK is
 *
 * \param [in,out]  ts       pointer to time stepper to initialize
 * \param [in]      dom      pointer to the domain
 * \param [in]      method   Low storage RK tyoe we are using
 */
void
bfam_ts_lsrk_init(bfam_ts_lsrk_t* ts, bfam_domain_t* dom,
    bfam_ts_lsrk_method_t method);

#endif
