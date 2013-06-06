#ifndef BFAM_TIMESTEP_H
#define BFAM_TIMESTEP_H

#include <bfam_base.h>
#include <bfam_domain.h>

/**
 * structure comtaining the necessary features of a time step routine
 */
typedef struct bfam_ts
{
  bfam_domain_t* domain;         /**< my domain */
} bfam_ts_t;

/** initialize a time step routine
 *
 * \warning It is the callers responsibility to ensure that
 *          \a dom is freed after the time stepper is
 *
 * \param [in,out]  ts       pointer to time stepper to initialize
 * \param [in]      dom      pointer to the domain
 */
void
bfam_ts_init(bfam_ts_t* ts, bfam_domain_t* dom);

/** free a time step routine
 *
 * \param [in,out]  ts       pointer to time stepper to free
 */
void
bfam_ts_free(bfam_ts_t* ts);

#endif
