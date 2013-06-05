#ifndef BFAM_TIMESTEP_LSRK_H
#define BFAM_TIMESTEP_LSRK_H

#include <bfam_base.h>
#include <bfam_domain.h>
#include <bfam_timestep.h>

/**
 * structure comtaining the necessary features of a low memory RK scheme
 */
typedef struct bfam_timestep_lsrk
{
  bfam_timestep_t* p_timestep;         /**< parent timestepper */
} bfam_timestep_lsrk_t;

#endif
