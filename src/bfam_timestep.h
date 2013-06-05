#ifndef BFAM_TIMESTEP_H
#define BFAM_TIMESTEP_H

#include <bfam_base.h>
#include <bfam_domain.h>

/**
 * structure comtaining the necessary features of a time step routine
 */
typedef struct bfam_timestep
{
  bfam_domain_t* domain;         /**< my domain */
} bfam_timestep_t;

#endif
