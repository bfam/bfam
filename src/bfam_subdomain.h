#ifndef BFAM_SUBDOMAIN_H
#define BFAM_SUBDOMAIN_H

#include <bfam_mpicomm.h>
#include <bfam_critbit.h>

/**
 * base structure for all subdomains types. Any new subdomain should have this
 * as its first member with the name base, i.e.,
 * \code{.c}
 * typedef struct new_subdomain_type
 * {
 *   bfam_subdomain_t base;
 *   ...
 * }
 */
typedef struct bfam_subdomain
{
  char*           name;     /**< Name of the subdomain */
  bfam_mpicomm_t* comm;     /**< communicator for this subdomain */
  int             hasWorkl; /**< boolean for whether or not I can do work on this processor */

  /* Function pointers that domain will need to call */
  /**< start communication */
  void (*start_communication) (struct bfam_subdomain *thisSubdomain);
  /**< end communication */
  void (*end_communication)   (struct bfam_subdomain *thisSubdomain);
  /**< start I/O */
  void (*start_io)            (struct bfam_subdomain *thisSubdomain);
  /**< end I/O */
  void (*start_end)           (struct bfam_subdomain *thisSubdomain);
  /**< do the work that can be done without communication */
  void (*do_internal_RHS)     (struct bfam_subdomain *thisSubdomain);
  /**< do the work that requires communication */
  void (*do_external_RHS)     (struct bfam_subdomain *thisSubdomain);
  /**< update solution */
  void (*update_fields)       (struct bfam_subdomain *thisSubdomain);
} bfam_subdomain_t;

#endif
