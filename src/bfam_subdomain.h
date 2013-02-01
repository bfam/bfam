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
  // comm
  char*           name;     /**< Name of the subdomain */
  bfam_mpicomm_t* comm;     /**< communicator for this subdomain */
  int             hasWorkl; /**< boolean for whether or not I can do work on this processor */
} bfam_subdomain_t;

/**
 * structure containing a domain (which is a collection of subdomains!)
 */
typedef struct bfam_domain 
{
  bfam_subdomain_t* subdomains; /**< array of pointers to subdomains */
  int numSubdomains;         /**< number of subdomains that are currently in the domain */
  int sizeSubdomains;        /**< total number of subdomains the domain can
                               hold, i.e.  size of the array*/
  bfam_mpicomm_t * comm;     /**< communicator for the whole domain */
  bfam_critbit0_tree_t name2num; /**< critbit map for convertings subdomain names to numbers */
} bfam_domain_t;

// /** Create a new subdomain
//  * 
//  * \param [in,out] domain pointer to the domain
//  * \param [in]     domComm pointer to the communicator for the domain
//  */
void bfam_domain_new(bfam_domain_t *domain, bfam_mpicomm_t *domComm);

// /** free up the memory allocated by the domain
//  * 
//  * \param [in,out] domain to clean up
//  */
void bfam_domain_free(bfam_domain_t *domain);

#endif
