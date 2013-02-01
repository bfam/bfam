#ifndef BFAM_DOMAIN_H
#define BFAM_DOMAIN_H

#include <bfam_mpicomm.h>
#include <bfam_critbit.h>
#include <bfam_subdomain.h>

/**
 * structure containing a domain (which is a collection of subdomains!)
 */
typedef struct bfam_domain 
{
  bfam_subdomain_t** subdomains;  /**< array of pointers to subdomains */
  int numSubdomains;             /**< number of subdomains that are currently in the domain */
  int sizeSubdomains;            /**< total number of subdomains the domain can
                                      hold, i.e.  size of the array*/
  bfam_mpicomm_t * comm;         /**< communicator for the whole domain */
  bfam_critbit0_tree_t name2num; /**< critbit map for convertings subdomain names to numbers */
} bfam_domain_t;


/* Domain based functions */
/** create a domain
 * 
 * \param [out] newDomain the newly created domain
 * \param [in]  domComm   pointer to the communicator for the domain
 */
bfam_domain_t* bfam_domain_new(bfam_mpicomm_t *domComm);

/** initializes a domain
 * 
 * \param [in,out] domain pointer to the domain
 * \param [in]     domComm pointer to the communicator for the domain
 */
void bfam_domain_init(bfam_domain_t *domain, bfam_mpicomm_t *domComm);

/** free up the memory allocated by the domain
 * 
 * \param [in,out] domain domain to clean up
 */
void bfam_domain_free(bfam_domain_t *domain);

/** Add subdomain
 * 
 * \param [in,out] thisDomain domain to add subdomain to
 * \param [in]     newSubdomain subdomain to add to the domain
 */
void bfam_domain_add_subdomain(bfam_domain_t *thisDomain, bfam_subdomain_t *newSubdomain);

#endif
