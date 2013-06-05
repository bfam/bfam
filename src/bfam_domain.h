#ifndef BFAM_DOMAIN_H
#define BFAM_DOMAIN_H

#include <bfam_base.h>
#include <bfam_mpicomm.h>
#include <bfam_critbit.h>
#include <bfam_subdomain.h>

/**
 * structure containing a domain (which is a collection of subdomains!)
 */
typedef struct bfam_domain 
{
  bfam_subdomain_t** subdomains;  /**< array of pointers to subdomains */
  bfam_locidx_t numSubdomains;   /**< number of subdomains that are
                                      currently in the domain */
  bfam_locidx_t sizeSubdomains;  /**< total number of subdomains the domain
                                      can hold, i.e.  size of the array*/
  bfam_mpicomm_t * comm;         /**< communicator for the whole domain */
  bfam_critbit0_tree_t name2num; /**< critbit map for convertings
                                      subdomain names to numbers */
} bfam_domain_t;

typedef enum bfam_domain_match
{
  BFAM_DOMAIN_AND,
  BFAM_DOMAIN_OR,
} bfam_domain_match_t;

/* Domain based functions */
/** create a domain
 *
 * \param [in]  domComm   pointer to the communicator for the domain
 *
 * \return the newly created domain
 */
bfam_domain_t* bfam_domain_new(bfam_mpicomm_t *domComm);

/** initializes a domain
 *
 * \param [in,out] domain pointer to the domain
 * \param [in]     domComm pointer to the communicator for the domain
 */
void
bfam_domain_init(bfam_domain_t *domain, bfam_mpicomm_t *domComm);

/** Clean up domain
 *
 * frees any mememory allocated by the domain and calls free command on all
 * subdomains
 *
 * \param [in,out] domain domain to clean up
 */
void
bfam_domain_free(bfam_domain_t *domain);

/** Add subdomain
 *
 * \param [in,out] thisDomain domain to add subdomain to
 * \param [in]     newSubdomain subdomain to add to the domain
 */
void
bfam_domain_add_subdomain(bfam_domain_t *thisDomain,
    bfam_subdomain_t *newSubdomain);

/** Get subdomains with tags passed in
 *
 * \param [in]  thisDomain    domain to search for subdomains in
 * \param [in]  matchType     type of match, \c BFAM_DOMAIN_OR will
 *                            match subdomains with any of the tags
 *                            and \c BFAM_DOMAIN_AND will match subdomains
 *                            with all of the tags.
 * \param [in]  numTags       number of tags to try to match
 * \param [in]  tags          array of the tags to match
 * \param [in]  numEntries    number of entries in the \a subdomains array
 * \param [out] subdomains    array of pointers to be filled with matching
 *                            subdomains
 * \param [out] numSubdomains number of matching subdomains
 *
 */
void
bfam_domain_get_subdomains(bfam_domain_t *thisDomain,
    bfam_domain_match_t match, size_t numTags, const char **tags,
    bfam_locidx_t numEntries, bfam_subdomain_t **subdomains,
    bfam_locidx_t *numSubdomains);

#endif
