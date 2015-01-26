#ifndef BFAM_DOMAIN_H
#define BFAM_DOMAIN_H

#include <bfam_base.h>
#include <bfam_critbit.h>
#include <bfam_dictionary.h>
#include <bfam_subdomain.h>

#define BFAM_DEFAULT_SUBDOMAIN_ROOT 0

/**
 * structure containing a domain (which is a collection of subdomains!)
 */
typedef struct bfam_domain
{
  bfam_subdomain_t **subdomains; /**< array of pointers to subdomains */
  bfam_locidx_t numSubdomains;   /**< number of subdomains that are
                                      currently in the domain */
  bfam_locidx_t sizeSubdomains;  /**< total number of subdomains the domain
                                      can hold, i.e.  size of the array*/
  MPI_Comm comm;                 /**< communicator for the whole domain */
  bfam_dictionary_t name2num;    /**< dictionary map for convertings
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
bfam_domain_t *bfam_domain_new(MPI_Comm domComm);

/** initializes a domain
 *
 * \param [in,out] domain pointer to the domain
 * \param [in]     domComm pointer to the communicator for the domain
 */
void bfam_domain_init(bfam_domain_t *domain, MPI_Comm domComm);

/** Clean up domain
 *
 * frees any mememory allocated by the domain and calls free command on all
 * subdomains
 *
 * \param [in,out] domain domain to clean up
 */
void bfam_domain_free(bfam_domain_t *domain);

/** Add subdomain
 *
 * \param [in,out] thisDomain domain to add subdomain to
 * \param [in]     newSubdomain subdomain to add to the domain
 *
 * \return subdomain id
 */
bfam_locidx_t bfam_domain_add_subdomain(bfam_domain_t *thisDomain,
                                        bfam_subdomain_t *newSubdomain);

/** Get subdomain
 *
 * \param [in]     id Id of subdomain to get
 *
 * \return subdomain pointer
 */
bfam_subdomain_t *bfam_domain_get_subdomain_by_num(bfam_domain_t *thisDomain,
                                                   bfam_locidx_t id);

/** Get subdomains with tags passed in
 *
 * \param [in]  thisDomain    domain to search for subdomains in
 * \param [in]  matchType     type of match, \c BFAM_DOMAIN_OR will
 *                            match subdomains with any of the tags
 *                            and \c BFAM_DOMAIN_AND will match subdomains
 *                            with all of the tags.
 * \param [in]  tags          \c NULL terminated array of the tags to match
 * \param [in]  numEntries    number of entries in the \a subdomains array
 * \param [out] subdomains    array of pointers to be filled with matching
 *                            subdomains
 * \param [out] numSubdomains number of matching subdomains
 *
 */
void bfam_domain_get_subdomains(bfam_domain_t *thisDomain,
                                bfam_domain_match_t match, const char **tags,
                                bfam_locidx_t numEntries,
                                bfam_subdomain_t **subdomains,
                                bfam_locidx_t *numSubdomains);

/** Get subdomains with tags passed in
 *
 * \param [in]  thisDomain    domain to search for subdomains in
 * \param [in]  matchType     type of match, \c BFAM_DOMAIN_OR will
 *                            match subdomains with any of the tags
 *                            and \c BFAM_DOMAIN_AND will match subdomains
 *                            with all of the tags.
 * \param [in]  tags          critbit tree of tags
 * \param [in]  numEntries    number of entries in the \a subdomains array
 * \param [out] subdomains    array of pointers to be filled with matching
 *                            subdomains
 * \param [out] numSubdomains number of matching subdomains
 *
 */
void bfam_domain_get_subdomains_critbit(bfam_domain_t *thisDomain,
                                        bfam_domain_match_t match,
                                        bfam_critbit0_tree_t *tags,
                                        bfam_locidx_t numEntries,
                                        bfam_subdomain_t **subdomains,
                                        bfam_locidx_t *numSubdomains);

/** Add tag to subdomains matching the tags passed in.
 *
 * \param [in,out]  thisDomain    domain to search for subdomains in
 * \param [in]      matchType     type of match, \c BFAM_DOMAIN_OR will
 *                                match subdomains with any of the tags
 *                                and \c BFAM_DOMAIN_AND will match subdomains
 *                                with all of the tags.
 * \param [in]      mtags         \c NULL terminated array of the tags to match
 * \param [in]      tag           tag to add to the subdomains
 *
 */
void bfam_domain_add_tag(bfam_domain_t *thisDomain, bfam_domain_match_t match,
                         const char **mtags, const char *tag);

/** Add fields to subdomains matching the tags passed in.
 *
 * \param [in,out]  thisDomain    domain to search for subdomains in
 * \param [in]      matchType     type of match, \c BFAM_DOMAIN_OR will
 *                                match subdomains with any of the tags
 *                                and \c BFAM_DOMAIN_AND will match subdomains
 *                                with all of the tags.
 * \param [in]      mtags         \c NULL terminated array of the tags to match
 * \param [in]      tags          tags to add to the subdomains
 *
 */
void bfam_domain_add_tags(bfam_domain_t *thisDomain, bfam_domain_match_t match,
                          const char **mtags, const char **tags);

/** Add fields to subdomains matching the tags passed in.
 *
 * \param [in,out]  thisDomain    domain to search for subdomains in
 * \param [in]      matchType     type of match, \c BFAM_DOMAIN_OR will
 *                                match subdomains with any of the tags
 *                                and \c BFAM_DOMAIN_AND will match subdomains
 *                                with all of the tags.
 * \param [in]      tags          \c NULL terminated array of the tags to match
 * \param [in]      field         field to add to the subdomains
 *
 */
void bfam_domain_add_field(bfam_domain_t *thisDomain, bfam_domain_match_t match,
                           const char **tags, const char *field);

/** Add fields to subdomains matching the tags passed in.
 *
 * \param [in,out]  thisDomain    domain to search for subdomains in
 * \param [in]      matchType     type of match, \c BFAM_DOMAIN_OR will
 *                                match subdomains with any of the tags
 *                                and \c BFAM_DOMAIN_AND will match subdomains
 *                                with all of the tags.
 * \param [in]      tags          \c NULL terminated array of the tags to match
 * \param [in]      fields        fields to add to the subdomains
 *
 */
void bfam_domain_add_fields(bfam_domain_t *thisDomain,
                            bfam_domain_match_t match, const char **tags,
                            const char **fields);

/** Add fields to subdomains matching the tags passed in.
 *
 * Here the tags are stored in a critbit.
 *
 * \param [in,out]  thisDomain    domain to search for subdomains in
 * \param [in]      matchType     type of match, \c BFAM_DOMAIN_OR will
 *                                match subdomains with any of the tags
 *                                and \c BFAM_DOMAIN_AND will match subdomains
 *                                with all of the tags.
 * \param [in]      tags          critbit of the tags to match
 * \param [in]      field         field to add to the subdomains
 *
 */
void bfam_domain_add_field_critbit(bfam_domain_t *thisDomain,
                                   bfam_domain_match_t match,
                                   bfam_critbit0_tree_t *tags,
                                   const char *field);

/** Add fields to subdomains matching the tags passed in
 *
 * \param [in,out]  thisDomain    domain to search for subdomains in
 * \param [in]      matchType     type of match, \c BFAM_DOMAIN_OR will
 *                                match subdomains with any of the tags
 *                                and \c BFAM_DOMAIN_AND will match subdomains
 *                                with all of the tags.
 * \param [in]      tags          critbit of the tags to match
 * \param [in]      fields        fields to add to the subdomains
 *
 */
void bfam_domain_add_fields_critbit(bfam_domain_t *thisDomain,
                                    bfam_domain_match_t match,
                                    bfam_critbit0_tree_t *tags,
                                    const char **fields);

/** Add minus fields to subdomains matching the tags passed in.
 *
 * \param [in,out]  thisDomain    domain to search for subdomains in
 * \param [in]      matchType     type of match, \c BFAM_DOMAIN_OR will
 *                                match subdomains with any of the tags
 *                                and \c BFAM_DOMAIN_AND will match subdomains
 *                                with all of the tags.
 * \param [in]      tags          \c NULL terminated array of the tags to match
 * \param [in]      field         field to add to the subdomains
 *
 */
void bfam_domain_add_minus_field(bfam_domain_t *thisDomain,
                                 bfam_domain_match_t match, const char **tags,
                                 const char *field);

/** Add minus fields to subdomains matching the tags passed in.
 *
 * \param [in,out]  thisDomain    domain to search for subdomains in
 * \param [in]      matchType     type of match, \c BFAM_DOMAIN_OR will
 *                                match subdomains with any of the tags
 *                                and \c BFAM_DOMAIN_AND will match subdomains
 *                                with all of the tags.
 * \param [in]      tags          \c NULL terminated array of the tags to match
 * \param [in]      fields        fields to add to the subdomains
 *
 */
void bfam_domain_add_minus_fields(bfam_domain_t *thisDomain,
                                  bfam_domain_match_t match, const char **tags,
                                  const char **fields);

/** Add minus fields to subdomains matching the tags passed in.
 *
 * Here the tags are stored in a critbit.
 *
 * \param [in,out]  thisDomain    domain to search for subdomains in
 * \param [in]      matchType     type of match, \c BFAM_DOMAIN_OR will
 *                                match subdomains with any of the tags
 *                                and \c BFAM_DOMAIN_AND will match subdomains
 *                                with all of the tags.
 * \param [in]      tags          critbit of the tags to match
 * \param [in]      field         field to add to the subdomains
 *
 */
void bfam_domain_add_minus_field_critbit(bfam_domain_t *thisDomain,
                                         bfam_domain_match_t match,
                                         bfam_critbit0_tree_t *tags,
                                         const char *field);

/** Add minus fields to subdomains matching the tags passed in
 *
 * \param [in,out]  thisDomain    domain to search for subdomains in
 * \param [in]      matchType     type of match, \c BFAM_DOMAIN_OR will
 *                                match subdomains with any of the tags
 *                                and \c BFAM_DOMAIN_AND will match subdomains
 *                                with all of the tags.
 * \param [in]      tags          critbit of the tags to match
 * \param [in]      fields        fields to add to the subdomains
 *
 */
void bfam_domain_add_minus_fields_critbit(bfam_domain_t *thisDomain,
                                          bfam_domain_match_t match,
                                          bfam_critbit0_tree_t *tags,
                                          const char **fields);

/** Add plus fields to subdomains matching the tags passed in.
 *
 * \param [in,out]  thisDomain    domain to search for subdomains in
 * \param [in]      matchType     type of match, \c BFAM_DOMAIN_OR will
 *                                match subdomains with any of the tags
 *                                and \c BFAM_DOMAIN_AND will match subdomains
 *                                with all of the tags.
 * \param [in]      tags          \c NULL terminated array of the tags to match
 * \param [in]      field         field to add to the subdomains
 *
 */
void bfam_domain_add_plus_field(bfam_domain_t *thisDomain,
                                bfam_domain_match_t match, const char **tags,
                                const char *field);

/** Add plus fields to subdomains matching the tags passed in.
 *
 * \param [in,out]  thisDomain    domain to search for subdomains in
 * \param [in]      matchType     type of match, \c BFAM_DOMAIN_OR will
 *                                match subdomains with any of the tags
 *                                and \c BFAM_DOMAIN_AND will match subdomains
 *                                with all of the tags.
 * \param [in]      tags          \c NULL terminated array of the tags to match
 * \param [in]      fields        fields to add to the subdomains
 *
 */
void bfam_domain_add_plus_fields(bfam_domain_t *thisDomain,
                                 bfam_domain_match_t match, const char **tags,
                                 const char **fields);

/** Add plus fields to subdomains matching the tags passed in.
 *
 * Here the tags are stored in a critbit.
 *
 * \param [in,out]  thisDomain    domain to search for subdomains in
 * \param [in]      matchType     type of match, \c BFAM_DOMAIN_OR will
 *                                match subdomains with any of the tags
 *                                and \c BFAM_DOMAIN_AND will match subdomains
 *                                with all of the tags.
 * \param [in]      tags          critbit of the tags to match
 * \param [in]      field         field to add to the subdomains
 *
 */
void bfam_domain_add_plus_field_critbit(bfam_domain_t *thisDomain,
                                        bfam_domain_match_t match,
                                        bfam_critbit0_tree_t *tags,
                                        const char *field);

/** Add fields to subdomains matching the tags passed in
 *
 * \param [in,out]  thisDomain    domain to search for subdomains in
 * \param [in]      matchType     type of match, \c BFAM_DOMAIN_OR will
 *                                match subdomains with any of the tags
 *                                and \c BFAM_DOMAIN_AND will match subdomains
 *                                with all of the tags.
 * \param [in]      tags          critbit of the tags to match
 * \param [in]      fields        fields to add to the subdomains
 *
 */
void bfam_domain_add_plus_fields_critbit(bfam_domain_t *thisDomain,
                                         bfam_domain_match_t match,
                                         bfam_critbit0_tree_t *tags,
                                         const char **fields);

/** Initialize a field.
 *
 * \param [in,out]  thisDomain    domain to search for subdomains in
 * \param [in]      matchType     type of match, \c BFAM_DOMAIN_OR will
 *                                match subdomains with any of the tags
 *                                and \c BFAM_DOMAIN_AND will match subdomains
 *                                with all of the tags.
 * \param [in]      tags          \c NULL terminated array of the tags to match
 * \param [in]      field         field to initialize
 * \param [in]      time          time to pass to initialization
 * \param [in]      init_field    field initialization function
 * \param [in]      arg           user pointer to pass to init function
 *
 */
void bfam_domain_init_field(bfam_domain_t *thisDomain,
                            bfam_domain_match_t match, const char **tags,
                            const char *field, bfam_real_t time,
                            bfam_subdomain_init_field_t init_field, void *arg);

/** Initialize a field.
 *
 * Here the tags are stored in a critbit.
 *
 * \param [in,out]  thisDomain    domain to search for subdomains in
 * \param [in]      matchType     type of match, \c BFAM_DOMAIN_OR will
 *                                match subdomains with any of the tags
 *                                and \c BFAM_DOMAIN_AND will match subdomains
 *                                with all of the tags.
 * \param [in]      tags          critbit of the tags to match
 * \param [in]      field         field to initialize
 * \param [in]      time          time to pass to initialization
 * \param [in]      init_field    field initialization function
 * \param [in]      arg           user pointer to pass to init function
 *
 */
void bfam_domain_init_field_critbit(bfam_domain_t *thisDomain,
                                    bfam_domain_match_t match,
                                    bfam_critbit0_tree_t *tags,
                                    const char *field, bfam_real_t time,
                                    bfam_subdomain_init_field_t init_field,
                                    void *arg);

#endif
