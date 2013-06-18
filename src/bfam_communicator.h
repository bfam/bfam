#ifndef BFAM_COMMUNICATOR_H
#define BFAM_COMMUNICATOR_H

#include <bfam_domain.h>

typedef struct bfam_communicator
{
  MPI_Comm comm;
  bfam_locidx_t numSubdomains;
} bfam_communicator_t;

/** create a communicator
 *
 * \param [in] domain     domain to output to vtk files
 * \param [in] match      type of match, \c BFAM_DOMAIN_OR will
 *                        match subdomains with any of the tags
 *                        and \c BFAM_DOMAIN_AND will match subdomains
 *                        with all of the tags.
 * \param [in] tags       \c NULL terminated array of the tags to match
 *                        glue grids doing the communication
 * \param [in] comm       MPI communicator
 *
 * \return the newly created communicator
 */
bfam_communicator_t*
bfam_communicator_new(bfam_domain_t *domain, bfam_domain_match_t match,
    const char **tags, MPI_Comm comm);

/** initializes a communicator
 *
 * \param [in,out] communicator pointer to the communicator
 * \param [in]     domain       domain to output to vtk files
 * \param [in]     match        type of match, \c BFAM_DOMAIN_OR will
 *                              match subdomains with any of the tags
 *                              and \c BFAM_DOMAIN_AND will match subdomains
 *                              with all of the tags.
 * \param [in]     tags         \c NULL terminated array of the tags to match
 *                              glue grids doing the communication
 * \param [in]     comm         MPI communicator
 */
void
bfam_communicator_init(bfam_communicator_t* communicator,
    bfam_domain_t *domain, bfam_domain_match_t match, const char **tags,
    MPI_Comm comm);

/** Clean up communicator
 *
 * frees any memory allocated by the communicator
 *
 * \param [in,out] communicator communicator to clean up
 */
void
bfam_communicator_free(bfam_communicator_t *communicator);

/** Start communication
 *
 * \param [in,out] communicator communicator
 */
void
bfam_communicator_start(bfam_communicator_t *communicator);

/** Finish communication
 *
 * \param [in,out] communicator communicator
 */
void
bfam_communicator_finish(bfam_communicator_t *communicator);

#endif
