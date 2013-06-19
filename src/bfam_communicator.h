#ifndef BFAM_COMMUNICATOR_H
#define BFAM_COMMUNICATOR_H

#include <bfam_domain.h>

typedef struct bfam_comm_subdata
{
  bfam_subdomain_t* subdomain;

  size_t send_sz;
  void*  send_buf;

  size_t recv_sz;
  void*  recv_buf;
} bfam_comm_subdata_t;

typedef struct bfam_comm_procdata
{
  bfam_locidx_t rank;

  size_t send_sz;
  void*  send_buf;

  size_t recv_sz;
  void*  recv_buf;
} bfam_comm_procdata_t;

typedef struct bfam_communicator
{
  MPI_Comm comm;
  bfam_locidx_t num_subs;

  bfam_locidx_t  num_procs;

  MPI_Request *send_request;
  MPI_Request *recv_request;

  MPI_Status  *send_status;
  MPI_Status  *recv_status;

  void*   send_buf;
  void*   recv_buf;

  bfam_comm_procdata_t* proc_data;
  bfam_comm_subdata_t*  sub_data;

  int tag;
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
    const char **tags, MPI_Comm comm, int tag);

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
    MPI_Comm comm, int tag);

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
