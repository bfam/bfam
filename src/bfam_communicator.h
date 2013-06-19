#ifndef BFAM_COMMUNICATOR_H
#define BFAM_COMMUNICATOR_H

#include <bfam_domain.h>

/**
 * structure for storing subdomain specific data for the communicator
 */
typedef struct bfam_comm_subdata
{
  bfam_subdomain_t* subdomain; /**< pointer to my local subdomain */

  size_t send_sz;  /**< amount of data this subdomain can send */
  void*  send_buf; /**< pointer for local send buffer */

  size_t recv_sz;  /**< amount of data this subdomain should receive */
  void*  recv_buf; /**< pointer for local recv buffer */
} bfam_comm_subdata_t;

/**
 * structure for storing processor specific data for the communicator
 */
typedef struct bfam_comm_procdata
{
  bfam_locidx_t rank; /**< neighboring processors rank */

  size_t send_sz;  /**< amount to send */
  void*  send_buf; /**< pointer to send buffer */

  size_t recv_sz;  /**< amount to recv */
  void*  recv_buf; /**< pointer to recv buffer */
} bfam_comm_procdata_t;

/**
 * structure for doing communication
 */
typedef struct bfam_communicator
{
  MPI_Comm comm;          /**< communicator used for communication */
  int tag;                /**< user specified tag for this communicator */
  bfam_locidx_t num_subs; /**< number of subdomains in the communicator */

  bfam_locidx_t  num_procs; /**< number of processors in the communicator */

  MPI_Request *send_request; /**< send request */
  MPI_Request *recv_request; /**< recv request */

  MPI_Status  *send_status; /**< send status */
  MPI_Status  *recv_status; /**< recv status */

  void*   send_buf; /**< full send buffer */
  void*   recv_buf; /**< full recv buffer */

  bfam_comm_procdata_t* proc_data; /**< array of structure with neighboring
                                        processor data */
  bfam_comm_subdata_t*  sub_data;  /**< array of structure with subdomains
                                        specific information */
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
