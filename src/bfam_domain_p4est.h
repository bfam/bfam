#ifndef BFAM_DOMAIN_P4EST_H
#define BFAM_DOMAIN_P4EST_H

#include <bfam_base.h>
#include <bfam_domain.h>

/**
 * structure containing a domain managed by p4est
 */
typedef struct bfam_domain_p4est
{
  bfam_domain_t         d;     /** parent domain */
  p4est_connectivity_t *conn;  /** connectivity for p4est */
  p4est_t              *p4est; /** forest of quadtrees */
} bfam_domain_p4est_t;


/* Domain managed by p4est based functions */
/** create a p4est managed domain
 *
 * \warning It is the callers responsibility to ensure that
 *          \a domComm and \a conn are freed after this domain is.
 *
 * \param [in]  domComm   pointer to the communicator for the domain
 * \param [in]  conn      pointer to the p4est connectivity for the domain
 *
 * \return the newly created p4est managed domain
 */
bfam_domain_p4est_t*
bfam_domain_p4est_new(MPI_Comm domComm,
                      p4est_connectivity_t *conn);

/** initializes a domain
 *
 * \param [in,out] domain  pointer to the p4est managed domain
 * \param [in]     domComm pointer to the communicator for the domain
 * \param [in]     conn    pointer to the p4est connectivity for the domain
 */
void
bfam_domain_p4est_init(bfam_domain_p4est_t *domain, MPI_Comm domComm,
                       p4est_connectivity_t *conn);

/** Clean up domain
 *
 * frees any memory allocated by the domain and calls free command on all
 * subdomains
 *
 * \param [in,out] domain domain to clean up
 */
void
bfam_domain_p4est_free(bfam_domain_p4est_t *domain);

/** Takes an initialized domain and generates a DG quad mesh
 *
 * \param [in,out] domain        pointer to the initialized p4est managed
 *                               domain
 * \param [in]     numSubdomains number of volume subdomains to generate
 * \param [in]     subdomainID   array of length \c p4est->local_num_quadrants
 *                               which indicates the subdomain id for each
 *                               element
 * \param [in]     N             array of orders for the volume subdomains
 */
void
bfam_domain_p4est_split_dgx_quad_subdomains(bfam_domain_p4est_t *domain,
    bfam_locidx_t numSubdomains, bfam_locidx_t *subdomainID, int *N);

#endif
