#ifndef BFAM_DOMAIN_PXEST_H
#define BFAM_DOMAIN_PXEST_H

#include <bfam_base.h>
#include <bfam_domain.h>

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wcast-align"
#include <p4est.h>
#include <p4est_bits.h>
#include <p4est_extended.h>
#include <p4est_ghost.h>
#include <p4est_lnodes.h>
#include <p4est_mesh.h>
#include <p4est_vtk.h>
#pragma clang diagnostic pop

/**
 * structure containing a domain managed by p4est
 */
typedef struct bfam_domain_pxest_2
{
  bfam_domain_t         base;  /** parent domain */
  p4est_connectivity_t *conn;  /** connectivity for p4est */
  p4est_t              *pxest; /** forest of quadtrees */
} bfam_domain_pxest_t_2;

/* Domain managed by pxest based functions */
/** create a pxest managed domain
 *
 * \warning It is the callers responsibility to ensure that
 *          \a domComm and \a conn are freed after this domain is.
 *
 * \param [in]  domComm   pointer to the communicator for the domain
 * \param [in]  conn      pointer to the pxest connectivity for the domain
 *
 * \return the newly created pxest managed domain
 */
bfam_domain_pxest_t_2*
bfam_domain_pxest_new_2(MPI_Comm domComm, p4est_connectivity_t *conn);

/** initializes a domain
 *
 * \param [in,out] domain  pointer to the pxest managed domain
 * \param [in]     domComm pointer to the communicator for the domain
 * \param [in]     conn    pointer to the pxest connectivity for the domain
 */
void
bfam_domain_pxest_init_2(bfam_domain_pxest_t_2 *domain, MPI_Comm domComm,
                         p4est_connectivity_t *conn);

/** Clean up domain
 *
 * frees any memory allocated by the domain and calls free command on all
 * subdomains
 *
 * \param [in,out] domain domain to clean up
 */
void
bfam_domain_pxest_free_2(bfam_domain_pxest_t_2 *domain);

/** Takes an initialized domain and generates a DG hex mesh
 *
 * \param [in,out] domain        pointer to the initialized pxest managed
 *                               domain
 * \param [in]     numSubdomains number of volume subdomains to generate
 * \param [in]     subdomainID   array of length \c pxest->local_num_quadrants
 *                               which indicates the subdomain id for each
 *                               element
 * \param [in]     N             array of orders for the volume subdomains
 */
void
bfam_domain_pxest_split_dgx_subdomains_2(bfam_domain_pxest_t_2 *domain,
    bfam_locidx_t numSubdomains, bfam_locidx_t *subdomainID, int *N);

#endif
