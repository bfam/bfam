#ifndef BFAM_DOMAIN_PXEST_H
#define BFAM_DOMAIN_PXEST_H

#include <bfam_base.h>
#include <bfam_domain.h>

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wcast-align"
#include <p4est_to_p8est.h>
#include <p8est.h>
#include <p8est_bits.h>
#include <p8est_extended.h>
#include <p8est_ghost.h>
#include <p8est_iterate.h>
#include <p8est_nodes.h>
#include <p8est_lnodes.h>
#include <p8est_mesh.h>
#include <p8est_vtk.h>
#pragma clang diagnostic pop

/**
 * structure containing a domain managed by p4est
 */
typedef struct bfam_domain_pxest_3
{
  bfam_domain_t base;         /** parent domain */
  bfam_dictionary_t *N2N;     /** Dictionary of projection operators */
  p8est_connectivity_t *conn; /** connectivity for p8est */
  p8est_t *pxest;             /** forest of quadtrees */
} bfam_domain_pxest_t_3;

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
bfam_domain_pxest_t_3 *bfam_domain_pxest_new_3(MPI_Comm domComm,
                                               p8est_connectivity_t *conn);

/** create a pxest managed domain
 *
 * \warning It is the callers responsibility to ensure that
 *          \a domComm and \a conn are freed after this domain is.
 *
 * \param [in]  domComm         pointer to the communicator for the domain
 * \param [in]  conn            pointer to the pxest connectivity for the domain
 * \param [in]  min_quadrants   Minimum initial quadrants per processor.
 *                              Makes the refinement pattern mpisize-specific.
 * \param [in]  min_level       The forest is refined at least to this level.
 *                              May be negative or 0, then it has no effect.
 * \param [in]  fill_uniform    If true, fill the forest with a uniform mesh
 *                              instead of the coarsest possible one.
 *                              The latter is partition-specific so that
 *                              is usually not a good idea.
*
 * \return the newly created pxest managed domain
 */
bfam_domain_pxest_t_3 *bfam_domain_pxest_new_ext_3(MPI_Comm domComm,
                                                   p8est_connectivity_t *conn,
                                                   p4est_locidx_t min_quadrants,
                                                   int min_level,
                                                   int fill_uniform);

/** initializes a domain
 *
 * \param [in,out] domain  pointer to the pxest managed domain
 * \param [in]     domComm pointer to the communicator for the domain
 * \param [in]     conn    pointer to the pxest connectivity for the domain
 */
void bfam_domain_pxest_init_3(bfam_domain_pxest_t_3 *domain, MPI_Comm domComm,
                              p8est_connectivity_t *conn);

/** initializes a domain
 *
 * \param [in,out] domain          pointer to the pxest managed domain
 * \param [in]     domComm         pointer to the communicator for the domain
 * \param [in]     conn            pointer to the pxest connectivity for the
 *                                 domain
 * \param [in]     min_quadrants   Minimum initial quadrants per processor.
 *                                 Makes the refinement pattern
 *                                 mpisize-specific.
 * \param [in]     min_level       The forest is refined at least to this level.
 *                                 May be negative or 0, then it has no effect.
 * \param [in]     fill_uniform    If true, fill the forest with a uniform mesh
 *                                 instead of the coarsest possible one.
 *                                 The latter is partition-specific so that
 *                                 is usually not a good idea.
 */
void bfam_domain_pxest_init_ext_3(bfam_domain_pxest_t_3 *domain,
                                  MPI_Comm domComm, p8est_connectivity_t *conn,
                                  p4est_locidx_t min_quadrants, int min_level,
                                  int fill_uniform);

/** Clean up domain
 *
 * frees any memory allocated by the domain and calls free command on all
 * subdomains
 *
 * \param [in,out] domain domain to clean up
 */
void bfam_domain_pxest_free_3(bfam_domain_pxest_t_3 *domain);

/** Fill a \c glueID based on tree ids.
 *
 * This fills a \c glueID array for the quadrants based on glue ids given for
 * the trees.
 *
 * \param [in]  pxest        pointer to a p8est
 * \param [in]  tree_to_glue an array indicating the glue ids for faces in the
 *                           connectivity structure of the p4est mesh.  The ids
 *                           for the tree faces are stored in -x +x -y +y -z +z
 *                           order for each tree with the face index running
 *                           the fastest.
 * \param [out] quad_to_glue an array indicating the glue ids for faces in the
 *                            p4est mesh.  The ids for the quadrant faces are
 *                            stored in -x +x -y +y -z +z order for each
 *                            quadrant with the face index running the fastest.
 */
void bfam_domain_pxest_quad_to_glueid_3(p4est_t *pxest,
                                        const bfam_locidx_t *tree_to_glueid,
                                        bfam_locidx_t *quad_to_glueid);

/** Takes an initialized domain and generates a DG hex mesh
 *
 * \param [in,out] domain        pointer to the initialized pxest managed
 *                               domain
 * \param [in]     numSubdomains number of volume subdomains to generate
 * \param [in]     subdomainID   array of length \c pxest->local_num_quadrants
 *                               which indicates the subdomain id for each
 *                               element
 * \param [in]     roots         array of roots for the volume subdomains
 * \param [in]     N             array of orders for the volume subdomains
 * \param [in]     glueID        array indicating what glue subdomains exist.
 *                               If the number is negative than no glue
 *                               subdomain will be created and if the number
 *                               is positive a subdomain will be created.
 *                               This is of length \c
 *                               pxest->local_num_quadrants*NumberOfFaces
 *                               if \c NULL it will be ignored.
 * \param [in] nodes_transform   user callback function to allow the user to
 *                               further transform the nodal locations
 * \param [in] user_args         user argument for nodes_transform
 */
void bfam_domain_pxest_split_dgx_subdomains_3(
    bfam_domain_pxest_t_3 *domain, bfam_locidx_t numSubdomains,
    bfam_locidx_t *subdomainID, bfam_locidx_t *roots, int *N,
    bfam_locidx_t *glueID, bfam_dgx_nodes_transform_t nodes_transform,
    void *user_args);

/** Adapt the mesh.
 *
 * This will adapt the mesh based on \c hadapt and \c padapt in the \c
 * bfam_subdomain_dgx structure.
 *
 * \param [in,out] domain        pointer to the initialized pxest managed
 *                               domain
 * \param [in] nodes_transform   user callback function to allow the user to
 *                               further transform the nodal locations
 * \param [in] user_args         user argument for nodes_transform
 */
void bfam_domain_pxest_adapt_3(bfam_domain_pxest_t_3 *domain,
                               bfam_dgx_nodes_transform_t nodes_transform,
                               void *user_args);

/** Callback for p4est user data init function.
 */
void bfam_domain_pxest_init_callback_3(p4est_t *p4est,
                                       p4est_topidx_t which_tree,
                                       p4est_quadrant_t *quadrant);

#endif
