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
#include <p4est_iterate.h>
#include <p4est_nodes.h>
#include <p4est_lnodes.h>
#include <p4est_mesh.h>
#include <p4est_vtk.h>
#pragma clang diagnostic pop

/**
 * structure containing a domain managed by p4est
 */
typedef struct bfam_domain_pxest_2
{
  bfam_domain_t base;         /** parent domain */
  bfam_dictionary_t *N2N;     /** Dictionary of projection operators */
  p4est_connectivity_t *conn; /** connectivity for p4est */
  p4est_t *pxest;             /** forest of quadtrees */
  bfam_dictionary_t *dgx_ops; /** Dictionary of dgx operators operators */
} bfam_domain_pxest_t_2;

typedef struct
{
  uint8_t *dst_to_adapt_flags;
  int8_t *dst_to_dst_chld_id;
  bfam_locidx_t *dst_to_src_subd_id;
  bfam_locidx_t *dst_to_src_elem_id;
  bfam_locidx_t num_dst;

  int8_t *coarse_dst_to_src_chld_id;
  bfam_locidx_t *coarse_dst_to_src_subd_id;
  bfam_locidx_t *coarse_dst_to_src_elem_id;
  bfam_locidx_t num_coarse_dst;
} bfam_domain_pxest_transfer_maps_t;

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
bfam_domain_pxest_t_2 *bfam_domain_pxest_new_2(MPI_Comm domComm,
                                               p4est_connectivity_t *conn);

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
 * \return the newly created pxest managed domain
 */
bfam_domain_pxest_t_2 *bfam_domain_pxest_new_ext_2(MPI_Comm domComm,
                                                   p4est_connectivity_t *conn,
                                                   p4est_locidx_t min_quadrants,
                                                   int min_level,
                                                   int fill_uniform);

/** initializes a domain
 *
 * \param [in,out] domain  pointer to the pxest managed domain
 * \param [in]     domComm pointer to the communicator for the domain
 * \param [in]     conn    pointer to the pxest connectivity for the domain
 */
void bfam_domain_pxest_init_2(bfam_domain_pxest_t_2 *domain, MPI_Comm domComm,
                              p4est_connectivity_t *conn);

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
void bfam_domain_pxest_init_ext_2(bfam_domain_pxest_t_2 *domain,
                                  MPI_Comm domComm, p4est_connectivity_t *conn,
                                  p4est_locidx_t min_quadrants, int min_level,
                                  int fill_uniform);

/** Clean up domain
 *
 * frees any memory allocated by the domain and calls free command on all
 * subdomains
 *
 * \param [in,out] domain domain to clean up
 */
void bfam_domain_pxest_free_2(bfam_domain_pxest_t_2 *domain);

/** Fill a \c glueID based on tree ids.
 *
 * This fills a \c glueID array for the quadrants based on glue ids given for
 * the trees.
 *
 * \param [in]  pxest        pointer to a p4est
 * \param [in]  tree_to_glue an array indicating the glue ids for faces in the
 *                           connectivity structure of the p4est mesh.  The ids
 *                           for the tree faces are stored in -x +x -y +y order
 *                           for each tree with the face index running the
 *                           fastest.
 * \param [out] quad_to_glue an array indicating the glue ids for faces in the
 *                            p4est mesh.  The ids for the quadrant faces are
 *                            stored in -x +x -y +y order for each quadrant
 *                            with the face index running the fastest.
 */
void bfam_domain_pxest_quad_to_glueid_2(p4est_t *pxest,
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
 * \param [in] glue_order        user callback function to allow the user to
 *                               set the order of the glue grids
 * \param [in] go_user_args      user argument for glue_order
 */
void bfam_domain_pxest_split_dgx_subdomains_2(
    bfam_domain_pxest_t_2 *domain, bfam_locidx_t numSubdomains,
    bfam_locidx_t *subdomainID, bfam_locidx_t *roots, int *N,
    bfam_locidx_t *glueID, bfam_glue_order_t glue_order, void *go_user_args);

/** Given a domain and a refined p4est generates a subdomain spliting
 *
 * \param [in]     domain        pointer to the initialized pxest managed
 *                               domain to split (if null the order of the
 *                               elements will be taken from pxest)
 * \param [in]     pxest         pointer to a p4est
 * \param [in]     pflags        flags on what type of p adaptivity
 * \param [out]    numSubdomains number of volume subdomains to generate
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
 */
void bfam_domain_pxest_compute_split_2(bfam_domain_pxest_t_2 *domain,
                                       p4est_t *pxest, uint8_t pflags,
                                       bfam_locidx_t *num_subdomains,
                                       bfam_locidx_t **subdomain_id,
                                       bfam_locidx_t **roots, int **N,
                                       bfam_locidx_t **glue_id);

/** Mark pxest quadrants with refinement info
 *
 * \param [in,out] domain        pointer to the initialized pxest managed
 *                               domain
 */
void bfam_domain_pxest_mark_elements_2(bfam_domain_pxest_t_2 *domain);

/** Takes an initialized domain and generates a DG hex mesh
 *
 * \param [in,out] domain        pointer to the pxest managed domain to create a
 *                               grid for
 * \param [in] nodes_transform   user callback function to allow the user to
 *                               further transform the nodal locations
 * \param [in] nt_user_args      user argument for nodes_transform
 *
 */

void bfam_domain_pxest_create_mesh_2(bfam_domain_pxest_t_2 *domain,
                                     bfam_dgx_nodes_transform_t nodes_transform,
                                     void *nt_user_args);

/** Generate mesh transfer maps
 *
 * \param [in,out] maps     pointer to the transfer maps to be filled
 * \param [in] domain_dst   destination domain
 * \param [in] domain_src   source domain
 *
 */
void
bfam_domain_pxest_transfer_maps_init_2(bfam_domain_pxest_transfer_maps_t *maps,
                                       bfam_domain_pxest_t_2 *domain_dst,
                                       bfam_domain_pxest_t_2 *domain_src);
/** Generate mesh transfer maps
 *
 * \param [in,out] maps     pointer to the transfer maps to be freed
 */
void
bfam_domain_pxest_transfer_maps_free_2(bfam_domain_pxest_transfer_maps_t *maps);

/** Adapt the mesh.
 *
 * This will adapt the mesh based on \c hadapt and \c padapt in the \c
 * bfam_subdomain_dgx structure.
 *
 * \param [in,out] domain        pointer to the initialized pxest managed
 *                               domain
 * \param [in] glue_order        user callback function to allow the user to
 *                               set the order of the glue grids
 * \param [in] go_user_args      user argument for glue_order
 */
void bfam_domain_pxest_adapt_2(bfam_domain_pxest_t_2 *domain,
                               bfam_glue_order_t glue_order,
                               void *go_user_args);

/** Callback for p4est user data init function.
 */
void bfam_domain_pxest_init_callback_2(p4est_t *p4est,
                                       p4est_topidx_t which_tree,
                                       p4est_quadrant_t *quadrant);

/* Callbacks for pxest quadrants */
int bfam_domain_pxest_quadrant_coarsen_2(p4est_t *p4est,
                                         p4est_topidx_t which_tree,
                                         p4est_quadrant_t *quadrants[]);
int bfam_domain_pxest_quadrant_refine_2(p4est_t *p4est,
                                        p4est_topidx_t which_tree,
                                        p4est_quadrant_t *quadrant);
void bfam_domain_pxest_quadrant_init_2(p4est_t *p4est,
                                       p4est_topidx_t which_tree,
                                       p4est_quadrant_t *quadrant);
void bfam_domain_pxest_quadrant_replace_2(p4est_t *p4est,
                                          p4est_topidx_t which_tree,
                                          int num_outgoing,
                                          p4est_quadrant_t *outgoing[],
                                          int num_incoming,
                                          p4est_quadrant_t *incoming[]);
#endif
