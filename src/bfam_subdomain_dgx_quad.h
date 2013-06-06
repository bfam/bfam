#ifndef BFAM_SUBDOMAIN_DGX_QUAD_H
#define BFAM_SUBDOMAIN_DGX_QUAD_H

#include <bfam_base.h>
#include <bfam_subdomain.h>

/*
 * Just like in p4est the quads are stored in z-order so the
 * face and corner orderings are:
 *
 *            2           3
 *            +-----------+
 *            |     3     |
 *            |           |
 *            | 0       1 |
 *            |           |
 *            |     2     |
 *            +-----------+
 *            0           1
 *
 * This subdomain is composed of nodal (tensor product LGL points) quads.
 * The naming convention follows that of:
 *
 * @book{Hesthaven:2007:NDG:1557392,
 *   author = {Hesthaven, Jan S. and Warburton, Tim},
 *   title = {Nodal Discontinuous {Galerkin} Methods: Algorithms, Analysis, and
 *            Applications},
 *   year = {2007},
 *   isbn = {0387720650, 9780387720654},
 *   edition = {1st},
 *   publisher = {Springer Publishing Company, Incorporated},
 * }
 *
 * where it makes sense.
 *
 */

typedef struct bfam_subdomain_dgx_quad
{
  bfam_subdomain_t base;
  int              N;   /* 1D Polynomial Order */
  int              Np;  /* Number of points in the element */
  int              Nfp; /* Number of points on a face of the element */
  int              Nfaces; /* Number of faces in each element */
  int              Ncorners; /* Number of corners in each element */
  bfam_real_t     *r;   /* 1D LGL Nodal Point in [-1,1] */
  bfam_real_t     *w;   /* 1D LGL Weights */

  bfam_locidx_t    K;   /* Number of elements in the subdomain */

  bfam_real_t     *x;   /* x-coordinates of the nodes */
  bfam_real_t     *y;   /* y-coordinates of the nodes */
  bfam_real_t     *z;   /* z-coordinates of the nodes */
} bfam_subdomain_dgx_quad_t;

/** create a dg quad subdomain.
 *
 * \param [in] name name of this subdomain
 * \param [in] N    polynomial order of elements in each dimension
 * \param [in] Nv   number of vertices in the subdomain
 * \param [in] VX   array of x-coordinates for the vertices
 * \param [in] VY   array of y-coordinates for the vertices
 * \param [in] VZ   array of z-coordinates for the vertices
 * \param [in] K    number of elements in the subdomain
 * \param [in] EToV Mapping such that \c EToV[k*4+c] gives the vertex number
 *                  for corner \c c of element \c k.
 * \param [in] EToE Mapping such that \c EToV[k*4+c] gives the vertex number
 *                  for corner \c c of element \c k.
 *
 * \return Initialized dg quad subdomain
 *
 */
bfam_subdomain_dgx_quad_t*
bfam_subdomain_dgx_quad_new(const char             *name,
                            const int               N,
                            const bfam_locidx_t     Nv,
                            const bfam_long_real_t *VX,
                            const bfam_long_real_t *VY,
                            const bfam_long_real_t *VZ,
                            const bfam_locidx_t     K,
                            const bfam_locidx_t    *EToV,
                            const bfam_locidx_t    *EToE,
                            const int8_t           *EToF);

/** initializes a dg quad subdomain
 *
 * \param [in,out] subdomain pointer to the subdomain to initialize
 * \param [in]     name      name of this subdomain
 * \param [in]     N         polynomial order of elements in each dimension
 * \param [in]     Nv        number of vertices in the subdomain
 * \param [in]     VX        array of x-coordinates for the vertices
 * \param [in]     VY        array of y-coordinates for the vertices
 * \param [in]     VZ        array of z-coordinates for the vertices
 * \param [in]     K         number of elements in the subdomain
 * \param [in]     EToV      Mapping such that \c EToV[k*4+c] gives the vertex
 *                           number for corner \c c of element \c k.
 * \param [in]     EToE      Mapping such that \c EToV[k*4+c] gives the vertex
 *                           number for corner \c c of element \c k.
 *
 */
void
bfam_subdomain_dgx_quad_init(bfam_subdomain_dgx_quad_t *subdomain,
                             const char                *name,
                             const int                  N,
                             const bfam_locidx_t        Nv,
                             const bfam_long_real_t    *VX,
                             const bfam_long_real_t    *VY,
                             const bfam_long_real_t    *VZ,
                             const bfam_locidx_t        K,
                             const bfam_locidx_t       *EToV,
                             const bfam_locidx_t       *EToE,
                             const int8_t              *EToF);

/** free up the memory allocated by the subdomain
 *
 * \param [in,out] subdomain subdomain to clean up
 */
void
bfam_subdomain_dgx_quad_free(bfam_subdomain_t *subdomain);

#endif
