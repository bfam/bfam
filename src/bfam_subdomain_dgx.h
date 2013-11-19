#ifndef BFAM_SUBDOMAIN_DGX_H
#define BFAM_SUBDOMAIN_DGX_H

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
 * This subdomain is composed of nodal (tensor product LGL points)
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

struct bfam_subdomain_dgx;

typedef struct bfam_subdomain_dgx
{
  bfam_subdomain_t base;

  int              dim;        /* dimensionality of this dgx */

  int              N;          /* 1D Polynomial Order */
  int              Np;         /* Number of points in the element */
  int              Nfp;        /* Number of points on a face of the element */
  int             *Ngeo;       /* geometry based quantities:
                                * corners, faces, edges, etc.
                                */
  int              Nh;         /* Number of interpolations to glue */
  int              No;         /* Number of orientations */

  bfam_real_t     *r;          /* 1D LGL Nodal Point in [-1,1] */
  bfam_real_t     *w;          /* 1D LGL Weights */
  bfam_real_t     *wi;         /* inverse of 1D LGL Weights */
  bfam_real_t     *exact_mass; /* exact mass matrix for this grid */


  bfam_real_t      *Dr;      /* 1D LGL differentiation matrix */

  bfam_long_real_t *V;       /* 1D Vandermonde matrix for this N */

  bfam_locidx_t     K;       /* Number of elements in the subdomain */

  bfam_locidx_t    *vmapM;   /* Mapping into the volume for the minus side of
                                the face mesh */
  bfam_locidx_t    *vmapP;   /* Mapping into the volume for the plus  side of
                                the face mesh */

  int             **fmask;   /* face mask */

  /* glue quantities
   * all glue quantities are postfixed with _m for minus side or _p for plus
   * side. All _p quantities are arrays to allow for multiple plus side
   * subdomains
   */
  int             num_p;    /* number of plus sides */
  int               N_m;    /* 1D Polynomial Order on the minus side */
  int              *N_p;    /* 1D Polynomial Order on the plus  side */

  bfam_locidx_t     rank_m; /* Rank of the subdomain on the minus side */
  bfam_locidx_t    *rank_p; /* Rank of the subdomain on the plus  side;
                             * if set to -1 then the subdomain is not
                             * connected on the plus side (i.e., it is
                             * a boundary subdomain).
                             */

  bfam_locidx_t     id_m;   /* Sort Id of the subdomain on the minus side */
  bfam_locidx_t    *id_p;   /* Sort Id of the subdomain on the plus  side */

  bfam_locidx_t     s_m;    /* Id of the subdomain on the minus side */
  bfam_locidx_t    *s_p;    /* Id of the subdomain on the plus  side */

  struct bfam_subdomain_dgx *sub_m;  /* Local neighboring subdomain */

  int                      Nh_m;    /* number of interpolation operators from
                                     * minus side
                                     */

  bfam_real_t     **interpolation_m; /* array of interpolation operators;
                                      * the first is for non-hanging faces
                                      * the rest are for the hanging faces;
                                      * if the operator is NULL it is assumed
                                      * to be the identity operator.
                                      */

  bfam_real_t     **projection_m;   /* array of projection operators;
                                     * the first is for non-hanging faces
                                     * the rest are for the hanging faces;
                                     * if the operator is NULL it is assumed
                                     * to be the identity operator.
                                     */

  bfam_real_t     **massprojection_m; /* array of mass projection operators;
                                       * the first is for non-hanging faces
                                       * the rest are for the hanging faces;
                                       * if the operator is NULL it is assumed
                                       * to be the identity operator.
                                       */

  bfam_locidx_t   **EToE_p; /* Element     number on neighboring glue */
  bfam_locidx_t    *EToE_m; /* Element     number on local subdomain */
  int8_t           *EToF_m; /* Face        number on local subdomain */
  int8_t           *EToH_m; /* Hanging     number on local subdomain */
  int8_t           *EToO_m; /* Orientation number on local subdomain */
} bfam_subdomain_dgx_t;

/** create a dgx subdomain.
 *
 * \param [in] id     unique id number for this subdomain
 * \param [in] name   name of this subdomain
 * \param [in] N      polynomial order of elements in each dimension
 * \param [in] Nv     number of vertices in the subdomain
 * \param [in] num_Vi number of coordinate for vertices to store
 * \param [in] Vi     array of array of coordinates for the vertives
 * \param [in] K      number of elements in the subdomain
 * \param [in] EToV   Mapping such that \c EToV[k*4+c] gives the vertex number
 *                    for corner \c c of element \c k.
 * \param [in] EToE   Mapping such that \c EToE[k*4+f] gives the element number
 *                    connected to face \c f of element \c k.
 * \param [in] EToF   Mapping such that \c EToF[k*4+f] gives the
 *                    face+orientation number connected to face \c f of element
 *                    \c k.
 * \param [in] dim    number of (computational) dimensions
 *
 * \return Initialized dgx subdomain
 *
 */

bfam_subdomain_dgx_t*
bfam_subdomain_dgx_new_(const bfam_locidx_t      id,
                       const char              *name,
                       const int                N,
                       const bfam_locidx_t      Nv,
                       const int                num_Vi,
                       const bfam_long_real_t **Vi,
                       const bfam_locidx_t      K,
                       const bfam_locidx_t     *EToV,
                       const bfam_locidx_t     *EToE,
                       const int8_t            *EToF,
                       const int                dim);

/** initializes a dgx subdomain
 *
 * \param [in,out] subdomain pointer to the subdomain to initialize
 * \param [in]     id        unique id number for this subdomain
 * \param [in]     name      name of this subdomain
 * \param [in]     N         polynomial order of elements in each dimension
 * \param [in]     Nv        number of vertices in the subdomain
 * \param [in]     num_Vi    number of coordinate for vertices to store
 * \param [in]     Vi        array of array of coordinates for the vertives
 * \param [in]     K         number of elements in the subdomain
 * \param [in]     EToV      Mapping such that \c EToV[k*4+c] gives the vertex
 *                           number for corner \c c of element \c k.
 * \param [in]     EToE      Mapping such that \c EToE[k*4+f] gives the element
 *                           number connected to face \c f of element \c k.
 * \param [in]     EToF      Mapping such that \c EToF[k*4+f] gives the
 *                           face+orientation number connected to face \c f of
 *                           element \c k.
 * \param [in]     dim       number of (computational) dimensions
 */
void
bfam_subdomain_dgx_init_(bfam_subdomain_dgx_t *subdomain,
                        const bfam_locidx_t        id,
                        const char                *name,
                        const int                  N,
                        const bfam_locidx_t        Nv,
                        const int                  num_Vi,
                        const bfam_long_real_t   **Vi,
                        const bfam_locidx_t        K,
                        const bfam_locidx_t       *EToV,
                        const bfam_locidx_t       *EToE,
                        const int8_t              *EToF,
                        const int                  dim);

#define X(dg_dim) \
bfam_subdomain_dgx_t*                                         \
bfam_subdomain_dgx_new_##dg_dim(const bfam_locidx_t      id,  \
                             const char              *name,   \
                             const int                N,      \
                             const bfam_locidx_t      Nv,     \
                             const int                num_Vi, \
                             const bfam_long_real_t **Vi,     \
                             const bfam_locidx_t      K,      \
                             const bfam_locidx_t     *EToV,   \
                             const bfam_locidx_t     *EToE,   \
                             const int8_t            *EToF,   \
                             const int                dim);   \
void                                                              \
bfam_subdomain_dgx_init_##dg_dim(bfam_subdomain_dgx_t *subdomain, \
                        const bfam_locidx_t        id,            \
                        const char                *name,          \
                        const int                  N,             \
                        const bfam_locidx_t        Nv,            \
                        const int                  num_Vi,        \
                        const bfam_long_real_t   **Vi,            \
                        const bfam_locidx_t        K,             \
                        const bfam_locidx_t       *EToV,          \
                        const bfam_locidx_t       *EToE,          \
                        const int8_t              *EToF,          \
                        const int                  dim);
BFAM_LIST_OF_BGX_DIMENSIONS
#undef X

#endif
