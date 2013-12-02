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

typedef struct bfam_subdomain_dgx_glue_data
{
  bfam_subdomain_glue_data_t base;

  bfam_locidx_t    *EToEp; /* Element     number on connected subdomain */
  bfam_locidx_t    *EToEm; /* Element     number on local subdomain */
  int8_t           *EToFm; /* Face        number on local subdomain */
  int8_t           *EToHm; /* Hanging     number on local subdomain */
  int8_t           *EToOm; /* Orientation number on local subdomain */
  bfam_locidx_t   **mapOm; /* mapping the orientation */
  int          num_orient; /* number of orientations */

  /* The following pointers should only be \ne NULL on the minus side */
  int                  num_interp;   /* number of interpolation operators */

  bfam_real_t     **interpolation;   /* array of interpolation operators;
                                      * the first is for non-hanging faces
                                      * the rest are for the hanging faces;
                                      * if the operator is NULL it is assumed
                                      * to be the identity operator.
                                      */

  bfam_real_t     **projection;     /* array of projection operators;
                                     * the first is for non-hanging faces
                                     * the rest are for the hanging faces;
                                     * if the operator is NULL it is assumed
                                     * to be the identity operator.
                                     */

  bfam_real_t     **massprojection;   /* array of mass projection operators;
                                       * the first is for non-hanging faces
                                       * the rest are for the hanging faces;
                                       * if the operator is NULL it is assumed
                                       * to be the identity operator.
                                       */

  bfam_real_t     *exact_mass;         /* exact mass matrix for this grid */

} bfam_subdomain_dgx_glue_data_t;

typedef struct bfam_subdomain_dgx
{
  bfam_subdomain_t base;

  int              dim;        /* dimensionality of this dgx */

  int              N;          /* 1D Polynomial Order */
  int              Np;         /* Number of points in the element */
  int             *Ngp;        /* Number of geometry points */

  int              numg;       /* number of geometry types */
  int             *Ng;         /* geometry based quantities:
                                * faces, edges, corners, etc.
                                */
  //JK int              Nh;         /* Number of interpolations to glue */
  //JK int              No;         /* Number of orientations */

  bfam_real_t     *r;          /* 1D LGL Nodal Point in [-1,1] */
  bfam_real_t     *w;          /* 1D LGL Weights */
  bfam_real_t     *wi;         /* inverse of 1D LGL Weights */


  bfam_real_t      *Dr;      /* 1D LGL differentiation matrix */

  bfam_long_real_t *V;       /* 1D Vandermonde matrix for this N */

  bfam_locidx_t     K;       /* Number of elements in the subdomain */

  bfam_locidx_t    *vmapM;   /* Mapping into the volume for the minus side of
                                the face mesh */
  bfam_locidx_t    *vmapP;   /* Mapping into the volume for the plus  side of
                                the face mesh */

  int            ***gmask;   /* geometry mask: same order as Ng */

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

/** create a dgx glue subdomain.
 *
 * \param [in]     id          unique id number for this subdomain
 * \param [in]     name        name of this subdomain
 * \param [in]     N_m         minus side polynomial order of elements in each
 *                             dimension
 * \param [in]     N_p         plus  side polynomial order of elements in each
 *                             dimension
 * \param [in]     rank_m      minus side processor rank
 * \param [in]     rank_p      plus  side processor rank
 * \param [in]     id_m        minus side subdomain id
 * \param [in]     id_p        plus  side subdomain id
 * \param [in]     subdomain_m minus side subdomain pointer
 * \param [in]     ktok_m      map: element number -> minus side element number
 * \param [in]     K           number of elements in the glue grid
 * \param [in,out] mapping     face mapping (might get sorted)
 * \param [in]     inDIM        dimension of the subdomain
 *
 * \return Initialized dg quad glue subdomain
 *
 */
bfam_subdomain_dgx_t*
bfam_subdomain_dgx_glue_new_(const bfam_locidx_t              id,
                             const char                      *name,
                             const int                        N_m,
                             const int                        N_p,
                             const bfam_locidx_t              rank_m,
                             const bfam_locidx_t              rank_p,
                             const bfam_locidx_t              id_m,
                             const bfam_locidx_t              id_p,
                             bfam_subdomain_dgx_t            *sub_m,
                             bfam_locidx_t                   *ktok_m,
                             const bfam_locidx_t              K,
                             bfam_subdomain_face_map_entry_t *mapping,
                             const int                        inDIM);

/** initializes a dg glue subdomain.
 *
 * \param [out]    subdomain    pointer to the subdomain to initialize
 * \param [in]     id           unique id number for this subdomain
 * \param [in]     name         name of this subdomain
 * \param [in]     N_m          minus side polynomial order of elements in each
 *                              dimension
 * \param [in]     N_p          plus  side polynomial order of elements in each
 *                              dimension
 * \param [in]     rank_m       minus side processor rank
 * \param [in]     rank_p       plus  side processor rank
 * \param [in]     id_m         minus side subdomain id
 * \param [in]     id_p         plus  side subdomain id
 * \param [in]     subdomain_m  minus side subdomain pointer
 * \param [in]     ktok_m       map: element number -> minus side element number
 * \param [in]     K            number of elements in the glue grid
 * \param [in,out] mapping      face mapping (might get sorted)
 * \param [in]     inDIM        dimension of the subdomain
 *
 */
void
bfam_subdomain_dgx_glue_init_(bfam_subdomain_dgx_t  *subdomain,
                              const bfam_locidx_t              id,
                              const char                      *name,
                              const int                        N_m,
                              const int                        N_p,
                              const bfam_locidx_t              rank_m,
                              const bfam_locidx_t              rank_p,
                              const bfam_locidx_t              id_m,
                              const bfam_locidx_t              id_p,
                              bfam_subdomain_dgx_t            *sub_m,
                              bfam_locidx_t                   *ktok_m,
                              const bfam_locidx_t              K,
                              bfam_subdomain_face_map_entry_t *mapping,
                              const int                        inDIM);

/** free up the memory allocated by the subdomain
 *
 * \param [in,out] subdomain subdomain to clean up
 */
void
bfam_subdomain_dgx_free_(bfam_subdomain_t *subdomain);

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
                        const int                  dim);          \
void                                                           \
bfam_subdomain_dgx_free_##dg_dim(bfam_subdomain_t *subdomain); \
bfam_subdomain_dgx_t*                                                          \
bfam_subdomain_dgx_glue_new_##dg_dim(const bfam_locidx_t              id,      \
                                      const char                      *name,   \
                                      const int                        N_m,    \
                                      const int                        N_p,    \
                                      const bfam_locidx_t              rank_m, \
                                      const bfam_locidx_t              rank_p, \
                                      const bfam_locidx_t              id_m,   \
                                      const bfam_locidx_t              id_p,   \
                                      bfam_subdomain_dgx_t            *sub_m,  \
                                      bfam_locidx_t                   *ktok_m, \
                                      const bfam_locidx_t              K,      \
                                      bfam_subdomain_face_map_entry_t *mapping,\
                                      const int                        inDIM); \
void                                                                       \
bfam_subdomain_dgx_glue_init_##dg_dim(bfam_subdomain_dgx_t  *subdomain,    \
                              const bfam_locidx_t              id,         \
                              const char                      *name,       \
                              const int                        N_m,        \
                              const int                        N_p,        \
                              const bfam_locidx_t              rank_m,     \
                              const bfam_locidx_t              rank_p,     \
                              const bfam_locidx_t              id_m,       \
                              const bfam_locidx_t              id_p,       \
                              bfam_subdomain_dgx_t            *sub_m,      \
                              bfam_locidx_t                   *ktok_m,     \
                              const bfam_locidx_t              K,          \
                              bfam_subdomain_face_map_entry_t *mapping,    \
                              const int                        inDIM);
BFAM_LIST_OF_DGX_DIMENSIONS
#undef X

#endif
