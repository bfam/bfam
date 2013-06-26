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
  int              Nh;  /* Number of interpolations to glue */
  int              No;  /* Number of orientations */
  bfam_real_t     *r;   /* 1D LGL Nodal Point in [-1,1] */
  bfam_real_t     *w;   /* 1D LGL Weights */

  bfam_real_t     *Dr;  /* 1D LGL differentiation matrix */

  bfam_locidx_t    K;   /* Number of elements in the subdomain */

  int            **fmask; /* face mask */
} bfam_subdomain_dgx_quad_t;

/** create a dg quad subdomain.
 *
 * \param [in] id   unique id number for this subdomain
 * \param [in] name name of this subdomain
 * \param [in] N    polynomial order of elements in each dimension
 * \param [in] Nv   number of vertices in the subdomain
 * \param [in] VX   array of x-coordinates for the vertices
 * \param [in] VY   array of y-coordinates for the vertices
 * \param [in] VZ   array of z-coordinates for the vertices
 * \param [in] K    number of elements in the subdomain
 * \param [in] EToV Mapping such that \c EToV[k*4+c] gives the vertex number
 *                  for corner \c c of element \c k.
 * \param [in] EToE Mapping such that \c EToE[k*4+f] gives the element number
 *                  connected to face \c f of element \c k.
 * \param [in] EToF Mapping such that \c EToF[k*4+f] gives the face+orientation
 *                  number connected to face \c f of element \c k.
 *
 * \return Initialized dg quad subdomain
 *
 */
bfam_subdomain_dgx_quad_t*
bfam_subdomain_dgx_quad_new(const bfam_locidx_t     id,
                            const char             *name,
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
 * \param [in]     id        unique id number for this subdomain
 * \param [in]     name      name of this subdomain
 * \param [in]     N         polynomial order of elements in each dimension
 * \param [in]     Nv        number of vertices in the subdomain
 * \param [in]     VX        array of x-coordinates for the vertices
 * \param [in]     VY        array of y-coordinates for the vertices
 * \param [in]     VZ        array of z-coordinates for the vertices
 * \param [in]     K         number of elements in the subdomain
 * \param [in]     EToV      Mapping such that \c EToV[k*4+c] gives the vertex
 *                           number for corner \c c of element \c k.
 * \param [in]     EToE      Mapping such that \c EToE[k*4+f] gives the element
 *                           number connected to face \c f of element \c k.
 * \param [in]     EToF      Mapping such that \c EToF[k*4+f] gives the
 *                           face+orientation number connected to face \c f of
 *                           element \c k.
 */
void
bfam_subdomain_dgx_quad_init(bfam_subdomain_dgx_quad_t *subdomain,
                             const bfam_locidx_t        id,
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

typedef struct bfam_subdomain_dgx_quad_glue
{
  bfam_subdomain_t base;

  int               N_m;    /* 1D Polynomial Order on the minus side */
  int               N_p;    /* 1D Polynomial Order on the plus  side */

  bfam_locidx_t     rank_m; /* Rank of the subdomain on the minus side */
  bfam_locidx_t     rank_p; /* Rank of the subdomain on the plus  side;
                             * if set to -1 then the subdomain is not
                             * connected on the plus side (i.e., it is
                             * a boundary subdomain).
                             */

  bfam_locidx_t     id_m;   /* Sort Id of the subdomain on the minus side */
  bfam_locidx_t     id_p;   /* Sort Id of the subdomain on the plus  side */

  bfam_locidx_t     s_m;    /* Id of the subdomain on the minus side */
  bfam_locidx_t     s_p;    /* Id of the subdomain on the plus  side */

  bfam_subdomain_dgx_quad_t *sub_m;  /* Local neighboring subdomain */

  int               N;        /* 1D Polynomial Order of the glue grid */
  int               Np;       /* Number of points in the element */
  int               Nfp;      /* Number of points on a face of the element */
  int               Nfaces;   /* Number of faces in each element */
  int               Ncorners; /* Number of corners in each element */

  bfam_real_t      *r;   /* 1D LGL Nodal Point in [-1,1] */
  bfam_real_t      *w;   /* 1D LGL Weights */

  bfam_locidx_t     K;   /* Number of elements in the subdomain */

  int               num_interp; /* number of interpolation operators */

  bfam_real_t     **interpolation; /* array of interpolation operators;
                                    * the first is for non-hanging faces
                                    * the rest are for the hanging faces;
                                    * if the operator is NULL it is assumed
                                    * to be the identity operator.
                                    */
  bfam_locidx_t    *EToEp; /* Element     number on neighboring glue */
  bfam_locidx_t    *EToEm; /* Element     number on local subdomain */
  int8_t           *EToFm; /* Face        number on local subdomain */
  int8_t           *EToHm; /* Hanging     number on local subdomain */
  int8_t           *EToOm; /* Orientation number on local subdomain */
} bfam_subdomain_dgx_quad_glue_t;

/** create a dg quad glue subdomain.
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
 *
 * \return Initialized dg quad glue subdomain
 *
 */
bfam_subdomain_dgx_quad_glue_t*
bfam_subdomain_dgx_quad_glue_new(const bfam_locidx_t              id,
                                 const char                      *name,
                                 const int                        N_m,
                                 const int                        N_p,
                                 const bfam_locidx_t              rank_m,
                                 const bfam_locidx_t              rank_p,
                                 const bfam_locidx_t              id_m,
                                 const bfam_locidx_t              id_p,
                                 bfam_subdomain_dgx_quad_t       *sub_m,
                                 bfam_locidx_t                   *ktok_m,
                                 const bfam_locidx_t              K,
                                 bfam_subdomain_face_map_entry_t *mapping);

/** initializes a dg quad glue subdomain.
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
 *
 */
void
bfam_subdomain_dgx_quad_glue_init(bfam_subdomain_dgx_quad_glue_t  *subdomain,
                                  const bfam_locidx_t              id,
                                  const char                      *name,
                                  const int                        N_m,
                                  const int                        N_p,
                                  const bfam_locidx_t              rank_m,
                                  const bfam_locidx_t              rank_p,
                                  const bfam_locidx_t              id_m,
                                  const bfam_locidx_t              id_p,
                                  bfam_subdomain_dgx_quad_t       *sub_m,
                                  bfam_locidx_t                   *ktok_m,
                                  const bfam_locidx_t              K,
                                  bfam_subdomain_face_map_entry_t *mapping);

/** free up the memory allocated by the subdomain
 *
 * \param [in,out] subdomain subdomain to clean up
 */
void
bfam_subdomain_dgx_quad_glue_free(bfam_subdomain_t *subdomain);

#endif
