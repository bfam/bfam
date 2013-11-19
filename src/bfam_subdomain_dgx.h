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

  int            dim;          /* dimensionality of this dgx */

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


  bfam_real_t     *Dr;       /* 1D LGL differentiation matrix */

  bfam_long_real_t *V;       /* 1D Vandermonde matrix for this N */

  bfam_locidx_t    K;        /* Number of elements in the subdomain */

  bfam_locidx_t   *vmapM;    /* Mapping into the volume for the minus side of
                                the face mesh */
  bfam_locidx_t   *vmapP;    /* Mapping into the volume for the plus  side of
                                the face mesh */

  int            **fmask;    /* face mask */

  /* glue quantities
   * all glue quantities are postfixed with _m for minus side or _p for plus
   * side. All _p quantities are arrays to allow for multiple plus side
   * subdomains
   */
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

  bfam_locidx_t    *EToEp_m; /* Element     number on neighboring glue */
  bfam_locidx_t    *EToEm_m; /* Element     number on local subdomain */
  int8_t           *EToFm_m; /* Face        number on local subdomain */
  int8_t           *EToHm_m; /* Hanging     number on local subdomain */
  int8_t           *EToOm_m; /* Orientation number on local subdomain */
} bfam_subdomain_dgx_t;


#endif
