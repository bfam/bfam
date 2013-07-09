#ifndef BFAM_SUBDOMAIN_SBP_H
#define BFAM_SUBDOMAIN_SBP_H

#include <bfam_base.h>
#include <bfam_subdomain.h>

/*
 * Just like in p4est the are stored in z-order so the
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
 */



/* Nb = buf_sz
 *
 *  (0,N[1])                                                        (N[0],N[1])
 *    +------------------------------------------------------------------+
 *    |                                                                  |
 *    |                                                                  |
 *    |  {0,Nl[1]+Nb[2]+Nb[3]}    {Nl[0]+Nb[0]+Nb[1]},Nl[1]+Nb[2]+Nb[3]} |
 *    |    +------------------------------------------------+            |
 *    |    |                                                |            |
 *    |    |  (gx[0],gx[1]+Nl[1]) (gx[0]+Nl[0],gx[1]+Nl[0]) |            |
 *    |    |  {Nb[0],Nl[1]+Nb[2]} {Nl[0]+Nb[0],Nl[1]+Nb[2]} |            |
 *    |    |        +-------------------------+             |            |
 *    |    |        |                         |             |            |
 *    |    |        |                         |             |            |
 *    |    |        |                         |             |            |
 *    |    |        +-------------------------+             |            |
 *    |    |  {Nb[0],Nb[2]}       {Nl[0]+Nb[0],Nl[1]}       |            |
 *    |    |  (gx[0],gx[1])       (gx[0]+Nl[0],gx[1])       |            |
 *    |    |                                                |            |
 *    |    +------------------------------------------------+            |
 *    |  {0,0}                           {Nl[0]+Nb[0]+Nb[1]},0}          |
 *    |                                                                  |
 *    +------------------------------------------------------------------+
 *  (0,0)                                                           (N[0],0)
 */
typedef struct bfam_subdomain_sbp
{
  bfam_subdomain_t base;
  int              p;   /**< 1D SBP finite difference order */
  int            dim;   /**< number of dimensions */
  bfam_locidx_t* sub_ix;  /**< local subdomain indices */
  bfam_locidx_t* sub_N;   /**< number pieces this subdomain is split in each dim */
  bfam_locidx_t  loc_id;
  bfam_locidx_t  num_id;

  bfam_gloidx_t   *N;   /**< "global" subdomain grid size (size N +1) */

  bfam_locidx_t   *Nl;  /**< "local"  subdomain grid size (size Nl+1) */
  bfam_locidx_t   *buf_sz;  /**< buffer size in each dimension (+/-) */
  /* total grid size for malloc is Prod_{i}(Nl[i] + 1 + buf_sz[2*i] + buf_sz[2*i+1]) */

  bfam_gloidx_t   *gx;  /**< my global starting indices
                             (corresponds to my data not ghost cells) */
} bfam_subdomain_sbp_t;

/** create a sbp subdomain.
 *
 * \param [in] id    unique id number for this subdomain
 * \param [in] sub_ix local id in each dimension
 * \param [in] sub_N number of pieces in each dimension
 * \param [in] name  name of this subdomain
 * \param [in] dim   dimension of the problem
 * \param [in] N     global grid top indices (size dim)
 * \param [in] Nl    local grid top indices (size dim)
 * \param [in] buf_sz    number of points in buffer for each edge (size 2*dim)
 * \param [in] gx    global data index (size dim)
 * \param [in] c_x   pointer to corners x coordinates in z-order (size 2^dim)
 * \param [in] c_y   pointer to corners y coordinates in z-order (size 2^dim)
 * \param [in] c_z   pointer to corners z coordinates in z-order (size 2^dim)
 *
 * \note z pointer can be null, in which case z is set to zero
 *
 * \return Initialized sbp subdomain
 *
 */
bfam_subdomain_sbp_t*
bfam_subdomain_sbp_new(const bfam_locidx_t     id,
                            const bfam_locidx_t *sub_ix,
                            const bfam_locidx_t *sub_N,
                            const char             *name,
                            const int               dim,
                            const bfam_gloidx_t    *N,
                            const bfam_locidx_t    *Nl,
                            const bfam_locidx_t    *buf_sz,
                            const bfam_gloidx_t    *gx,
                            const bfam_long_real_t *c_x,
                            const bfam_long_real_t *c_y,
                            const bfam_long_real_t *c_z);

/** initializes a sbp subdomain
 *
 * \param [in,out]   subdomain pointer to the subdomain to initialize
 * \param [in] id    unique id number for this subdomain
 * \param [in] sub_ix local id in each dimension
 * \param [in] sub_N number of pieces in each dimension
 * \param [in] name  name of this subdomain
 * \param [in] dim   dimension of the problem
 * \param [in] N     global grid top indices (size dim)
 * \param [in] Nl    local grid top indices (size dim)
 * \param [in] buf_sz    number of points in buffer for each edge (size 2*dim)
 * \param [in] gx    global data index (size dim)
 * \param [in] c_x   pointer to corners x coordinates in z-order (size 2^dim)
 * \param [in] c_y   pointer to corners y coordinates in z-order (size 2^dim)
 * \param [in] c_z   pointer to corners z coordinates in z-order (size 2^dim)
 *
 * \note z pointer can be null, in which case z is set to zero
 *
 */
void
bfam_subdomain_sbp_init(bfam_subdomain_sbp_t *subdomain,
                            const bfam_locidx_t     id,
                            const bfam_locidx_t* sub_ix,
                            const bfam_locidx_t* sub_N,
                            const char             *name,
                            const int               dim,
                            const bfam_gloidx_t    *N,
                            const bfam_locidx_t    *Nl,
                            const bfam_locidx_t    *buf_sz,
                            const bfam_gloidx_t    *gx,
                            const bfam_long_real_t *c_x,
                            const bfam_long_real_t *c_y,
                            const bfam_long_real_t *c_z);

/** free up the memory allocated by the subdomain
 *
 * \param [in,out] subdomain subdomain to clean up
 */
void
bfam_subdomain_sbp_free(bfam_subdomain_t *subdomain);

/**
 * intra glue grid
 */
typedef struct bfam_subdomain_sbp_intra_glue
{
  bfam_subdomain_t base;

  bfam_locidx_t     rank_m; /* Rank of the subdomain on the minus side */
  bfam_locidx_t     rank_p; /* Rank of the subdomain on the plus  side;
                             * if set to -1 then the subdomain is not
                             * connected on the plus side (i.e., it is
                             * a boundary subdomain).
                             */

  bfam_locidx_t     id_m;   /* Subdomain ID on the minus side */
  bfam_locidx_t     id_p;   /* Subdomain ID on the plus  side */

  bfam_locidx_t     loc_m;   /* local subdomain ID on the minus side */
  bfam_locidx_t     loc_p;   /* local subdomain ID on the plus side */

  bfam_locidx_t     fce_m;   /* face / corner / edge ID on the minus side */
  bfam_locidx_t     fce_p;   /* face / corner / edge ID on the plus side */

  bfam_subdomain_sbp_t *sub_m;  /* Local neighboring subdomain */

  int               face;  /* face of the subdomain being handled */

  bfam_locidx_t     *send_ix;   /* idices to be sent in Z-order */
  bfam_locidx_t     *recv_ix;   /* idices to be recv in Z-order */

  size_t            field_size_m;
  size_t            field_size_p;
} bfam_subdomain_sbp_intra_glue_t;

/** create a sbp intra glue subdomain.
 *
 * \param [in]     id          unique id number for this subdomain
 * \param [in]     name        name of this subdomain
 * \param [in]     rank_m      minus side processor rank
 * \param [in]     rank_p      plus  side processor rank
 * \param [in]     subdomain_m minus side subdomain pointer
 * \param [in]     face        face being handled
 *
 * \return Initialized dg quad glue subdomain
 *
 */
bfam_subdomain_sbp_intra_glue_t*
bfam_subdomain_sbp_intra_glue_new(const bfam_locidx_t              id,
                                 const char                      *name,
                                 const bfam_locidx_t              rank_m,
                                 const bfam_locidx_t              rank_p,
                                 bfam_subdomain_sbp_t            *sub_m,
                                 const int                        face);

/** create a sbp intra glue subdomain.
 *
 * \param [out]    subdomain    pointer to the subdomain to initialize
 *
 */
void
bfam_subdomain_sbp_intra_glue_init(bfam_subdomain_sbp_intra_glue_t  *subdomain,
                                   const bfam_locidx_t              id,
                                   const char                       *name,
                                   const bfam_locidx_t              rank_m,
                                   const bfam_locidx_t              rank_p,
                                   bfam_subdomain_sbp_t            *sub_m,
                                   const int                         face);

/** free up the memory allocated by the subdomain
 *
 * \param [in,out] subdomain subdomain to clean up
 */
void
bfam_subdomain_sbp_intra_glue_free(bfam_subdomain_t *subdomain);

/**
 * inter glue grid
 */
typedef struct bfam_subdomain_sbp_inter_glue
{
  bfam_subdomain_t base;

  bfam_locidx_t     rank_m; /* Rank of the subdomain on the minus side */
  bfam_locidx_t     rank_p; /* Rank of the subdomain on the plus  side;
                             * if set to -1 then the subdomain is not
                             * connected on the plus side (i.e., it is
                             * a boundary subdomain).
                             */

  bfam_locidx_t     id_m;   /* Subdomain ID on the minus side */
  bfam_locidx_t     id_p;   /* Subdomain ID on the plus  side */

  bfam_locidx_t     loc_m;   /* local subdomain ID on the minus side */
  bfam_locidx_t     loc_p;   /* local subdomain ID on the plus side */

  bfam_locidx_t     fce_m;   /* face / corner / edge ID on the minus side */
  bfam_locidx_t     fce_p;   /* face / corner / edge ID on the plus side */

  bfam_subdomain_sbp_t *sub_m;  /* Local neighboring subdomain */

  int               face;   /* face of the subdomain being handled */

  bfam_locidx_t     *ix;    /* face index I work with */
  int               orient; /* orientation of my neighbor */

  size_t            field_size_m;
  size_t            field_size_p;
} bfam_subdomain_sbp_inter_glue_t;

/** create a sbp inter glue subdomain.
 *
 * \param [in]     id          unique id number for this subdomain
 * \param [in]     name        name of this subdomain
 * \param [in]     rank_m      minus side processor rank
 * \param [in]     rank_p      plus  side processor rank
 * \param [in]     subdomain_m minus side subdomain pointer
 * \param [in]     ix          global face indices being handled by this glue
 * \param [in]     face        face being handled
 * \param [in]     orient      orientation code for neigh (see bfam_base.h)
 * \param [in]     id_p        plus sides id
 * \param [in]     loc_p       plus sides local id
 * \param [in]     face_p      plus sides face
 *
 * \return Initialized dg quad glue subdomain
 *
 */
bfam_subdomain_sbp_inter_glue_t*
bfam_subdomain_sbp_inter_glue_new(const bfam_locidx_t              id,
                                 const char                      *name,
                                 const bfam_locidx_t              rank_m,
                                 const bfam_locidx_t              rank_p,
                                 bfam_subdomain_sbp_t            *sub_m,
                                 const bfam_gloidx_t             *ix,
                                 const int                        face,
                                 const int                        orient,
                                 bfam_locidx_t                    id_p,
                                 bfam_locidx_t                    loc_p,
                                 bfam_locidx_t                    face_p);

/** create a sbp inter glue subdomain.
 *
 * \param [out]    subdomain    pointer to the subdomain to initialize
 * \param [in]     id           unique id number for this subdomain
 * \param [in]     name         name of this subdomain
 * \param [in]     rank_m       minus side processor rank
 * \param [in]     rank_p       plus  side processor rank
 * \param [in]     subdomain_m  minus side subdomain pointer
 * \param [in]     ix           global face indices being handled by this glue
 * \param [in]     face         face being handled
 * \param [in]     orient       orientation code for neigh (see bfam_base.h)
 * \param [in]     id_p        plus sides id
 * \param [in]     loc_p       plus sides local id
 * \param [in]     face_p      plus sides face
 *
 */
void
bfam_subdomain_sbp_inter_glue_init(bfam_subdomain_sbp_inter_glue_t  *subdomain,
                                   const bfam_locidx_t              id,
                                   const char                       *name,
                                   const bfam_locidx_t              rank_m,
                                   const bfam_locidx_t              rank_p,
                                   bfam_subdomain_sbp_t            *sub_m,
                                   const bfam_gloidx_t             *ix,
                                   const int                        face,
                                   const int                        orient,
                                   bfam_locidx_t                    id_p,
                                   bfam_locidx_t                    loc_p,
                                   bfam_locidx_t                    face_p);

/** free up the memory allocated by the subdomain
 *
 * \param [in,out] subdomain subdomain to clean up
 */
void
bfam_subdomain_sbp_inter_glue_free(bfam_subdomain_t *subdomain);

#endif
