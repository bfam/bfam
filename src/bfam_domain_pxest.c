#include <bfam_base.h>

#ifdef BFAM_PXEST_DIMENSION

#if   BFAM_PXEST_DIMENSION==2
#include <bfam_domain_pxest_2.h>
#elif BFAM_PXEST_DIMENSION==3
#include <bfam_domain_pxest_3.h>
#else
#error "Bad Dimension"
#endif

#include <bfam_log.h>
#include <bfam_subdomain_dgx.h>



#define DIM BFAM_PXEST_DIMENSION

#if   DIM==2
#define BFAM_PXEST_RELATIVE_ORIENTATIONS 2
#define BDIM 1
#define BFAM_PXEST_CONNECT P4EST_CONNECT_FULL
#define BFAM_PXEST_ORIENTATION(n,nf,o) (o)
#define BFAM_PXEST_BOHTONH(bo,h) (((bo)+(h))%(2))
#define BFAM_PXEST_BOHTONH_INV(bo,h) BFAM_PXEST_BOHTONH(bo,h)
#elif DIM==3
#define BFAM_PXEST_RELATIVE_ORIENTATIONS 4
#define BFAM_PXEST_ORIENTATION(n,nf,o) BFAM_P8EST_ORIENTATION(n,nf,o)
#define BDIM 2
#define BFAM_PXEST_CONNECT P8EST_CONNECT_FULL
#define BFAM_PXEST_BOHTONH(bo,h) (bfam_p8est_perm_to_order[(bo)][(h)])
#define BFAM_PXEST_BOHTONH_INV(bo,h) (bfam_p8est_perm_to_order_inv[(bo)][(h)])
#else
#error "Bad Dimension"
#endif

#define              bfam_domain_pxest_t \
  BFAM_APPEND_EXPAND(bfam_domain_pxest_t_,DIM)
#define              bfam_domain_pxest_new \
  BFAM_APPEND_EXPAND(bfam_domain_pxest_new_,DIM)
#define              bfam_domain_pxest_init \
  BFAM_APPEND_EXPAND(bfam_domain_pxest_init_,DIM)
#define              bfam_domain_pxest_free \
  BFAM_APPEND_EXPAND(bfam_domain_pxest_free_,DIM)
#define              bfam_subdomain_dgx_new \
  BFAM_APPEND_EXPAND(bfam_subdomain_dgx_new_,DIM)
#define              bfam_subdomain_dgx_glue_new \
  BFAM_APPEND_EXPAND(bfam_subdomain_dgx_glue_new_,BDIM)
#define              bfam_domain_pxest_quad_to_glueid \
  BFAM_APPEND_EXPAND(bfam_domain_pxest_quad_to_glueid_,DIM)
#define              bfam_domain_pxest_split_dgx_subdomains \
  BFAM_APPEND_EXPAND(bfam_domain_pxest_split_dgx_subdomains_,DIM)

bfam_domain_pxest_t*
bfam_domain_pxest_new(MPI_Comm domComm, p4est_connectivity_t *conn)
{
  bfam_domain_pxest_t* newDomain = bfam_malloc(sizeof(bfam_domain_pxest_t));
  bfam_domain_pxest_init(newDomain,domComm,conn);
  return newDomain;
}

void
bfam_domain_pxest_init(bfam_domain_pxest_t *domain, MPI_Comm domComm,
                       p4est_connectivity_t *conn)
{
  bfam_domain_init(&domain->base,domComm);

  domain->conn = conn;
  domain->pxest = p4est_new(domComm,conn,0,NULL,NULL);
}

void
bfam_domain_pxest_free(bfam_domain_pxest_t *domain)
{
  /* Memory we don't manage */
  domain->conn = NULL;

  /* Memory we do manage */
  p4est_destroy(domain->pxest);
  domain->pxest = NULL;

  bfam_domain_free(&domain->base);
}

static void
bfam_domain_pxest_dgx_print_stats(bfam_domain_pxest_t *domain)
{
  bfam_domain_t * dbase = &domain->base;
  bfam_subdomain_t **subdomains =
    bfam_malloc(dbase->numSubdomains*sizeof(bfam_subdomain_t**));

  bfam_locidx_t numSubdomains = 0;

  const char *volume[] = {"_volume", NULL};

  bfam_domain_get_subdomains(dbase, BFAM_DOMAIN_AND, volume,
      dbase->numSubdomains, subdomains, &numSubdomains);

  const size_t GRID_PTS = 0;
  const size_t ELEMENTS = 1;
  const size_t NUM_VALS = 2;

  bfam_gloidx_t vals_loc[NUM_VALS];
  bfam_gloidx_t vals_glo[NUM_VALS];

  for(size_t i = 0; i < NUM_VALS; ++i)
  {
    vals_loc[i] = 0;
    vals_glo[i] = 0;
  }

  for(bfam_locidx_t s = 0; s < numSubdomains; ++s)
  {
    bfam_subdomain_dgx_t *sub = (bfam_subdomain_dgx_t*) subdomains[s];

    vals_loc[GRID_PTS] += sub->K*sub->Np;
    vals_loc[ELEMENTS] += sub->K;
  }

  BFAM_MPI_CHECK(MPI_Reduce(vals_loc, vals_glo, NUM_VALS,
        BFAM_GLOIDX_MPI, MPI_SUM, 0, dbase->comm));

  BFAM_ROOT_INFO("Domain Stats --- Elements: %jd", (intmax_t) vals_glo[ELEMENTS]);
  BFAM_ROOT_INFO("Domain Stats --- Grid Pts: %jd", (intmax_t) vals_glo[GRID_PTS]);

  bfam_free(subdomains);
}

#ifdef BFAM_DEBUG
static void
bfam_domain_pxest_tree_to_glueid_check(p4est_connectivity_t *conn,
    const p4est_locidx_t* tree_to_glueid)
{
  for(p4est_topidx_t t = 0; t < conn->num_trees; ++t)
  {
    for(int f = 0; f < P4EST_FACES; ++f)
    {
      p4est_topidx_t nt =      conn->tree_to_tree[P4EST_FACES*t + f];
      int            nf = (int)conn->tree_to_face[P4EST_FACES*t + f];

      bfam_locidx_t  id = tree_to_glueid[P4EST_FACES* t+ f];
      bfam_locidx_t nid = tree_to_glueid[P4EST_FACES*nt+nf];

      BFAM_ABORT_IF_NOT(id == nid,
          "tree_to_glueid invalid for (%jd,%d)->%jd!=%jd<-(%jd,%d)",
          (intmax_t) t, f, (intmax_t) id, (intmax_t) nid, (intmax_t) nt, nf);
    }
  }
}
#endif

typedef struct bfam_domain_glueid_iter_data
{
  const bfam_locidx_t *tree_to_glueid;
  bfam_locidx_t *quad_to_glueid;
} bfam_domain_glueid_iter_data_t;

static void
set_quad_to_glueid(p4est_iter_face_info_t *info, void *arg)
{
  bfam_domain_glueid_iter_data_t *data = (bfam_domain_glueid_iter_data_t*) arg;

  int                 limit = (int) info->sides.elem_count;
  p4est_iter_face_side_t *fside;

  if(!info->tree_boundary)
    return;

  for(int i = 0; i < limit; ++i)
  {
    fside = p4est_iter_fside_array_index_int(&info->sides, i);

    int face = (int) fside->face;
    p4est_topidx_t treeid = fside->treeid;

    p4est_tree_t *tree = p4est_tree_array_index (info->p4est->trees, treeid);
    p4est_locidx_t offset = tree->quadrants_offset;


    BFAM_ASSERT(treeid >= 0 && treeid < info->p4est->connectivity->num_trees);

    if(!fside->is_hanging)
    {
      p4est_locidx_t qid = fside->is.full.quadid + offset;
      if(!fside->is.full.is_ghost)
      {
        BFAM_ASSERT(qid >= 0 && qid < info->p4est->local_num_quadrants);
        BFAM_LDEBUG("Glueid adding face (%jd %jd):(%jd %jd)",
            (intmax_t)treeid, (intmax_t)face, (intmax_t)qid, (intmax_t)face);
        data->quad_to_glueid[P4EST_FACES*qid+face] =
          data->tree_to_glueid[P4EST_FACES*treeid+face];
      }
    }
    else
    {
      for(int h = 0; h < P4EST_HALF; ++h)
      {
        p4est_locidx_t qid = fside->is.hanging.quadid[h] + offset;
        if(!fside->is.hanging.is_ghost[h])
        {
          BFAM_ASSERT(qid >= 0 && qid < info->p4est->local_num_quadrants);
          BFAM_LDEBUG("Glueid adding face (%jd %jd):(%jd %jd)",
            (intmax_t)treeid, (intmax_t)face, (intmax_t)qid, (intmax_t)face);
          data->quad_to_glueid[P4EST_FACES*qid+face] =
            data->tree_to_glueid[P4EST_FACES*treeid+face];
        }
      }
    }
  }

  return;
}


void
bfam_domain_pxest_quad_to_glueid(p4est_t *pxest,
    const bfam_locidx_t* tree_to_glueid, bfam_locidx_t* quad_to_glueid)
{
#ifdef BFAM_DEBUG
  bfam_domain_pxest_tree_to_glueid_check(pxest->connectivity, tree_to_glueid);
#endif

  for(p4est_locidx_t i = 0; i < pxest->local_num_quadrants*P4EST_FACES; ++i)
    quad_to_glueid[i] = -1;

  bfam_domain_glueid_iter_data_t data = {tree_to_glueid,quad_to_glueid};
  p4est_iterate(pxest, NULL, &data, NULL, set_quad_to_glueid,
#if DIM==3
      NULL,
#endif
      NULL);
}

static bfam_locidx_t
bfam_domain_pxest_num_parallel_faces(p4est_mesh_t *mesh)
{
  const p4est_locidx_t K = mesh->local_num_quadrants;

  bfam_locidx_t numParallelFaces = 0;

  for(p4est_locidx_t k = 0; k < K; ++k)
  {
    for (int f = 0; f < P4EST_FACES; ++f)
    {
      const p4est_locidx_t ck = mesh->quad_to_quad[P4EST_FACES * k + f];
      const int            cf = mesh->quad_to_face[P4EST_FACES * k + f];

      if(cf >= 0)
      {
        /*
         * neighbor is same or double size
         */
        if(ck >= mesh->local_num_quadrants)
          ++numParallelFaces;
      }
      else
      {
        /*
         * two neighbors half size
         */
        p4est_locidx_t *cks;
        cks = sc_array_index(mesh->quad_to_half, ck);

        for(int h = 0; h < P4EST_HALF; ++h)
          if(cks[h] >= mesh->local_num_quadrants)
            ++numParallelFaces;
      }
    }
  }

  BFAM_LDEBUG("Counted %jd parallel faces.", (intmax_t)numParallelFaces);

  return numParallelFaces;
}


/** Build the parallel face mapping array.
 *
 * This array is used to determine the order in which data is sent
 * and received around the forest.  The order can be obtained by sorting
 * the mapping using \c bfam_subdomain_face_send_cmp
 * and \c bfam_domain_pxest_face_recv_cmp comparison
 * functions.
 *
 * \param [in]  mesh             p4est mesh to build the mapping for.
 * \param [in]  glueID           glue id of each face in the mesh
 * \param [in]  numParallelFaces the number of parallel faces in the
 *                               p4est mesh.
 * \param [out] mapping          the mapping array that will be filled.
 *
 */
static void
bfam_domain_pxest_parallel_face_mapping(p4est_ghost_t *ghost, p4est_mesh_t *mesh,
    bfam_locidx_t *glueID,
    bfam_locidx_t numParallelFaces, bfam_subdomain_face_map_entry_t *mapping)
{
  const p4est_locidx_t K = mesh->local_num_quadrants;

  for(p4est_locidx_t k = 0, sk = 0; k < K; ++k)
  {
    for (int f = 0; f < P4EST_FACES; ++f)
    {
      const p4est_locidx_t ck = mesh->quad_to_quad[P4EST_FACES * k + f];
      const int            cf = mesh->quad_to_face[P4EST_FACES * k + f];

      const bfam_locidx_t  glueid = (glueID) ? glueID[P4EST_FACES * k + f] : -1;

      if(cf >= 0)
      {
        /*
         * neighbor is same or double size
         */
        if(ck >= mesh->local_num_quadrants)
        {
          p4est_locidx_t ghostid = ck-mesh->local_num_quadrants;
          p4est_quadrant_t *ghostquad = p4est_quadrant_array_index (&ghost->ghosts, (size_t) ghostid);
          p4est_locidx_t ghostp  = mesh->ghost_to_proc[ghostid];
          p4est_locidx_t ghostk  = ghostquad->p.piggy3.local_num;
          int8_t         ghostf = cf;
          int8_t         ghosth = 0;

          if(ghostf >= BFAM_PXEST_RELATIVE_ORIENTATIONS*P4EST_FACES)
          {
            ghostf -= BFAM_PXEST_RELATIVE_ORIENTATIONS*P4EST_FACES;

            ghosth = ghostf/(BFAM_PXEST_RELATIVE_ORIENTATIONS*P4EST_FACES) + 1;
            ghostf = ghostf%(BFAM_PXEST_RELATIVE_ORIENTATIONS*P4EST_FACES);

            BFAM_ASSERT(ghosth >= 0 && ghosth <= P4EST_HALF);
          }

          int8_t o = ghostf/P4EST_FACES;

          ghostf   = ghostf%P4EST_FACES;

          const int8_t bo = BFAM_PXEST_ORIENTATION(f,ghostf,o);
          if(ghosth > 0)
            ghosth = BFAM_PXEST_BOHTONH(bo,ghosth-1)+1;

          mapping[sk].np = ghostp;
          mapping[sk].ns = 0;
          mapping[sk].nk = ghostk;
          mapping[sk].nf = ghostf;
          mapping[sk].nh = ghosth;
          mapping[sk].gi = ghostid;
          mapping[sk].i  = sk;
          mapping[sk].s  = 0;
          mapping[sk].k  = k;
          mapping[sk].f  = f;
          mapping[sk].h  = 0;
          mapping[sk].o  = bo;
          mapping[sk].id = glueid;
          ++sk;
        }
      }
      else
      {
        /*
         * two neighbors half size
         */
        p4est_locidx_t *cks;
        cks = sc_array_index(mesh->quad_to_half, ck);

        for(int8_t h = 0; h < P4EST_HALF; ++h)
        {
          if(cks[h] >= mesh->local_num_quadrants)
          {
            p4est_locidx_t ghostid = cks[h]-mesh->local_num_quadrants;
            p4est_quadrant_t *ghostquad = p4est_quadrant_array_index (&ghost->ghosts, (size_t) ghostid);
            p4est_locidx_t ghostp  = mesh->ghost_to_proc[ghostid];
            p4est_locidx_t ghostk  = ghostquad->p.piggy3.local_num;
            int8_t         ghostf  =
              ((BFAM_PXEST_RELATIVE_ORIENTATIONS*P4EST_FACES) + cf)%
              P4EST_FACES;
            int8_t         ghosth  = 0;
            int8_t         o       =
              ((BFAM_PXEST_RELATIVE_ORIENTATIONS*P4EST_FACES) + cf)/
              P4EST_FACES;

            const int8_t bo = BFAM_PXEST_ORIENTATION(f,ghostf,o);

            mapping[sk].np = ghostp;
            mapping[sk].ns = 0;
            mapping[sk].nk = ghostk;
            mapping[sk].nf = ghostf;
            mapping[sk].nh = ghosth;
            mapping[sk].gi = ghostid;
            mapping[sk].i  = sk;
            mapping[sk].s  = 0;
            mapping[sk].k  = k;
            mapping[sk].f  = f;
            mapping[sk].h  = BFAM_PXEST_BOHTONH_INV(bo,h) + 1;
            mapping[sk].o  = bo;
            mapping[sk].id = glueid;
            ++sk;
          }
        }
      }
    }
  }

  qsort(mapping, numParallelFaces, sizeof(bfam_subdomain_face_map_entry_t),
      bfam_subdomain_face_send_cmp);

#ifdef BFAM_DEBUG
  {
    BFAM_LDEBUG("mapping: %5s %5s %5s %5s %5s %5s %5s %5s %5s %5s %5s",
        "np", "nk", "nf", "nh", "gi", "i", "k", "f", "h", "o", "id");
    for(bfam_locidx_t i = 0; i < numParallelFaces; ++i)
    {
      BFAM_LDEBUG(
          "mapping: %5jd %5jd %5jd %5jd %5jd %5jd %5jd %5jd %5jd %5jd %5jd",
          (intmax_t)mapping[i].np,
          (intmax_t)mapping[i].nk,
          (intmax_t)mapping[i].nf,
          (intmax_t)mapping[i].nh,
          (intmax_t)mapping[i].gi,
          (intmax_t)mapping[i].i,
          (intmax_t)mapping[i].k,
          (intmax_t)mapping[i].f,
          (intmax_t)mapping[i].h,
          (intmax_t)mapping[i].o,
          (intmax_t)mapping[i].id);
    }
  }
#endif
}

/** Count the number of processor neighbors.
 *
 * \param [in] numParallelFaces number of parallel faces in the mapping
 * \param [in] mapping          parallel face mapping array
 *
 * \returns the number of processor neighbors.
 */
static bfam_locidx_t
bfam_domain_pxest_parallel_face_num_neighbors(bfam_locidx_t numParallelFaces,
    bfam_subdomain_face_map_entry_t *mapping)
{
  qsort(mapping, numParallelFaces, sizeof(bfam_subdomain_face_map_entry_t),
      bfam_subdomain_face_send_cmp);

  bfam_locidx_t numNeighbors = 0;

  if(numParallelFaces != 0)
  {
    numNeighbors = 1;
    for(bfam_locidx_t i = 1; i < numParallelFaces; ++i)
      if(mapping[i-1].np != mapping[i].np)
        ++numNeighbors;
  }

  return numNeighbors;
}

/** Determine the neighbor ranks and number of faces per neighbor.
 *
 * \param [in]  numParallelFaces number of parallel faces in the mapping
 * \param [in]  mapping          parallel face mapping array
 * \param [in]  numNeighbors     the number of processor neighbors
 * \param [out] numNeighborFaces number of faces per neighbor
 * \param [out] neighborRank     rank of the neighbor
 */
static void
bfam_domain_pxest_num_neighbor_faces(bfam_locidx_t numParallelFaces,
    bfam_subdomain_face_map_entry_t *mapping, bfam_locidx_t numNeighbors,
    bfam_locidx_t *numNeighborFaces, bfam_locidx_t *neighborRank)
{
  qsort(mapping, numParallelFaces, sizeof(bfam_subdomain_face_map_entry_t),
      bfam_subdomain_face_send_cmp);


  if(numParallelFaces != 0)
  {
    BFAM_ASSERT(numNeighbors > 0);

    numNeighborFaces[0] = 1;
    neighborRank[0] = mapping[0].np;

    for(bfam_locidx_t i = 1, sk = 0; i < numParallelFaces; ++i)
    {
      if(mapping[i-1].np == mapping[i].np)
      {
        ++numNeighborFaces[sk];
      }
      else
      {
        ++sk;
        numNeighborFaces[sk] = 1;
        neighborRank[sk] = mapping[i].np;
      }
    }
  }

#ifdef BFAM_DEBUG
  {
    /*
     * Sanity check on the number of parallel faces per processor
     */
    bfam_locidx_t faces = 0;
    for(bfam_locidx_t n = 0; n < numNeighbors; ++n)
      faces += numNeighborFaces[n];

    BFAM_LDEBUG(" XXX  NumNeighbors %jd", (intmax_t) numNeighbors);
    for(bfam_locidx_t n = 0; n < numNeighbors; ++n)
      BFAM_LDEBUG(" XXX     neighborFaces[%jd] = %jd", (intmax_t) n,
          (intmax_t) numNeighborFaces[n]);

    for(bfam_locidx_t n = 0; n < numNeighbors; ++n)
      BFAM_LDEBUG(" XXX     neighborRanks[%jd] = %jd", (intmax_t) n,
          (intmax_t) neighborRank[n]);


    BFAM_LDEBUG(" XXX  faces: %jd =?= %jd", (intmax_t) faces,
        (intmax_t) numParallelFaces);

    BFAM_ASSERT(faces == numParallelFaces);
  }
#endif
}

/** Send and receive \a entries values per parallel face.
 *
 * \param [in]  numParallelFaces number of parallel faces in the mapping
 * \param [in]  mapping          parallel face mapping array
 * \param [in]  numNeighbors     the number of processor neighbors
 * \param [in]  numNeighborFaces number of faces per neighbor
 * \param [in]  neighborRank     rank of the neighbor
 * \param [in]  entries          number of entries to send per face
 * \param [in]  sendBuffer       values to send
 * \param [out] recvBuffer       values received
 *
 */
static void
bfam_domain_pxest_parallel_face_send_recv(MPI_Comm comm,
    bfam_locidx_t numParallelFaces, bfam_subdomain_face_map_entry_t *mapping,
    bfam_locidx_t numNeighbors, bfam_locidx_t *numNeighborFaces,
    bfam_locidx_t *neighborRank, bfam_locidx_t entries,
    bfam_locidx_t *sendBuffer, bfam_locidx_t *recvBuffer)
{
  const int tag = 666;

  MPI_Request *sendRequest = bfam_malloc(2*numNeighbors*sizeof(MPI_Request));
  MPI_Request *recvRequest = sendRequest + numNeighbors;

  MPI_Status *status = bfam_malloc(2*numNeighbors*sizeof(MPI_Status));

  for(bfam_locidx_t n = 0; n < numNeighbors; ++n)
  {
    sendRequest[n] = MPI_REQUEST_NULL;
    recvRequest[n] = MPI_REQUEST_NULL;
  }

  /*
   * Post receives
   */
  for(bfam_locidx_t n = 0, sk = 0; n < numNeighbors; ++n)
  {
    const bfam_locidx_t count = entries*numNeighborFaces[n];
    BFAM_MPI_CHECK(MPI_Irecv(&recvBuffer[sk], count, BFAM_LOCIDX_MPI,
          neighborRank[n], tag, comm, &recvRequest[n]));
    sk += count;
  }

  /*
   * Post sends
   */
  for(bfam_locidx_t n = 0, sk = 0; n < numNeighbors; ++n)
  {
    const bfam_locidx_t count = entries*numNeighborFaces[n];
    BFAM_MPI_CHECK(MPI_Isend(&sendBuffer[sk], count, BFAM_LOCIDX_MPI,
          neighborRank[n], tag, comm, &sendRequest[n]));
    sk += count;
  }

  /*
   * Wait for the communication
   */
  BFAM_MPI_CHECK(MPI_Waitall(2*numNeighbors, sendRequest, status));

  bfam_free(status);
  bfam_free(sendRequest);
}


#ifdef BFAM_DEBUG
/** Debug function to test the mapping.
 *
 * \param [in] numParallelFaces number of parallel faces in the mapping
 * \param [in] mapping          parallel face mapping array
 * \param [in] numNeighbors     the number of processor neighbors
 * \param [in] numNeighborFaces number of faces per neighbor
 * \param [in] neighborRank     rank of the neighbor
 *
 * \note This function aborts if it found an error in the mapping
 *
 */
static void
bfam_domain_pxest_parallel_face_mapping_check(MPI_Comm comm,
    bfam_locidx_t numParallelFaces, bfam_subdomain_face_map_entry_t *mapping,
    bfam_locidx_t numNeighbors, bfam_locidx_t *numNeighborFaces,
    bfam_locidx_t *neighborRank)
{
  const int entries = 9;

  bfam_locidx_t *recvBuffer =
    bfam_malloc_aligned(numParallelFaces*entries*sizeof(bfam_locidx_t));
  bfam_locidx_t *sendBuffer =
    bfam_malloc_aligned(numParallelFaces*entries*sizeof(bfam_locidx_t));

  /*
   * Sort the mapping in send order
   */
  qsort(mapping, numParallelFaces, sizeof(bfam_subdomain_face_map_entry_t),
      bfam_subdomain_face_send_cmp);

  /*
   * Fill Send buffers
   */
  for(bfam_locidx_t i = 0; i < numParallelFaces; ++i)
  {
    sendBuffer[entries*i + 0] = mapping[i].ns;
    sendBuffer[entries*i + 1] = mapping[i].nk;
    sendBuffer[entries*i + 2] = mapping[i].nf;
    sendBuffer[entries*i + 3] = mapping[i].nh;

    sendBuffer[entries*i + 4] = mapping[i].s;
    sendBuffer[entries*i + 5] = mapping[i].k;
    sendBuffer[entries*i + 6] = mapping[i].f;
    sendBuffer[entries*i + 7] = mapping[i].h;

    sendBuffer[entries*i + 8] = mapping[i].id;
  }

  bfam_domain_pxest_parallel_face_send_recv(comm, numParallelFaces, mapping,
      numNeighbors, numNeighborFaces, neighborRank, entries, sendBuffer,
      recvBuffer);

  /*
   * Sort the mapping in recv order
   */
  qsort(mapping, numParallelFaces, sizeof(bfam_subdomain_face_map_entry_t),
      bfam_subdomain_face_recv_cmp);

  /*
   * Check the receive buffers
   */
  for(bfam_locidx_t i = 0; i < numParallelFaces; ++i)
  {
    BFAM_ASSERT(recvBuffer[entries*i + 0] == mapping[i].s);
    BFAM_ASSERT(recvBuffer[entries*i + 1] == mapping[i].k);
    BFAM_ASSERT(recvBuffer[entries*i + 2] == mapping[i].f);
    BFAM_ASSERT(recvBuffer[entries*i + 3] == mapping[i].h);

    BFAM_ASSERT(recvBuffer[entries*i + 4] == mapping[i].ns);
    BFAM_ASSERT(recvBuffer[entries*i + 5] == mapping[i].nk);
    BFAM_ASSERT(recvBuffer[entries*i + 6] == mapping[i].nf);
    BFAM_ASSERT(recvBuffer[entries*i + 7] == mapping[i].nh);

    BFAM_ASSERT(recvBuffer[entries*i + 8] == mapping[i].id);
  }

  bfam_free_aligned(sendBuffer);
  bfam_free_aligned(recvBuffer);
}
#endif

/** Get neighbor's subdomain information.
 *
 * \param [in]     comm             MPI Communicator
 * \param [in]     numParallelFaces number of parallel faces in the mapping
 * \param [in,out] mapping          parallel face mapping array (subdomain
 *                                  information is filled in the mapping).
 * \param [in]     numNeighbors     the number of processor neighbors
 * \param [in]     numNeighborFaces number of faces per neighbor
 * \param [in]     neighborRank     rank of the neighbor
 * \param [in]     subdomainID      array of subdomain ids for each local
 *                                  element
 * \param [in]     N                array of orders for each subdomain
 * \param [out]    ghostSubdomainID array of subdomain ids for each ghost
 *                                  element in the mesh
 * \param [out]     ghostN          array of orders for each ghost element in
 *                                  the mesh
 */
static void
bfam_domain_pxest_fill_ghost_subdomain_ids(MPI_Comm comm,
    bfam_locidx_t numParallelFaces, bfam_subdomain_face_map_entry_t *mapping,
    bfam_locidx_t numNeighbors, bfam_locidx_t *numNeighborFaces,
    bfam_locidx_t *neighborRank, bfam_locidx_t *subdomainID,
    int *N, bfam_locidx_t *ghostSubdomainID, bfam_locidx_t *ghostN)
{
  const int entries = 2;

  bfam_locidx_t *recvBuffer =
    bfam_malloc_aligned(numParallelFaces*entries*sizeof(bfam_locidx_t));
  bfam_locidx_t *sendBuffer =
    bfam_malloc_aligned(numParallelFaces*entries*sizeof(bfam_locidx_t));

  /*
   * Sort the mapping in send order
   */
  qsort(mapping, numParallelFaces, sizeof(bfam_subdomain_face_map_entry_t),
      bfam_subdomain_face_send_cmp);

  /*
   * Fill Send buffers
   */
  for(bfam_locidx_t i = 0; i < numParallelFaces; ++i)
  {
    bfam_locidx_t s = subdomainID[mapping[i].k];
    sendBuffer[entries*i + 0] = s;
    sendBuffer[entries*i + 1] = N[s];
  }

  bfam_domain_pxest_parallel_face_send_recv(comm, numParallelFaces, mapping,
      numNeighbors, numNeighborFaces, neighborRank, entries, sendBuffer,
      recvBuffer);

  /*
   * Sort the mapping in recv order
   */
  qsort(mapping, numParallelFaces, sizeof(bfam_subdomain_face_map_entry_t),
      bfam_subdomain_face_recv_cmp);

  /*
   * Fill mapping with subdomain id
   */
  for(bfam_locidx_t i = 0; i < numParallelFaces; ++i)
    mapping[i].s = subdomainID[mapping[i].k];

  /*
   * Fill the ghost information
   */
  for(bfam_locidx_t i = 0; i < numParallelFaces; ++i)
  {
                     mapping[i].ns  = recvBuffer[entries*i + 0];
    ghostSubdomainID[mapping[i].gi] = recvBuffer[entries*i + 0];
              ghostN[mapping[i].gi] = recvBuffer[entries*i + 1];
  }

  bfam_free_aligned(sendBuffer);
  bfam_free_aligned(recvBuffer);
}

/** Count the number of local inter-subdomain faces.
 *
 * Note that faces on nonconforming interfaces are always counted even if
 * they connect to the same subdomain.  This is because we will use a
 * glue grid to handle the nonconforming interfaces.
 *
 * \param [in]  mesh          p4est mesh of the elements
 * \param [in]  subdomainID   subdomain id of each element in the mesh
 * \param [in]  glueID        glue id of each face in the mesh
 *
 * \return number of total local inter-subdomain faces.
 *
 */
static bfam_locidx_t
bfam_domain_pxest_num_inter_subdomain_faces(p4est_mesh_t *mesh,
    bfam_locidx_t *subdomainID, bfam_locidx_t *glueID)
{
  const p4est_locidx_t K = mesh->local_num_quadrants;

  bfam_locidx_t numInterSubdomainFaces = 0;

  for(p4est_locidx_t k = 0; k < K; ++k)
  {
    const bfam_locidx_t  idk   = subdomainID[k];

    for (int f = 0; f < P4EST_FACES; ++f)
    {
      const p4est_locidx_t ck = mesh->quad_to_quad[P4EST_FACES * k + f];
      const int            cf = mesh->quad_to_face[P4EST_FACES * k + f];

      const bfam_locidx_t  glueid = (glueID) ? glueID[P4EST_FACES * k + f] : -1;

      if(cf >= 0)
      {
        /*
         * Neighbor is same or double size
         */
        if(ck < mesh->local_num_quadrants)
        {
          /*
           * Neighbor is on the same processor
           */
          const bfam_locidx_t  idnk   = subdomainID[ck];

          int hanging;

          if(cf >= BFAM_PXEST_RELATIVE_ORIENTATIONS*P4EST_FACES)
            hanging = 1;
          else
            hanging = 0;

          /*
           * Count intra and inter subdomain faces.
           *
           * Only count same subdomain to same subdomain if it is a hanging
           * face.
           */
          if((idnk == idk && hanging) || idnk != idk ||
              (glueid >= 0 && (ck != k || cf != f)))
            ++numInterSubdomainFaces;
        }
      }
      else
      {
        p4est_locidx_t *cks;
        cks = sc_array_index(mesh->quad_to_half, ck);
        for(int8_t h = 0; h < P4EST_HALF; ++h)
          if(cks[h] < mesh->local_num_quadrants)
            ++numInterSubdomainFaces;
      }
    }
  }

  return numInterSubdomainFaces;
}

/** Build the inter subdomain face mapping array.
 *
 * This array is used to determine the order in which data is sent
 * and received around the subdomains.  The order can be obtained by sorting
 * the mapping using \c bfam_subdomain_face_send_cmp
 * and \c bfam_subdomain_face_recv_cmp comparison
 * functions.
 *
 * \param [in]  rank                   local MPI rank.
 * \param [in]  mesh                   pxest mesh to build the mapping for.
 * \param [in]  subdomainID            subdomain id of each element in the mesh.
 * \param [in]  glueID                 glue id of each face in the mesh
 * \param [in]  numInterSubdomainFaces the number of inter subdomain faces in
 *                                     the pxest mesh.
 * \param [out] mapping                the mapping array that will be filled.
 *
 */
static void
bfam_domain_pxest_inter_subdomain_face_mapping(bfam_locidx_t rank,
    p4est_mesh_t *mesh, bfam_locidx_t *subdomainID,
    bfam_locidx_t *glueID,
    bfam_locidx_t numInterSubdomainFaces,
    bfam_subdomain_face_map_entry_t *mapping)
{
  const p4est_locidx_t K = mesh->local_num_quadrants;

  for(p4est_locidx_t k = 0, sk = 0; k < K; ++k)
  {
    const bfam_locidx_t  idk   = subdomainID[k];

    for (int f = 0; f < P4EST_FACES; ++f)
    {
      const p4est_locidx_t ck = mesh->quad_to_quad[P4EST_FACES * k + f];
      const int            cf = mesh->quad_to_face[P4EST_FACES * k + f];

      const bfam_locidx_t  glueid = (glueID) ? glueID[P4EST_FACES * k + f] : -1;

      if(cf >= 0)
      {
        /*
         * Neighbor is same or double size
         */
        if(ck < mesh->local_num_quadrants)
        {
          /*
           * Neighbor is on the same processor
           */
          const bfam_locidx_t  idnk   = subdomainID[ck];

          int hanging;

          if(cf >= BFAM_PXEST_RELATIVE_ORIENTATIONS*P4EST_FACES)
            hanging = 1;
          else
            hanging = 0;

          int8_t         nf = cf;
          int8_t         nh = 0;

          if(nf >= BFAM_PXEST_RELATIVE_ORIENTATIONS*P4EST_FACES)
          {
            nf -= BFAM_PXEST_RELATIVE_ORIENTATIONS*P4EST_FACES;

            nh = nf / (BFAM_PXEST_RELATIVE_ORIENTATIONS*P4EST_FACES) + 1;
            nf = nf % (BFAM_PXEST_RELATIVE_ORIENTATIONS*P4EST_FACES);
          }

          int8_t o = nf/P4EST_FACES;

          nf = nf%P4EST_FACES;
          const int8_t bo = BFAM_PXEST_ORIENTATION(f,nf,o);

          if(nh > 0)
          {
            nh = BFAM_PXEST_BOHTONH(bo,nh-1)+1;
          }

          /*
           * Count intra and inter subdomain faces.
           *
           * Only count same subdomain to same subdomain if it is a hanging
           * face.
           */
          if((idnk == idk && hanging) || idnk != idk ||
              (glueid >= 0 && (ck != k || cf != f)))
          {
            BFAM_ASSERT(sk < numInterSubdomainFaces);

            mapping[sk].np = rank;

            mapping[sk].ns = idnk;
            mapping[sk].nk = ck;
            mapping[sk].nf = nf;
            mapping[sk].nh = nh;

            mapping[sk].s  = idk;
            mapping[sk].k  = k;
            mapping[sk].f  = f;
            mapping[sk].h  = 0;
            mapping[sk].o  = bo;

            mapping[sk].i  = -1;
            mapping[sk].gi = -1;

            mapping[sk].id = glueid;
            ++sk;
          }
        }
      }
      else
      {
        p4est_locidx_t *cks;
        cks = sc_array_index(mesh->quad_to_half, ck);
        for(int8_t h = 0; h < P4EST_HALF; ++h)
        {
          if(cks[h] < mesh->local_num_quadrants)
          {
            BFAM_ASSERT(sk < numInterSubdomainFaces);

            const bfam_locidx_t  idnk   = subdomainID[cks[h]];

            int8_t nf =
              ((BFAM_PXEST_RELATIVE_ORIENTATIONS*P4EST_FACES) + cf)%P4EST_FACES;
            int8_t nh = 0;
            int8_t o  =
              ((BFAM_PXEST_RELATIVE_ORIENTATIONS*P4EST_FACES) + cf)/P4EST_FACES;

            const int8_t bo = BFAM_PXEST_ORIENTATION(f,nf,o);

            mapping[sk].np = rank;

            mapping[sk].ns = idnk;
            mapping[sk].nk = cks[h];
            mapping[sk].nf = nf;
            mapping[sk].nh = nh;

            mapping[sk].s  = idk;
            mapping[sk].k  = k;
            mapping[sk].f  = f;
            mapping[sk].h  = BFAM_PXEST_BOHTONH_INV(bo,h) + 1;
            mapping[sk].o  = bo;

            mapping[sk].i  = -1;
            mapping[sk].gi = -1;

            mapping[sk].id = glueid;

            ++sk;
          }
        }
      }
    }
  }
}

/** Count the number of boundary faces.
 *
 * \param [in]  mesh pxest mesh of the elements
 *
 * \return number of boundary faces.
 *
 */
static bfam_locidx_t
bfam_domain_pxest_num_boundary_faces(p4est_mesh_t *mesh)
{
  const p4est_locidx_t K = mesh->local_num_quadrants;

  bfam_locidx_t numBoundaryFaces = 0;

  for(p4est_locidx_t k = 0; k < K; ++k)
  {
    for (int f = 0; f < P4EST_FACES; ++f)
    {
      const p4est_locidx_t ck = mesh->quad_to_quad[P4EST_FACES * k + f];
      const int            cf = mesh->quad_to_face[P4EST_FACES * k + f];

      if(k == ck && f == cf)
        ++numBoundaryFaces;
    }
  }

  return numBoundaryFaces;
}

/** Build the boundary subdomain face mapping array.
 *
 * \param [in]  mesh             p4est mesh to build the mapping for.
 * \param [in]  subdomainID      subdomain id of each element in the mesh.
 * \param [in]  glueID           glue id of each face in the mesh
 * \param [in]  numBoundaryFaces the number of boundary subdomain faces in
 *                               the p4est mesh.
 * \param [out] mapping          the mapping array that will be filled.
 *
 */
static void
bfam_domain_pxest_boundary_subdomain_face_mapping(p4est_mesh_t *mesh,
    bfam_locidx_t *subdomainID, bfam_locidx_t *glueID,
    bfam_locidx_t numBoundaryFaces, bfam_subdomain_face_map_entry_t *mapping)
{
  const p4est_locidx_t K = mesh->local_num_quadrants;

  for(p4est_locidx_t k = 0, sk = 0; k < K; ++k)
  {
    const bfam_locidx_t  idk   = subdomainID[k];

    for (int f = 0; f < P4EST_FACES; ++f)
    {
      const p4est_locidx_t ck = mesh->quad_to_quad[P4EST_FACES * k + f];
      const int            cf = mesh->quad_to_face[P4EST_FACES * k + f];

      const bfam_locidx_t  glueid = (glueID) ? glueID[P4EST_FACES * k + f] : -1;

      if(k == ck && f == cf)
      {
        BFAM_ASSERT(sk < numBoundaryFaces);

        mapping[sk].np = -1;

        mapping[sk].ns = idk;
        mapping[sk].nk = k;
        mapping[sk].nf = f;
        mapping[sk].nh = 0;

        mapping[sk].s  = idk;
        mapping[sk].k  = k;
        mapping[sk].f  = f;
        mapping[sk].h  = 0;
        mapping[sk].o  = 0;

        mapping[sk].i  = -1;
        mapping[sk].gi = -1;

        mapping[sk].id = glueid;
        ++sk;
      }
    }
  }
}


void
bfam_domain_pxest_split_dgx_subdomains(bfam_domain_pxest_t *domain,
    bfam_locidx_t numSubdomains, bfam_locidx_t *subdomainID, int *N,
    bfam_locidx_t *glueID)
{
  BFAM_ROOT_LDEBUG("Begin splitting p4est domain into subdomains.");
  const int         HF = P4EST_HALF * P4EST_FACES;

  p4est_t       *pxest = domain->pxest;
  p4est_ghost_t *ghost = p4est_ghost_new(pxest, BFAM_PXEST_CONNECT);
  p4est_mesh_t  *mesh  = p4est_mesh_new(pxest, ghost, BFAM_PXEST_CONNECT);
  p4est_nodes_t *nodes = p4est_nodes_new(pxest, NULL);

  const p4est_locidx_t Nv = (p4est_locidx_t) nodes->indep_nodes.elem_count;

  bfam_long_real_t *VX = bfam_malloc_aligned(Nv * sizeof(bfam_long_real_t));
  bfam_long_real_t *VY = bfam_malloc_aligned(Nv * sizeof(bfam_long_real_t));
  bfam_long_real_t *VZ = bfam_malloc_aligned(Nv * sizeof(bfam_long_real_t));

  p4est_locidx_t *subK = bfam_calloc(numSubdomains,sizeof(p4est_locidx_t));
  p4est_locidx_t *subk = bfam_calloc(numSubdomains,sizeof(p4est_locidx_t));
  char          **name = bfam_malloc(numSubdomains*sizeof(char*));
  bfam_locidx_t **EToV = bfam_malloc(numSubdomains*sizeof(bfam_locidx_t*));
  bfam_locidx_t **EToE = bfam_malloc(numSubdomains*sizeof(bfam_locidx_t*));
  int8_t        **EToF = bfam_malloc(numSubdomains*sizeof(int8_t*));

  bfam_locidx_t *ktosubk =
    bfam_malloc(mesh->local_num_quadrants * sizeof(bfam_locidx_t));

  bfam_locidx_t *ghostSubdomainID =
    bfam_malloc(mesh->ghost_num_quadrants * sizeof(bfam_locidx_t));

  bfam_locidx_t *ghostN = bfam_malloc(mesh->ghost_num_quadrants * sizeof(int));

  bfam_locidx_t numParallelFaces = bfam_domain_pxest_num_parallel_faces(mesh);

  bfam_subdomain_face_map_entry_t *pfmapping
    = bfam_malloc_aligned(numParallelFaces*
        sizeof(bfam_subdomain_face_map_entry_t));

  bfam_domain_pxest_parallel_face_mapping(ghost, mesh, glueID, numParallelFaces,
      pfmapping);

  bfam_locidx_t numNeighbors =
    bfam_domain_pxest_parallel_face_num_neighbors(numParallelFaces, pfmapping);

  bfam_locidx_t *numNeighborFaces =
    bfam_malloc_aligned(numNeighbors*sizeof(bfam_locidx_t));

  bfam_locidx_t *neighborRank =
    bfam_malloc_aligned(numNeighbors*sizeof(bfam_locidx_t));

  bfam_domain_pxest_num_neighbor_faces(numParallelFaces, pfmapping,
      numNeighbors, numNeighborFaces, neighborRank);

#ifdef BFAM_DEBUG
  bfam_domain_pxest_parallel_face_mapping_check(pxest->mpicomm,
      numParallelFaces, pfmapping, numNeighbors, numNeighborFaces,
      neighborRank);
#endif

  bfam_domain_pxest_fill_ghost_subdomain_ids(pxest->mpicomm, numParallelFaces,
      pfmapping, numNeighbors, numNeighborFaces, neighborRank, subdomainID, N,
      ghostSubdomainID, ghostN);

  bfam_locidx_t numInterSubdomainFaces =
    bfam_domain_pxest_num_inter_subdomain_faces(mesh, subdomainID, glueID);

  BFAM_LDEBUG("numInterSubdomainFaces = %jd",
      (intmax_t) numInterSubdomainFaces);

  bfam_subdomain_face_map_entry_t *ifmapping
    = bfam_malloc_aligned(numInterSubdomainFaces*
        sizeof(bfam_subdomain_face_map_entry_t));

  int rank;
  BFAM_MPI_CHECK(MPI_Comm_rank(pxest->mpicomm, &rank));

  bfam_domain_pxest_inter_subdomain_face_mapping(rank, mesh, subdomainID,
      glueID, numInterSubdomainFaces, ifmapping);


  bfam_locidx_t numBoundaryFaces =
    bfam_domain_pxest_num_boundary_faces(mesh);

  BFAM_LDEBUG("numBoundaryFaces = %jd", (intmax_t) numBoundaryFaces);

  bfam_subdomain_face_map_entry_t *bfmapping
    = bfam_malloc_aligned(numBoundaryFaces*
        sizeof(bfam_subdomain_face_map_entry_t));

  bfam_domain_pxest_boundary_subdomain_face_mapping(mesh, subdomainID,
      glueID, numBoundaryFaces, bfmapping);

  for(p4est_locidx_t v = 0; v < Nv; ++v) {
    double xyz[3];
    p4est_quadrant_t *quad = p4est_quadrant_array_index (&nodes->indep_nodes, v);
    p4est_qcoord_to_vertex(pxest->connectivity,
                           quad->p.which_tree, quad->x, quad->y,
#if DIM==3
                           quad->z,
#endif
                           xyz);

    VX[v] = xyz[0];
    VY[v] = xyz[1];
    VZ[v] = xyz[2];
  }

  /*
   * Count the number of elements in each new subdomain
   */
  for(p4est_locidx_t k = 0; k < pxest->local_num_quadrants; ++k)
  {
    bfam_locidx_t id = (bfam_locidx_t) subdomainID[k];

    BFAM_ABORT_IF(id < 0 || id >= numSubdomains, "Bad Subdomain id: %jd",
        (intmax_t) id);

    ++subK[id];
  }

  for(bfam_locidx_t id = 0; id < numSubdomains; ++id)
  {
    name[id] = bfam_malloc(BFAM_BUFSIZ*sizeof(char));
    snprintf(name[id], BFAM_BUFSIZ, "dgx_dim_%1d_%05jd", (int) DIM,
        (intmax_t) id);

    EToV[id] = bfam_malloc(subK[id]*P4EST_CHILDREN*sizeof(bfam_locidx_t));
    EToE[id] = bfam_malloc(subK[id]*P4EST_FACES*sizeof(bfam_locidx_t));
    EToF[id] = bfam_malloc(subK[id]*P4EST_FACES*sizeof(int8_t));
  }

  const p4est_locidx_t K = mesh->local_num_quadrants;

  BFAM_ASSERT(K == pxest->local_num_quadrants);

  for(bfam_locidx_t id = 0; id < numSubdomains; ++id)
  {
    subk[id] = 0;
  }
  for(p4est_locidx_t k = 0; k < K; ++k)
  {
    const bfam_locidx_t  idk  = subdomainID[k];
    ktosubk[k] = subk[idk];
    ++subk[idk];
  }


  /*
   * Here we are decoding the p4est_mesh_t structure.  See p4est_mesh.h
   * for more details on how the data is stored.
   */

  /*
   * First build up the volume grids
   */
  for(bfam_locidx_t id = 0; id < numSubdomains; ++id)
  {
    subk[id] = 0;
  }
  for(p4est_locidx_t k = 0; k < K; ++k)
  {
    const bfam_locidx_t  idk  = subdomainID[k];

    for (int v = 0; v < P4EST_CHILDREN; ++v)
    {
      EToV[idk][P4EST_CHILDREN * subk[idk] + v] =
        nodes->local_nodes[P4EST_CHILDREN * k + v];
    }

    for (int f = 0; f < P4EST_FACES; ++f)
    {
      const p4est_locidx_t ck = mesh->quad_to_quad[P4EST_FACES * k + f];
      const int            cf = mesh->quad_to_face[P4EST_FACES * k + f];
      const bfam_locidx_t  glueid = (glueID) ? glueID[P4EST_FACES * k + f] : -1;

      int nk = k;
      int nf = f;

      if(glueid < 0 && cf >= 0 && cf < HF && ck < K && idk == subdomainID[ck])
      {
        /*
         * Neighbor is local to the subdomain and is the same size and does
         * not have a glue grid inbetween.
         */
        nk = ck;
        nf = cf;
      }

      EToE[idk][P4EST_FACES * subk[idk] + f] = ktosubk[nk];
      EToF[idk][P4EST_FACES * subk[idk] + f] = nf;
    }

    ++subk[idk];
  }


  bfam_subdomain_dgx_t **subdomains =
    bfam_malloc(numSubdomains*sizeof(bfam_subdomain_dgx_t**));
  for(bfam_locidx_t id = 0; id < numSubdomains; ++id)
  {
    const bfam_long_real_t *Vi[] = {VX, VY, VZ};
    subdomains[id] =
      bfam_subdomain_dgx_new(id,
                             -1,
                             name[id],
                             N[id],
                             Nv,
                             DIM,
                             Vi,
                             subK[id],
                             EToV[id],
                             EToE[id],
                             EToF[id],
                             DIM);

    bfam_subdomain_add_tag((bfam_subdomain_t *) subdomains[id], "_volume");
    bfam_domain_add_subdomain((bfam_domain_t    *) domain,
                              (bfam_subdomain_t *) subdomains[id]);
  }

  /*
   * Sort the local mapping
   */
  qsort(ifmapping, numInterSubdomainFaces,
      sizeof(bfam_subdomain_face_map_entry_t),
      bfam_subdomain_face_recv_cmp);

  /*
   * Setup the local glue grids
   */
  bfam_locidx_t numGlue = 0;

  for(bfam_locidx_t ifk = 0; ifk < numInterSubdomainFaces; )
  {
    /*
     * Count the number of element in the glue grid
     */
    bfam_locidx_t Kglue = 1;
    while(ifk + Kglue < numInterSubdomainFaces &&
        ifmapping[ifk+Kglue].np == ifmapping[ifk+Kglue-1].np &&
        ifmapping[ifk+Kglue].ns == ifmapping[ifk+Kglue-1].ns &&
        ifmapping[ifk+Kglue].s  == ifmapping[ifk+Kglue-1].s  &&
        ifmapping[ifk+Kglue].id == ifmapping[ifk+Kglue-1].id)
      ++Kglue;

    const bfam_locidx_t id_m = ifmapping[ifk].s;
    const bfam_locidx_t id_p = ifmapping[ifk].ns;
    const bfam_locidx_t glueid = ifmapping[ifk].id;

    int repeat = (id_m == id_p);
    repeat = 0;

    for(int r = 0; r <= repeat; ++r)
    {
      const bfam_locidx_t id = numSubdomains + numGlue;

      const bfam_locidx_t rank_m = rank;
      const bfam_locidx_t rank_p = rank;

      char glueName[BFAM_BUFSIZ];
      snprintf(glueName, BFAM_BUFSIZ,
          "dgx_dim_%1d_glue_%05jd_%05jd_%05jd_%05jd",
          (int) BDIM, (intmax_t) id, (intmax_t) id_m, (intmax_t) id_p,
          (intmax_t) glueid);

      /*
       * For subdomains that connect to themselves we need to distinguish
       * between them based on id.  So we have decided to use a minus sign
       * to distinguish between the two different glue grids.
       */
      // bfam_locidx_t sign_m = 1;
      bfam_locidx_t sign_p = 1;
      if(id_m == id_p)
      {
       // sign_m = (r) ? -1 :  1;
       sign_p = (r) ?  1 : -1;
      }

      bfam_subdomain_dgx_t *glue =
        bfam_subdomain_dgx_glue_new(id,
                                    glueid,
                                    glueName,
                                    N[id_m],
                                    N[id_p],
                                    rank_m,
                                    rank_p,
                                    sign_p * (id_m+1),
                                    sign_p * (id_p+1),
                                    subdomains[id_m],
                                    ktosubk,
                                    Kglue,
                                    ifmapping + ifk,
                                    DIM-1);

      bfam_subdomain_add_tag((bfam_subdomain_t *) glue, "_glue");
      bfam_subdomain_add_tag((bfam_subdomain_t *) glue, "_glue_local");
      char glue_num_tag[BFAM_BUFSIZ];
      snprintf(glue_num_tag,BFAM_BUFSIZ,"_glue_%d_%d",
          BFAM_MIN(id_m,id_p),BFAM_MAX(id_m,id_p));
      bfam_subdomain_add_tag((bfam_subdomain_t *) glue, glue_num_tag);

      char glue_id_tag[BFAM_BUFSIZ];
      snprintf(glue_id_tag,BFAM_BUFSIZ,"_glue_id_%d", glueid);
      bfam_subdomain_add_tag((bfam_subdomain_t *) glue, glue_id_tag);

      bfam_domain_add_subdomain((bfam_domain_t    *) domain,
                                (bfam_subdomain_t *) glue);

      numGlue += 1;
    }
    ifk += Kglue;
  }

  /*
   * Sort the boundary mapping
   */
  qsort(bfmapping, numBoundaryFaces,
      sizeof(bfam_subdomain_face_map_entry_t),
      bfam_subdomain_face_recv_cmp);

  /*
   * Setup the boundary glue grids
   */
  for(bfam_locidx_t bfk = 0; bfk < numBoundaryFaces; )
  {
    /*
     * Count the number of elements in the glue grid
     */
    bfam_locidx_t Kglue = 1;
    while(bfk + Kglue < numBoundaryFaces &&
        bfmapping[bfk+Kglue].np == bfmapping[bfk+Kglue-1].np &&
        bfmapping[bfk+Kglue].ns == bfmapping[bfk+Kglue-1].ns &&
        bfmapping[bfk+Kglue].s  == bfmapping[bfk+Kglue-1].s  &&
        bfmapping[bfk+Kglue].id == bfmapping[bfk+Kglue-1].id)
      ++Kglue;

    const bfam_locidx_t id = numSubdomains + numGlue;

    const bfam_locidx_t id_m = bfmapping[bfk].s;
    const bfam_locidx_t id_p = bfmapping[bfk].ns;

    const bfam_locidx_t rank_m = rank;
    const bfam_locidx_t rank_p = bfmapping[bfk].np;

    const bfam_locidx_t glueid = bfmapping[bfk].id;

    char glueName[BFAM_BUFSIZ];
    snprintf(glueName, BFAM_BUFSIZ,
        "dgx_dim_%1d_glue_b_%05jd_%05jd_%05jd_%05jd",
        (int) BDIM, (intmax_t) id, (intmax_t) id_m, (intmax_t) id_p,
        (intmax_t) glueid);


    bfam_subdomain_dgx_t *glue =
      bfam_subdomain_dgx_glue_new(id,
                                  glueid,
                                  glueName,
                                  N[id_m],
                                  N[id_p],
                                  rank_m,
                                  rank_p,
                                  id_m+1,
                                  id_p+1,
                                  subdomains[id_m],
                                  ktosubk,
                                  Kglue,
                                  bfmapping + bfk,
                                  DIM-1);

    bfam_subdomain_add_tag((bfam_subdomain_t *) glue, "_glue");
    bfam_subdomain_add_tag((bfam_subdomain_t *) glue, "_glue_boundary");

    char glue_id_tag[BFAM_BUFSIZ];
    snprintf(glue_id_tag,BFAM_BUFSIZ,"_glue_id_%d", glueid);
    bfam_subdomain_add_tag((bfam_subdomain_t *) glue, glue_id_tag);

    bfam_domain_add_subdomain((bfam_domain_t    *) domain,
                              (bfam_subdomain_t *) glue);

    numGlue += 1;
    bfk += Kglue;
  }

  /*
   * Sort the parallel mapping
   */
  qsort(pfmapping, numParallelFaces,
      sizeof(bfam_subdomain_face_map_entry_t),
      bfam_subdomain_face_recv_cmp);

  /*
   * Setup the parallel glue grids
   */
  for(bfam_locidx_t pfk = 0; pfk < numParallelFaces; )
  {
    /*
     * Count the number of elements in the glue grid
     */
    bfam_locidx_t Kglue = 1;
    while(pfk + Kglue < numParallelFaces &&
        pfmapping[pfk+Kglue].np == pfmapping[pfk+Kglue-1].np &&
        pfmapping[pfk+Kglue].ns == pfmapping[pfk+Kglue-1].ns &&
        pfmapping[pfk+Kglue].s  == pfmapping[pfk+Kglue-1].s  &&
        pfmapping[pfk+Kglue].id == pfmapping[pfk+Kglue-1].id)
      ++Kglue;

    const bfam_locidx_t id = numSubdomains + numGlue;

    const bfam_locidx_t id_m = pfmapping[pfk].s;
    const bfam_locidx_t id_p = pfmapping[pfk].ns;
    const bfam_locidx_t gid_p = pfmapping[pfk].gi;
    const bfam_locidx_t glueid = pfmapping[pfk].id;

    const bfam_locidx_t rank_m = rank;
    const bfam_locidx_t rank_p = pfmapping[pfk].np;

    char glueName[BFAM_BUFSIZ];
    snprintf(glueName, BFAM_BUFSIZ,
        "dgx_dim_%1d_glue_p_%05jd_%05jd_%05jd_%05jd",
        (int) BDIM, (intmax_t) id, (intmax_t) id_m, (intmax_t) id_p,
        (intmax_t) glueid);


    bfam_subdomain_dgx_t *glue =
      bfam_subdomain_dgx_glue_new(id,
                                  glueid,
                                  glueName,
                                  N[id_m],
                                  ghostN[gid_p],
                                  rank_m,
                                  rank_p,
                                  id_m+1,
                                  id_p+1,
                                  subdomains[id_m],
                                  ktosubk,
                                  Kglue,
                                  pfmapping + pfk,
                                  DIM-1);

    bfam_subdomain_add_tag((bfam_subdomain_t *) glue, "_glue");
    bfam_subdomain_add_tag((bfam_subdomain_t *) glue, "_glue_parallel");
      char glue_num_tag[BFAM_BUFSIZ];
      snprintf(glue_num_tag,BFAM_BUFSIZ,"_glue_%d_%d",
          BFAM_MIN(id_m,id_p),BFAM_MAX(id_m,id_p));
      bfam_subdomain_add_tag((bfam_subdomain_t *) glue, glue_num_tag);

    char glue_id_tag[BFAM_BUFSIZ];
    snprintf(glue_id_tag,BFAM_BUFSIZ,"_glue_id_%d", glueid);
    bfam_subdomain_add_tag((bfam_subdomain_t *) glue, glue_id_tag);

    bfam_domain_add_subdomain((bfam_domain_t    *) domain,
                              (bfam_subdomain_t *) glue);

    numGlue += 1;
    pfk += Kglue;
  }

  bfam_free(subdomains);

  for(bfam_locidx_t id = 0; id < numSubdomains; ++id)
  {
    bfam_free(name[id]);
    bfam_free(EToV[id]);
    bfam_free(EToE[id]);
    bfam_free(EToF[id]);
  }

  bfam_free(name);
  bfam_free(EToV);
  bfam_free(EToE);
  bfam_free(EToF);

  bfam_free(subK);
  bfam_free(subk);

  bfam_free(ktosubk);

  bfam_free_aligned(bfmapping);
  bfam_free_aligned(ifmapping);

  bfam_free_aligned(neighborRank);
  bfam_free_aligned(numNeighborFaces);
  bfam_free_aligned(pfmapping);
  bfam_free(ghostN);
  bfam_free(ghostSubdomainID);

  bfam_free_aligned(VX);
  bfam_free_aligned(VY);
  bfam_free_aligned(VZ);

  p4est_nodes_destroy(nodes);
  p4est_mesh_destroy(mesh);
  p4est_ghost_destroy(ghost);

  BFAM_ROOT_LDEBUG("End splitting pxest domain into subdomains.");
  bfam_domain_pxest_dgx_print_stats(domain);
}

#else

void
bfam_domain_pxest_stub() { }
#endif
