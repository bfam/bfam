#include <bfam_base.h>
#include <bfam_domain_p4est.h>
#include <bfam_log.h>
#include <bfam_subdomain_dgx_quad.h>

bfam_domain_p4est_t*
bfam_domain_p4est_new(MPI_Comm domComm,
                      p4est_connectivity_t *conn)
{
  bfam_domain_p4est_t* newDomain = bfam_malloc(sizeof(bfam_domain_p4est_t));
  bfam_domain_p4est_init(newDomain,domComm,conn);
  return newDomain;
}

void
bfam_domain_p4est_init(bfam_domain_p4est_t *domain, MPI_Comm domComm,
                       p4est_connectivity_t *conn)
{
  bfam_domain_init(&domain->base,domComm);

  domain->conn = conn;
  domain->p4est = p4est_new_ext(domComm,conn,0,0,0,0,NULL,NULL);
}

void
bfam_domain_p4est_free(bfam_domain_p4est_t *domain)
{
  /* Memory we don't manage */
  domain->conn = NULL;

  /* Memory we do manage */
  p4est_destroy(domain->p4est);
  domain->p4est = NULL;

  bfam_domain_free(&domain->base);
}

static bfam_locidx_t
bfam_domain_p4est_num_parallel_faces(p4est_mesh_t *mesh)
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
 * and \c bfam_domain_p4est_face_recv_cmp comparison
 * functions.
 *
 * \param [in]  mesh             p4est mesh to build the mapping for.
 * \param [in]  numParallelFaces the number of parallel faces in the
 *                               p4est mesh.
 * \param [out] mapping          the mapping array that will be filled.
 *
 */
static void
bfam_domain_p4est_parallel_face_mapping(p4est_mesh_t *mesh,
    bfam_locidx_t numParallelFaces, bfam_subdomain_face_map_entry_t *mapping)
{
  const p4est_locidx_t K = mesh->local_num_quadrants;

  for(p4est_locidx_t k = 0, sk = 0; k < K; ++k)
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
        {
          p4est_locidx_t ghostid = ck-mesh->local_num_quadrants;
          p4est_locidx_t ghostp  = mesh->ghost_to_proc[ghostid];
          p4est_locidx_t ghostk  = mesh->ghost_to_index[ghostid];
          int8_t         ghostf = cf;
          int8_t         ghosth = 0;

          if(ghostf >= 8)
          {
            ghostf -= 8;

            if(ghostf >= 8)
            {
              ghostf -= 8;
              ghosth  = 2;
            }
            else
            {
              ghosth  = 1;
            }
          }

          int8_t o = ghostf/4;

          ghostf = ghostf%4;

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
          mapping[sk].o  = o;
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
            p4est_locidx_t ghostp  = mesh->ghost_to_proc[ghostid];
            p4est_locidx_t ghostk  = mesh->ghost_to_index[ghostid];
            int8_t         ghostf  = (8 + cf)%4;
            int8_t         ghosth  = 0;
            int8_t         o       = (8 + cf)/4;

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
            mapping[sk].h  = h + 1;
            mapping[sk].o  = o;
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
    BFAM_LDEBUG("mapping: %5s %5s %5s %5s %5s %5s %5s %5s %5s %5s",
        "np", "nk", "nf", "nh", "gi", "i", "k", "f", "h", "o");
    for(bfam_locidx_t i = 0; i < numParallelFaces; ++i)
    {
      BFAM_LDEBUG("mapping: %5jd %5jd %5jd %5jd %5jd %5jd %5jd %5jd %5jd %5jd",
          (intmax_t)mapping[i].np,
          (intmax_t)mapping[i].nk,
          (intmax_t)mapping[i].nf,
          (intmax_t)mapping[i].nh,
          (intmax_t)mapping[i].gi,
          (intmax_t)mapping[i].i,
          (intmax_t)mapping[i].k,
          (intmax_t)mapping[i].f,
          (intmax_t)mapping[i].h,
          (intmax_t)mapping[i].o);
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
bfam_domain_p4est_parallel_face_num_neighbors(bfam_locidx_t numParallelFaces,
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
bfam_domain_p4est_num_neighbor_faces(bfam_locidx_t numParallelFaces,
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
bfam_domain_p4est_parallel_face_send_recv(MPI_Comm comm,
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
bfam_domain_p4est_parallel_face_mapping_check(MPI_Comm comm,
    bfam_locidx_t numParallelFaces, bfam_subdomain_face_map_entry_t *mapping,
    bfam_locidx_t numNeighbors, bfam_locidx_t *numNeighborFaces,
    bfam_locidx_t *neighborRank)
{
  const int entries = 8;

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
  }

  bfam_domain_p4est_parallel_face_send_recv(comm, numParallelFaces, mapping,
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
bfam_domain_p4est_fill_ghost_subdomain_ids(MPI_Comm comm,
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

  bfam_domain_p4est_parallel_face_send_recv(comm, numParallelFaces, mapping,
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
 *
 * \return number of total local inter-subdomain faces.
 *
 */
static bfam_locidx_t
bfam_domain_p4est_num_inter_subdomain_faces(p4est_mesh_t *mesh,
    bfam_locidx_t *subdomainID)
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

          if(cf >= 8)
            hanging = 1;
          else
            hanging = 0;

          /*
           * Count intra and inter subdomain faces.
           *
           * Only count same subdomain to same subdomain if it is a hanging
           * face.
           */
          if((idnk == idk && hanging) || idnk != idk)
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
 * \param [in]  mesh                   p4est mesh to build the mapping for.
 * \param [in]  subdomainID            subdomain id of each element in the mesh.
 * \param [in]  numInterSubdomainFaces the number of inter subdomain faces in
 *                                     the p4est mesh.
 * \param [out] mapping                the mapping array that will be filled.
 *
 */
static void
bfam_domain_p4est_inter_subdomain_face_mapping(bfam_locidx_t rank,
    p4est_mesh_t *mesh, bfam_locidx_t *subdomainID,
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

          if(cf >= 8)
            hanging = 1;
          else
            hanging = 0;

          int8_t         nf = cf;
          int8_t         nh = 0;

          if(nf >= 8)
          {
            nf -= 8;

            if(nf >= 8)
            {
              nf -= 8;
              nh  = 2;
            }
            else
            {
              nh = 1;
            }
          }

          int8_t o = nf/4;

          nf = nf%4;

          /*
           * Count intra and inter subdomain faces.
           *
           * Only count same subdomain to same subdomain if it is a hanging
           * face.
           */
          if((idnk == idk && hanging) || idnk != idk)
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
            mapping[sk].o  = o;

            mapping[sk].i  = -1;
            mapping[sk].gi = -1;
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

            int8_t         nf = (8 + cf)%4;
            int8_t         nh = 0;
            int8_t         o  = (8 + cf)/4;

            mapping[sk].np = rank;

            mapping[sk].ns = idnk;
            mapping[sk].nk = cks[h];
            mapping[sk].nf = nf;
            mapping[sk].nh = nh;

            mapping[sk].s  = idk;
            mapping[sk].k  = k;
            mapping[sk].f  = f;
            mapping[sk].h  = h + 1;
            mapping[sk].o  = o;

            mapping[sk].i  = -1;
            mapping[sk].gi = -1;

            ++sk;
          }
        }
      }
    }
  }
}

/** Count the number of boundary faces.
 *
 * \param [in]  mesh p4est mesh of the elements
 *
 * \return number of boundary faces.
 *
 */
static bfam_locidx_t
bfam_domain_p4est_num_boundary_faces(p4est_mesh_t *mesh)
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
 * \param [in]  numBoundaryFaces the number of boundary subdomain faces in
 *                               the p4est mesh.
 * \param [out] mapping          the mapping array that will be filled.
 *
 */
static void
bfam_domain_p4est_boundary_subdomain_face_mapping(p4est_mesh_t *mesh,
    bfam_locidx_t *subdomainID, bfam_locidx_t numBoundaryFaces,
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
        ++sk;
      }
    }
  }
}

void
bfam_domain_p4est_split_dgx_quad_subdomains(bfam_domain_p4est_t *domain,
    bfam_locidx_t numSubdomains, bfam_locidx_t *subdomainID, int *N)
{
  BFAM_ROOT_LDEBUG("Begin splitting p4est domain into subdomains.");
  const int         HF = P4EST_HALF * P4EST_FACES;

  p4est_t       *p4est = domain->p4est;
  p4est_ghost_t *ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FACE);
  p4est_mesh_t  *mesh  = p4est_mesh_new(p4est, ghost, P4EST_CONNECT_FACE);

  const p4est_locidx_t Nv = mesh->local_num_vertices;

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

  bfam_locidx_t numParallelFaces = bfam_domain_p4est_num_parallel_faces(mesh);

  bfam_subdomain_face_map_entry_t *pfmapping
    = bfam_malloc_aligned(numParallelFaces*
        sizeof(bfam_subdomain_face_map_entry_t));

  bfam_domain_p4est_parallel_face_mapping(mesh, numParallelFaces, pfmapping);

  bfam_locidx_t numNeighbors =
    bfam_domain_p4est_parallel_face_num_neighbors(numParallelFaces, pfmapping);

  bfam_locidx_t *numNeighborFaces =
    bfam_malloc_aligned(numNeighbors*sizeof(bfam_locidx_t));

  bfam_locidx_t *neighborRank =
    bfam_malloc_aligned(numNeighbors*sizeof(bfam_locidx_t));

  bfam_domain_p4est_num_neighbor_faces(numParallelFaces, pfmapping,
      numNeighbors, numNeighborFaces, neighborRank);

#ifdef BFAM_DEBUG
  bfam_domain_p4est_parallel_face_mapping_check(p4est->mpicomm,
      numParallelFaces, pfmapping, numNeighbors, numNeighborFaces,
      neighborRank);
#endif

  bfam_domain_p4est_fill_ghost_subdomain_ids(p4est->mpicomm, numParallelFaces,
      pfmapping, numNeighbors, numNeighborFaces, neighborRank, subdomainID, N,
      ghostSubdomainID, ghostN);

  bfam_locidx_t numInterSubdomainFaces =
    bfam_domain_p4est_num_inter_subdomain_faces(mesh, subdomainID);

  BFAM_LDEBUG("numInterSubdomainFaces = %jd", (intmax_t) numInterSubdomainFaces);

  bfam_subdomain_face_map_entry_t *ifmapping
    = bfam_malloc_aligned(numInterSubdomainFaces*
        sizeof(bfam_subdomain_face_map_entry_t));

  int rank;
  BFAM_MPI_CHECK(MPI_Comm_rank(p4est->mpicomm, &rank));

  bfam_domain_p4est_inter_subdomain_face_mapping(rank, mesh, subdomainID,
      numInterSubdomainFaces, ifmapping);


  bfam_locidx_t numBoundaryFaces =
    bfam_domain_p4est_num_boundary_faces(mesh);

  BFAM_LDEBUG("numBoundaryFaces = %jd", (intmax_t) numBoundaryFaces);

  bfam_subdomain_face_map_entry_t *bfmapping
    = bfam_malloc_aligned(numBoundaryFaces*
        sizeof(bfam_subdomain_face_map_entry_t));

  bfam_domain_p4est_boundary_subdomain_face_mapping(mesh, subdomainID,
      numBoundaryFaces, bfmapping);

  /*
   * Get vertex coordinates
   */
  for(p4est_locidx_t v = 0; v < mesh->local_num_vertices; ++v)
  {
    VX[v] = (bfam_long_real_t) mesh->vertices[3*v + 0];
    VY[v] = (bfam_long_real_t) mesh->vertices[3*v + 1];
    VZ[v] = (bfam_long_real_t) mesh->vertices[3*v + 2];
  }

  /*
   * Count the number of elements in each new subdomain
   */
  for(p4est_locidx_t k = 0; k < p4est->local_num_quadrants; ++k)
  {
    bfam_locidx_t id = (bfam_locidx_t) subdomainID[k];

    BFAM_ABORT_IF(id < 0 || id >= numSubdomains, "Bad Subdomain id: %jd",
        (intmax_t) id);

    ++subK[id];
  }

  for(bfam_locidx_t id = 0; id < numSubdomains; ++id)
  {
    name[id] = bfam_malloc(BFAM_BUFSIZ*sizeof(char));
    snprintf(name[id], BFAM_BUFSIZ, "dg_quad_%05jd", (intmax_t) id);

    EToV[id] = bfam_malloc(subK[id]*P4EST_CHILDREN*sizeof(bfam_locidx_t));
    EToE[id] = bfam_malloc(subK[id]*P4EST_FACES*sizeof(bfam_locidx_t));
    EToF[id] = bfam_malloc(subK[id]*P4EST_FACES*sizeof(int8_t));
  }

  const p4est_locidx_t K = mesh->local_num_quadrants;

  BFAM_ASSERT(K == p4est->local_num_quadrants);

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
        mesh->quad_to_vertex[P4EST_CHILDREN * k + v];
    }

    for (int f = 0; f < P4EST_FACES; ++f)
    {
      const p4est_locidx_t ck = mesh->quad_to_quad[P4EST_FACES * k + f];
      const int            cf = mesh->quad_to_face[P4EST_FACES * k + f];

      int nk = k;
      int nf = f;

      if(cf >= 0 && cf < HF && ck < K && idk == subdomainID[ck])
      {
        /*
         * Neighbor is local to the subdomain and is the same size.
         */
        nk = ck;
        nf = cf;
      }

      EToE[idk][P4EST_FACES * subk[idk] + f] = ktosubk[nk];
      EToF[idk][P4EST_FACES * subk[idk] + f] = nf;
    }

    ++subk[idk];
  }


  bfam_subdomain_dgx_quad_t **subdomains =
    bfam_malloc(numSubdomains*sizeof(bfam_subdomain_t**));
  for(bfam_locidx_t id = 0; id < numSubdomains; ++id)
  {
    subdomains[id] =
      bfam_subdomain_dgx_quad_new(id,
                                  name[id],
                                  N[id],
                                  Nv,
                                  VX,
                                  VY,
                                  VZ,
                                  subK[id],
                                  EToV[id],
                                  EToE[id],
                                  EToF[id]);

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
        ifmapping[ifk+Kglue].s  == ifmapping[ifk+Kglue-1].s)
      ++Kglue;

    const bfam_locidx_t id_m = ifmapping[ifk].s;
    const bfam_locidx_t id_p = ifmapping[ifk].ns;

    int repeat = (id_m == id_p);
    repeat = 0;

    for(int r = 0; r <= repeat; ++r)
    {
      const bfam_locidx_t id = numSubdomains + numGlue;

      const bfam_locidx_t rank_m = rank;
      const bfam_locidx_t rank_p = rank;

      char glueName[BFAM_BUFSIZ];
      snprintf(glueName, BFAM_BUFSIZ, "dg_quad_glue_%d_%05jd_%05jd_%05jd",
          r, (intmax_t) id, (intmax_t) id_m, (intmax_t) id_p);

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

      bfam_subdomain_dgx_quad_glue_t *glue =
        bfam_subdomain_dgx_quad_glue_new(id,
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
                                         ifmapping + ifk);

      bfam_subdomain_add_tag((bfam_subdomain_t *) glue, "_glue");
      bfam_subdomain_add_tag((bfam_subdomain_t *) glue, "_glue_local");
      char glue_num_tag[BFAM_BUFSIZ];
      snprintf(glue_num_tag,BFAM_BUFSIZ,"_glue_%d_%d",
          BFAM_MIN(id_m,id_p),BFAM_MAX(id_m,id_p));
      bfam_subdomain_add_tag((bfam_subdomain_t *) glue, glue_num_tag);

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
        bfmapping[bfk+Kglue].s  == bfmapping[bfk+Kglue-1].s)
      ++Kglue;

    const bfam_locidx_t id = numSubdomains + numGlue;

    const bfam_locidx_t id_m = bfmapping[bfk].s;
    const bfam_locidx_t id_p = bfmapping[bfk].ns;

    const bfam_locidx_t rank_m = rank;
    const bfam_locidx_t rank_p = bfmapping[bfk].np;

    char glueName[BFAM_BUFSIZ];
    snprintf(glueName, BFAM_BUFSIZ, "dg_quad_glue_b_%05jd_%05jd_%05jd",
        (intmax_t) id, (intmax_t) id_m, (intmax_t) id_p);


    bfam_subdomain_dgx_quad_glue_t *glue =
      bfam_subdomain_dgx_quad_glue_new(id,
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
                                       bfmapping + bfk);

    bfam_subdomain_add_tag((bfam_subdomain_t *) glue, "_glue");
    bfam_subdomain_add_tag((bfam_subdomain_t *) glue, "_glue_boundary");

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
        pfmapping[pfk+Kglue].s  == pfmapping[pfk+Kglue-1].s)
      ++Kglue;

    const bfam_locidx_t id = numSubdomains + numGlue;

    const bfam_locidx_t id_m = pfmapping[pfk].s;
    const bfam_locidx_t id_p = pfmapping[pfk].ns;
    const bfam_locidx_t gid_p = pfmapping[pfk].gi;

    const bfam_locidx_t rank_m = rank;
    const bfam_locidx_t rank_p = pfmapping[pfk].np;

    char glueName[BFAM_BUFSIZ];
    snprintf(glueName, BFAM_BUFSIZ, "dg_quad_glue_p_%05jd_%05jd_%05jd",
        (intmax_t) id, (intmax_t) id_m, (intmax_t) id_p);


    bfam_subdomain_dgx_quad_glue_t *glue =
      bfam_subdomain_dgx_quad_glue_new(id,
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
                                       pfmapping + pfk);

    bfam_subdomain_add_tag((bfam_subdomain_t *) glue, "_glue");
    bfam_subdomain_add_tag((bfam_subdomain_t *) glue, "_glue_parallel");
      char glue_num_tag[BFAM_BUFSIZ];
      snprintf(glue_num_tag,BFAM_BUFSIZ,"_glue_%d_%d",
          BFAM_MIN(id_m,id_p),BFAM_MAX(id_m,id_p));
      bfam_subdomain_add_tag((bfam_subdomain_t *) glue, glue_num_tag);

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

  p4est_mesh_destroy(mesh);
  p4est_ghost_destroy(ghost);

  BFAM_ROOT_LDEBUG("End splitting p4est domain into subdomains.");
}
