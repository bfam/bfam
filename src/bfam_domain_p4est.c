#include <bfam_base.h>
#include <bfam_domain_p4est.h>
#include <bfam_log.h>
#include <bfam_subdomain_dgx_quad.h>

bfam_domain_p4est_t*
bfam_domain_p4est_new(bfam_mpicomm_t *domComm,
                      p4est_connectivity_t *conn)
{
  bfam_domain_p4est_t* newDomain = bfam_malloc(sizeof(bfam_domain_p4est_t));
  bfam_domain_p4est_init(newDomain,domComm,conn);
  return newDomain;
}

void
bfam_domain_p4est_init(bfam_domain_p4est_t *domain, bfam_mpicomm_t *domComm,
                       p4est_connectivity_t *conn)
{
  bfam_domain_init(&domain->d,domComm);

  domain->conn = conn;
  domain->p4est = p4est_new_ext(domComm->comm,conn,0,0,0,0,NULL,NULL);
}

void
bfam_domain_p4est_free(bfam_domain_p4est_t *domain)
{
  /* Memory we don't manage */
  domain->conn = NULL;

  /* Memory we do manage */
  p4est_destroy(domain->p4est);
  domain->p4est = NULL;

  bfam_domain_free(&domain->d);
}

void
bfam_domain_p4est_split_dgx_quad_subdomains(bfam_domain_p4est_t *domain,
    bfam_locidx_t numSubdomains, bfam_locidx_t *subdomainID, int *N)
{
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
    bfam_malloc(mesh->local_num_quadrants * sizeof(p4est_locidx_t));

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
    p4est_locidx_t id = (p4est_locidx_t) subdomainID[k];

    BFAM_ABORT_IF(id < 0 || id >= numSubdomains, "Bad Subdomain id: %jd",
        (intmax_t) id);

    ++subK[id];
  }

  for(p4est_locidx_t id = 0; id < numSubdomains; ++id)
  {
    name[id] = bfam_malloc(BFAM_BUFSIZ*sizeof(char));
    snprintf(name[id], BFAM_BUFSIZ, "dg_quad_%05jd", (intmax_t) id);

    EToV[id] = bfam_malloc(subK[id]*P4EST_CHILDREN*sizeof(bfam_locidx_t));
    EToE[id] = bfam_malloc(subK[id]*P4EST_FACES*sizeof(bfam_locidx_t));
    EToF[id] = bfam_malloc(subK[id]*P4EST_FACES*sizeof(int8_t));
  }

  const p4est_locidx_t K = mesh->local_num_quadrants;

  BFAM_ASSERT(K == p4est->local_num_quadrants);

  for(p4est_locidx_t id = 0; id < numSubdomains; ++id)
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
  for(p4est_locidx_t id = 0; id < numSubdomains; ++id)
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


  for(p4est_locidx_t id = 0; id < numSubdomains; ++id)
  {
    bfam_subdomain_dgx_quad_t *subdomain =
      bfam_subdomain_dgx_quad_new(name[id],
                                  N[id],
                                  Nv,
                                  VX,
                                  VY,
                                  VZ,
                                  subK[id],
                                  EToV[id],
                                  EToE[id],
                                  EToF[id]);

    bfam_subdomain_add_tag((bfam_subdomain_t *) subdomain, "_volume");
    bfam_domain_add_subdomain((bfam_domain_t    *) domain,
                              (bfam_subdomain_t *) subdomain);
  }

  //   {
  //   }

  //   /*
  //   bfam_subdomain_dgx_quad_t subdomain =
  //     bfam_subdomain_dgx_quad_new(name,
  //                                 N[id],
  //                                 Nv,
  //                                 VX,
  //                                 VY,
  //                                 K,
  //                                 EToV,
  //                                 EToE,
  //                                 EToF);
  //                                 */

  for(p4est_locidx_t id = 0; id < numSubdomains; ++id)
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

  bfam_free_aligned(VX);
  bfam_free_aligned(VY);
  bfam_free_aligned(VZ);

  p4est_mesh_destroy(mesh);
  p4est_ghost_destroy(ghost);
}
