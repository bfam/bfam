#include <bfam_domain_p4est.h>
#include <bfam_log.h>

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
  p4est_t *p4est = domain->p4est;

  p4est_locidx_t *subK = bfam_calloc(numSubdomains,sizeof(p4est_locidx_t));
  char          **name = bfam_malloc(numSubdomains*sizeof(char*));
  bfam_locidx_t  *Nv   = bfam_calloc(numSubdomains,sizeof(bfam_locidx_t));
  bfam_locidx_t **EToV = bfam_malloc(numSubdomains*sizeof(bfam_locidx_t*));
  bfam_locidx_t **EToE = bfam_malloc(numSubdomains*sizeof(bfam_locidx_t*));
  int8_t        **EToF = bfam_malloc(numSubdomains*sizeof(int8_t*));

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
    snprintf(name[id], BFAM_BUFSIZ, "DG QUAD %05jd", (intmax_t) id);

    EToV[id] = bfam_malloc(subK[id]*4*sizeof(bfam_locidx_t));
    EToE[id] = bfam_malloc(subK[id]*4*sizeof(bfam_locidx_t));
    EToF[id] = bfam_malloc(subK[id]*4*sizeof(int8_t));
  }

  //   /*
  //    * Loop through all of the p4est trees
  //    */
  //   for (p4est_topidx_t t  = p4est->first_local_tree;
  //                       t <= p4est->last_local_tree;
  //                       ++t)
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

  bfam_free(Nv);
  bfam_free(subK);
}
