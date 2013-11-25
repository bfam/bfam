#include <bfam_base.h>
#include <bfam_domain_p8est.h>
#include <bfam_log.h>
#include <bfam_subdomain_dgx.h>

bfam_domain_p8est_t*
bfam_domain_p8est_new(MPI_Comm domComm,
                      p8est_connectivity_t *conn)
{
  bfam_domain_p8est_t* newDomain = bfam_malloc(sizeof(bfam_domain_p8est_t));
  bfam_domain_p8est_init(newDomain,domComm,conn);
  return newDomain;
}

void
bfam_domain_p8est_init(bfam_domain_p8est_t *domain, MPI_Comm domComm,
                       p8est_connectivity_t *conn)
{
  bfam_domain_init(&domain->base,domComm);

  domain->conn = conn;
  domain->p8est = p8est_new(domComm,conn,0,NULL,NULL);
}

void
bfam_domain_p8est_free(bfam_domain_p8est_t *domain)
{
  /* Memory we don't manage */
  domain->conn = NULL;

  /* Memory we do manage */
  p8est_destroy(domain->p8est);
  domain->p8est = NULL;

  bfam_domain_free(&domain->base);
}

void
bfam_domain_p8est_split_dgx_subdomains(bfam_domain_p8est_t *domain,
    bfam_locidx_t numSubdomains, bfam_locidx_t *subdomainID, int *N)
{
  BFAM_ABORT("Not implemented");
}
