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
