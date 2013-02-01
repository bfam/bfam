#include <bfam_subdomain.h>
#include <bfam_base.h>

void bfam_domain_new(bfam_domain_t *thisDomain, bfam_mpicomm_t *domComm)
{
  const int sizeSubdomains = 16;
  thisDomain->comm = domComm; // Perhaps we should duplicate it?
  thisDomain->numSubdomains = 0;
  thisDomain->sizeSubdomains = sizeSubdomains;
  thisDomain->subdomains = bfam_malloc(sizeSubdomains*sizeof(bfam_subdomain_t*));
}

void bfam_domain_free(bfam_domain_t *thisDomain)
{
  thisDomain->comm = NULL;
  thisDomain->numSubdomains = 0;
  thisDomain->sizeSubdomains = 0;
  bfam_free(thisDomain->subdomains);
  thisDomain->subdomains = NULL;
}
