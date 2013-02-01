#include <bfam_domain.h>
#include <bfam_base.h>

bfam_domain_t* bfam_domain_new(bfam_mpicomm_t *domComm)
{
  bfam_domain_t* newDomain = bfam_malloc(sizeof(bfam_domain_t));
  bfam_domain_init(newDomain,domComm);
  return newDomain;
}

void
bfam_domain_init(bfam_domain_t *thisDomain, bfam_mpicomm_t *domComm)
{
  const int sizeSubdomains = 16;
  thisDomain->comm = domComm; // Perhaps we should duplicate it?
  thisDomain->numSubdomains = 0;
  thisDomain->sizeSubdomains = sizeSubdomains;
  thisDomain->subdomains = bfam_malloc(sizeSubdomains*sizeof(bfam_subdomain_t*));
}

void
bfam_domain_free(bfam_domain_t *thisDomain)
{
  thisDomain->comm = NULL;
  thisDomain->numSubdomains = 0;
  thisDomain->sizeSubdomains = 0;
  bfam_free(thisDomain->subdomains);
  thisDomain->subdomains = NULL;
}

void
bfam_domain_add_subdomain(bfam_domain_t* thisDomain, bfam_subdomain_t* newSubdomain)
{
  // double size
  if(thisDomain->numSubdomains == thisDomain->sizeSubdomains)
  {
    bfam_subdomain_t** oldSubdomains = thisDomain->subdomains;
    thisDomain->sizeSubdomains = 2*thisDomain->sizeSubdomains;
    thisDomain->subdomains =
      bfam_malloc(thisDomain->sizeSubdomains*sizeof(bfam_subdomain_t*));
    for(int i = 0;i < thisDomain->numSubdomains;i++)
      thisDomain->subdomains[i] = oldSubdomains[i];
    bfam_free(oldSubdomains);
  }

  // Add next domain
  thisDomain->subdomains[thisDomain->numSubdomains] = newSubdomain;
  thisDomain->numSubdomains++;
}
