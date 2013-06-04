#include <bfam_domain.h>
#include <bfam_base.h>
#include <bfam_log.h>

const char keyValueSplitter = '\033';

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
  thisDomain->subdomains =
    bfam_malloc(sizeSubdomains*sizeof(bfam_subdomain_t*));
  thisDomain->name2num.root = NULL;
}

void
bfam_domain_free(bfam_domain_t *thisDomain)
{
  for(int i = 0;i < thisDomain->numSubdomains;i++)
  {
    thisDomain->subdomains[i]->free(thisDomain->subdomains[i]);
    bfam_free(thisDomain->subdomains[i]);
  }
  thisDomain->comm = NULL;
  thisDomain->numSubdomains = 0;
  thisDomain->sizeSubdomains = 0;
  bfam_free(thisDomain->subdomains);
  thisDomain->subdomains = NULL;
  bfam_critbit0_clear(&thisDomain->name2num);
}

void
bfam_domain_add_subdomain(bfam_domain_t* thisDomain,
    bfam_subdomain_t* newSubdomain)
{
  // double size
  if(thisDomain->numSubdomains == thisDomain->sizeSubdomains)
  {
    BFAM_ROOT_VERBOSE("Doubling domain size");
    thisDomain->sizeSubdomains = 2*thisDomain->sizeSubdomains;
    thisDomain->subdomains =
      bfam_realloc(thisDomain->subdomains,
          thisDomain->sizeSubdomains*sizeof(bfam_subdomain_t*));
  }

  // create the key value pair
  BFAM_ROOT_VERBOSE("adding subdomain %3d with name %s",
      thisDomain->numSubdomains,newSubdomain->name);
  int len = strlen(newSubdomain->name);
  char* keyValue = bfam_malloc(sizeof(char)*(len+2)+sizeof(int));
  strncpy(keyValue,newSubdomain->name,len);
  keyValue[len] = keyValueSplitter;
  keyValue[len+1] = '\0';

  // check if it's already there
  if(bfam_critbit0_contains(&(thisDomain->name2num),keyValue))
  {
    BFAM_ABORT("domain already contains subdomain \"%s\"",
        newSubdomain->name);
  }

  // Since not in the table add the rest of the key and add
  *((int *) (keyValue+len+1)) = thisDomain->numSubdomains;
  keyValue[len+1+sizeof(int)] = '\0';
  bfam_critbit0_insert(&(thisDomain->name2num),keyValue);

  // add block
  thisDomain->subdomains[thisDomain->numSubdomains] = newSubdomain;
  thisDomain->numSubdomains++;

  bfam_free(keyValue);
}
