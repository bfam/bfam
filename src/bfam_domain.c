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
  const bfam_locidx_t sizeSubdomains = 16;
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
  for(bfam_locidx_t i = 0;i < thisDomain->numSubdomains;i++)
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
  size_t len = strlen(newSubdomain->name);
  char* keyValue = bfam_malloc(sizeof(char)*(len+2)+sizeof(bfam_locidx_t));
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
  *((bfam_locidx_t *) (keyValue+len+1)) = thisDomain->numSubdomains;
  keyValue[len+1+sizeof(bfam_locidx_t)] = '\0';
  bfam_critbit0_insert(&(thisDomain->name2num),keyValue);

  // add block
  thisDomain->subdomains[thisDomain->numSubdomains] = newSubdomain;
  thisDomain->numSubdomains++;

  bfam_free(keyValue);
}

void
bfam_domain_get_subdomains(bfam_domain_t *thisDomain,
    bfam_domain_match_t matchType, const char **tags,
    bfam_locidx_t numEntries, bfam_subdomain_t **subdomains,
    bfam_locidx_t *numSubdomains)
{
  BFAM_ASSERT(subdomains != NULL);
  BFAM_ASSERT(numSubdomains != NULL);

  if(numEntries<=0)
    return;

  *numSubdomains = 0;
  for(bfam_locidx_t d = 0; d < thisDomain->numSubdomains; ++d)
  {
    bfam_subdomain_t *subdomain = thisDomain->subdomains[d];
    int matched = 0;
    switch(matchType)
    {
      case BFAM_DOMAIN_OR:
        matched = 0;
        for(size_t t = 0; !matched && tags[t]; ++t)
        {
          int hasTag = bfam_subdomain_has_tag(subdomain, tags[t]);
          matched = hasTag || matched;
        }
        break;
      case BFAM_DOMAIN_AND:
        matched = 1;
        for(size_t t = 0; matched && tags[t]; ++t)
          matched = matched && bfam_subdomain_has_tag(subdomain, tags[t]);
        break;
      default:
        BFAM_ABORT("Unsupported Match Type");
    }

    if(matched)
    {
      subdomains[*numSubdomains] = subdomain;
      ++(*numSubdomains);
    }

    if(*numSubdomains == numEntries)
      return;
  }

  return;
}

int
bfam_domain_get_subdomains_critbit_or(const char *elem,void *arg)
{
  int* matched = (int *)((void**)arg)[0];
  bfam_subdomain_t *subdomain = (bfam_domain_t*)((void**)arg)[1];

  (*matched) = (*matched) || bfam_subdomain_has_tag(subdomain, elem);
  return !(*matched);
}


int
bfam_domain_get_subdomains_critbit_and(const char *elem,void *arg)
{
  int* matched = (int *)((void**)arg)[0];
  bfam_subdomain_t *subdomain = (bfam_domain_t*)((void**)arg)[1];

  (*matched) = (*matched) && bfam_subdomain_has_tag(subdomain, elem);
  return 1;
}


void
bfam_domain_get_subdomains_critbit(bfam_domain_t *thisDomain,
    bfam_domain_match_t matchType, bfam_critbit0_tree_t *tags,
    bfam_locidx_t numEntries, bfam_subdomain_t **subdomains,
    bfam_locidx_t *numSubdomains)
{
  BFAM_ASSERT(subdomains != NULL);
  BFAM_ASSERT(numSubdomains != NULL);


  if(numEntries<=0)
    return;

  *numSubdomains = 0;
  for(bfam_locidx_t d = 0; d < thisDomain->numSubdomains; ++d)
  {
    bfam_subdomain_t *subdomain = thisDomain->subdomains[d];
    int matched = 0;
    void *arg[2];
    arg[0] = &matched;
    arg[1] = subdomain;
    switch(matchType)
    {
      case BFAM_DOMAIN_OR:
        matched = 0;
        bfam_critbit0_allprefixed(tags, "",
            &bfam_domain_get_subdomains_critbit_or, arg);
        break;
      case BFAM_DOMAIN_AND:
        matched = 1;
        bfam_critbit0_allprefixed(tags, "",
            &bfam_domain_get_subdomains_critbit_and, arg);
        break;
      default:
        BFAM_ABORT("Unsupported Match Type");
    }

    if(matched)
    {
      subdomains[*numSubdomains] = subdomain;
      ++(*numSubdomains);
    }

    if(*numSubdomains == numEntries)
      return;
  }

  return;
}
