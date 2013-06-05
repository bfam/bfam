#include <bfam_subdomain.h>
#include <bfam_base.h>

void
bfam_subdomain_init(bfam_subdomain_t *thisSubdomain,const char* name)
{
  int len = strlen(name);
  thisSubdomain->name = bfam_malloc((len+1)*sizeof(char));
  strncpy(thisSubdomain->name,name,len+1);
  
  thisSubdomain->comm = NULL;
  thisSubdomain->hasWork = (0!=0);

  thisSubdomain->tags.root = NULL;

  thisSubdomain->start_communication = NULL;
  thisSubdomain->end_communication   = NULL;
  thisSubdomain->start_io            = NULL;
  thisSubdomain->start_end           = NULL;
  thisSubdomain->do_internal_RHS     = NULL;
  thisSubdomain->do_external_RHS     = NULL;
  thisSubdomain->update_fields       = NULL;
  thisSubdomain->free                = bfam_subdomain_free;
}

void
bfam_subdomain_free(bfam_subdomain_t *thisSubdomain)
{
  bfam_free(thisSubdomain->name);
  bfam_critbit0_clear(&thisSubdomain->tags);
  thisSubdomain->comm = NULL;
  thisSubdomain->hasWork = (0!=0);
  thisSubdomain->tags.root = NULL;
  thisSubdomain->start_communication = NULL;
  thisSubdomain->end_communication   = NULL;
  thisSubdomain->start_io            = NULL;
  thisSubdomain->start_end           = NULL;
  thisSubdomain->do_internal_RHS     = NULL;
  thisSubdomain->do_external_RHS     = NULL;
  thisSubdomain->update_fields       = NULL;
}

void
bfam_subdomain_add_tag(bfam_subdomain_t *thisSubdomain, const char* tag)
{
  int r = bfam_critbit0_insert(&thisSubdomain->tags, tag);

  BFAM_ABORT_IF(!r, "Out of memory when adding tag: %s to subdomain %s",
      tag, thisSubdomain->name);
}

int
bfam_subdomain_delete_tag(bfam_subdomain_t *thisSubdomain, const char* tag)
{
  return bfam_critbit0_delete(&thisSubdomain->tags, tag);
}

int
bfam_subdomain_has_tag(bfam_subdomain_t *thisSubdomain, const char* tag)
{
  return bfam_critbit0_contains(&thisSubdomain->tags, tag);
}
