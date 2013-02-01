#include <bfam_subdomain.h>
#include <bfam_base.h>

void
bfam_subdomain_init(bfam_subdomain_t *thisSubdomain,const char* name)
{
  int len = strlen(name);
  thisSubdomain->name = bfam_malloc(len*sizeof(char));
  strncpy(thisSubdomain->name,name,len);
  
  thisSubdomain->comm = NULL;
  thisSubdomain->hasWork = (0!=0);

  thisSubdomain->start_communication = NULL;
  thisSubdomain->end_communication   = NULL;
  thisSubdomain->start_io            = NULL;
  thisSubdomain->start_end           = NULL;
  thisSubdomain->do_internal_RHS     = NULL;
  thisSubdomain->do_external_RHS     = NULL;
  thisSubdomain->update_fields       = NULL;
}

void
bfam_subdomain_free(bfam_subdomain_t *thisSubdomain)
{
  bfam_free(thisSubdomain->name);
  thisSubdomain->comm = NULL;
  thisSubdomain->hasWork = (0!=0);
  thisSubdomain->start_communication = NULL;
  thisSubdomain->end_communication   = NULL;
  thisSubdomain->start_io            = NULL;
  thisSubdomain->start_end           = NULL;
  thisSubdomain->do_internal_RHS     = NULL;
  thisSubdomain->do_external_RHS     = NULL;
  thisSubdomain->update_fields       = NULL;
}
