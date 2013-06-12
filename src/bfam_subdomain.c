#include <bfam_subdomain.h>
#include <bfam_base.h>


void
bfam_subdomain_init(bfam_subdomain_t *thisSubdomain, bfam_locidx_t id,
    const char* name)
{
  thisSubdomain->id = id;
  int len = strlen(name);
  thisSubdomain->name = bfam_malloc((len+1)*sizeof(char));
  strncpy(thisSubdomain->name,name,len+1);

  bfam_dictionary_init(&thisSubdomain->fields);

  thisSubdomain->comm = NULL;

  thisSubdomain->tags.root = NULL;

  thisSubdomain->free                = bfam_subdomain_free;

  thisSubdomain->vtk_write_vtu_piece = NULL;
  thisSubdomain->field_add           = NULL;
  thisSubdomain->field_init          = NULL;
}

void
bfam_subdomain_free(bfam_subdomain_t *thisSubdomain)
{
  thisSubdomain->id = -1;
  bfam_free(thisSubdomain->name);
  bfam_critbit0_clear(&thisSubdomain->tags);
  bfam_dictionary_clear(&thisSubdomain->fields);
  thisSubdomain->comm = NULL;
  thisSubdomain->tags.root = NULL;

  thisSubdomain->vtk_write_vtu_piece = NULL;
  thisSubdomain->field_add           = NULL;
  thisSubdomain->field_init          = NULL;
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

int
bfam_subdomain_field_add(bfam_subdomain_t *thisSubdomain, const char* name)
{
  if(thisSubdomain->field_add)
    return thisSubdomain->field_add(thisSubdomain, name);
  else
    return 0;
}

void
bfam_subdomain_field_init(bfam_subdomain_t *thisSubdomain, const char* name,
      bfam_real_t time, bfam_subdomain_init_field_t init_field, void *arg)
{
  if(thisSubdomain->field_init)
    thisSubdomain->field_init(thisSubdomain, name, time, init_field, arg);
  else
    BFAM_ABORT("Field Init Not Implemented for %s!!", thisSubdomain->name);
}
