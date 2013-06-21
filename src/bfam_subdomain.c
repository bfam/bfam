#include <bfam_subdomain.h>
#include <bfam_base.h>

/*
 * Compare function which sorts a bfam_subdomain_face_map_entry_t array
 * in sending order.
 */
int
bfam_subdomain_face_send_cmp(const void* a, const void*b)
{
  const bfam_subdomain_face_map_entry_t *la = a;
  const bfam_subdomain_face_map_entry_t *lb = b;

  if(la->np < lb->np)
    return -1;
  else if(la->np > lb->np)
    return  1;
  else if(la->ns < lb->ns)
    return -1;
  else if(la->ns > lb->ns)
    return  1;
  else if(la->s  < lb->s)
    return -1;
  else if(la->s  > lb->s)
    return  1;
  else if(la->nk < lb->nk)
    return -1;
  else if(la->nk > lb->nk)
    return  1;
  else if(la->nf < lb->nf)
    return -1;
  else if(la->nf > lb->nf)
    return  1;
  else if(la->nh < lb->nh)
    return -1;
  else if(la->nh > lb->nh)
    return  1;
  else
    return  0;

  BFAM_ABORT("We should never reach here.");
}

/*
 * Compare function which sorts a bfam_subdomain_face_map_entry_t array
 * in receiving order.
 */
int
bfam_subdomain_face_recv_cmp(const void* a, const void*b)
{
  const bfam_subdomain_face_map_entry_t *la = a;
  const bfam_subdomain_face_map_entry_t *lb = b;

  if(la->np < lb->np)
    return -1;
  else if(la->np > lb->np)
    return  1;
  else if(la->s  < lb->s)
    return -1;
  else if(la->s  > lb->s)
    return  1;
  else if(la->ns < lb->ns)
    return -1;
  else if(la->ns > lb->ns)
    return  1;
  else if(la->k  < lb->k)
    return -1;
  else if(la->k  > lb->k)
    return  1;
  else if(la->f  < lb->f)
    return -1;
  else if(la->f  > lb->f)
    return  1;
  else if(la->h  < lb->h)
    return -1;
  else if(la->h  > lb->h)
    return  1;
  else
    return  0;

  BFAM_ABORT("We should never reach here.");
}

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
  thisSubdomain->field_plus_add      = NULL;
  thisSubdomain->field_minus_add     = NULL;
  thisSubdomain->field_init          = NULL;

  thisSubdomain->glue_comm_info      = NULL;
  thisSubdomain->glue_put_send_buffer = NULL;
  thisSubdomain->glue_get_recv_buffer = NULL;
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

  thisSubdomain->vtk_write_vtu_piece  = NULL;
  thisSubdomain->field_add            = NULL;
  thisSubdomain->field_plus_add       = NULL;
  thisSubdomain->field_minus_add      = NULL;
  thisSubdomain->field_init           = NULL;

  thisSubdomain->glue_comm_info       = NULL;
  thisSubdomain->glue_put_send_buffer = NULL;
  thisSubdomain->glue_get_recv_buffer = NULL;
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

int
bfam_subdomain_field_plus_add(bfam_subdomain_t *thisSubdomain,
    const char* name)
{
  if(thisSubdomain->field_plus_add)
    return thisSubdomain->field_plus_add(thisSubdomain, name);
  else
    return 0;
}

int
bfam_subdomain_field_minus_add(bfam_subdomain_t *thisSubdomain,
    const char* name)
{
  if(thisSubdomain->field_minus_add)
    return thisSubdomain->field_minus_add(thisSubdomain, name);
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
