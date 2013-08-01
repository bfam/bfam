#include <bfam_subdomain_dummy.h>

static int
bfam_subdomain_dummy_field_add(bfam_subdomain_t *subdomain, const char *name)
{
  bfam_subdomain_dummy_t *s = (bfam_subdomain_dummy_t*) subdomain;

  if(bfam_dictionary_get_value_ptr(&s->base.fields,name))
    return 1;

  bfam_real_t *field = bfam_malloc_aligned(sizeof(bfam_real_t));

  int rval = bfam_dictionary_insert_ptr(&s->base.fields, name, field);

  BFAM_ASSERT(rval != 1);

  if(rval == 0)
    bfam_free_aligned(field);

  return rval;
}


bfam_subdomain_dummy_t*
bfam_subdomain_dummy_new(const bfam_locidx_t     id,
                         const char             *name,
                         const int               N)
{
  bfam_subdomain_dummy_t* newSubdomain =
    bfam_malloc(sizeof(bfam_subdomain_dummy_t));

  bfam_subdomain_dummy_init(newSubdomain, id, name, N);

  return newSubdomain;
}

void
bfam_subdomain_dummy_init(bfam_subdomain_dummy_t       *subdomain,
                             const bfam_locidx_t        id,
                             const char                *name,
                             const int                  N)
{
  bfam_subdomain_init(&subdomain->base, id, name);
  bfam_subdomain_add_tag(&subdomain->base, "_subdomain_dummy");
  subdomain->base.field_add = bfam_subdomain_dummy_field_add;
  subdomain->base.free = bfam_subdomain_dummy_free;
  subdomain->N = N;
}

void
bfam_subdomain_dummy_free(bfam_subdomain_t *thisSubdomain)
{
  bfam_subdomain_free(thisSubdomain);
}

bfam_long_real_t
bfam_subdomain_dummy_exact(bfam_subdomain_dummy_t *subdom, bfam_long_real_t t)
{
  return pow(t,(bfam_long_real_t) subdom->N);
}
