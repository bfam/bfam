#include <bfam_subdomain_dummy.h>

bfam_subdomain_dummy_t*
bfam_subdomain_dummy_new(const char             *name,
                         const int               N)
{
  bfam_subdomain_dummy_t* newSubdomain =
    bfam_malloc(sizeof(bfam_subdomain_dummy_t));

  bfam_subdomain_dummy_init(newSubdomain, name, N);

  return newSubdomain;
}

void
bfam_subdomain_dummy_init(bfam_subdomain_dummy_t       *subdomain,
                             const char                *name,
                             const int                  N)
{
  bfam_subdomain_init(&subdomain->base, name);
  bfam_subdomain_add_tag(&subdomain->base, "_subdomain_dummy");
  subdomain->base.free = bfam_subdomain_dummy_free;
  subdomain->N = N;
}

void
bfam_subdomain_dummy_free(bfam_subdomain_t *thisSubdomain)
{
  bfam_subdomain_free(thisSubdomain);
}
