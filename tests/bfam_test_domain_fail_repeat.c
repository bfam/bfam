#include <bfam.h>

static void
test_insert()
{
  // domain test
  bfam_domain_t* domain;

  // Test new domain
  domain = bfam_domain_new(MPI_COMM_WORLD);

  static const char *elems[] =
      {"a", "aa", "b", "bb", "ab", "ba", "aba", "bab",
        "a1", "a2", "a3", "a4", "a5", "a6", "a7",
        "a9", "a10", "a11", "a12", "a13", "a14", "a15",
        "bab",NULL};

  for (unsigned i = 0; elems[i]; ++i)
  {
    bfam_subdomain_t* newSub = bfam_malloc(sizeof(bfam_subdomain_t));
    bfam_subdomain_init(newSub,i,-1,elems[i]);
    bfam_domain_add_subdomain(domain,newSub);
  }

  // free a domain
  bfam_domain_free(domain);
  bfam_free(domain);
}

int
main (int argc, char *argv[])
{
  test_insert();

  return EXIT_SUCCESS;
}
