#include <bfam.h>

static void
test_create()
{
  // domain test
  bfam_domain_t  domain1;
  bfam_domain_t* domain2;

  // Test init domain
  bfam_domain_init(&domain1,NULL);

  // Test new domain
  domain2 = bfam_domain_new(NULL);

  // free a domain
  bfam_domain_free(&domain1);
  bfam_domain_free( domain2);
  bfam_free(domain2);
}

static void
test_insert()
{
  // domain test
  bfam_domain_t* domain;

  // Test new domain
  domain = bfam_domain_new(NULL);

  static const char *elems[] =
      {"a", "aa", "b", "bb", "ab", "ba", "aba", "bab", NULL};

  for (unsigned i = 0; elems[i]; ++i)
  {
    printf("Inserting %s\n",elems[i]);
    bfam_subdomain_t* newSub = bfam_malloc(sizeof(bfam_subdomain_t));
    bfam_subdomain_init(newSub,elems[i]);
    bfam_domain_add_subdomain(domain,newSub);
  }

  // free a domain
  bfam_domain_free(domain);
  bfam_free(domain);
}

int
main (int argc, char *argv[])
{
  test_create();
  test_insert();

  return EXIT_SUCCESS;
}
