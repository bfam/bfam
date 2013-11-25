#include <bfam.h>

static void
test_create()
{
  // domain test
  bfam_domain_t  domain1;
  bfam_domain_t* domain2;

  // Test init domain
  bfam_domain_init(&domain1,MPI_COMM_WORLD);

  // Test new domain
  domain2 = bfam_domain_new(MPI_COMM_WORLD);

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
  domain = bfam_domain_new(MPI_COMM_WORLD);

  static const char *elems[] =
      {"a", "aa", "b", "bb", "ab", "ba", "aba", "bab",
        "a1", "a2", "a3", "a4", "a5", "a6", "a7",
        "a9", "a10", "a11", "a12", "a13", "a14", "a15",
        NULL};

  for (unsigned i = 0; elems[i]; ++i)
  {
    bfam_subdomain_t* newSub = bfam_malloc(sizeof(bfam_subdomain_t));
    bfam_subdomain_init(newSub,i,elems[i]);
    bfam_domain_add_subdomain(domain,newSub);

    BFAM_ABORT_IF_NOT(newSub->id == (bfam_locidx_t)i,
        "Id number not set correctly");

    bfam_subdomain_add_tag(newSub, "testing 1 2 3");
    bfam_subdomain_add_tag(newSub, elems[i]);
    bfam_subdomain_add_tag(newSub, "try");

    BFAM_ABORT_IF(bfam_subdomain_has_tag(newSub, "bad tag"),
        "Error finding tag: bad tag");
    BFAM_ABORT_IF_NOT(bfam_subdomain_has_tag(newSub, "try"),
        "Error finding tag: try");
    BFAM_ABORT_IF_NOT(bfam_subdomain_has_tag(newSub, elems[i]),
        "Error finding tag: %s", elems[i]);

    BFAM_ABORT_IF_NOT(bfam_subdomain_delete_tag(newSub, "try"),
        "Error removing tag: try");
    BFAM_ABORT_IF(bfam_subdomain_has_tag(newSub, "try"),
        "Error finding tag: try");
  }

  bfam_subdomain_t **matchedSubdomains =
    bfam_malloc(domain->numSubdomains*sizeof(bfam_subdomain_t*));
  bfam_locidx_t numMatchedSubdomains;

  const char *tags1[] = {"a11", "a10", "b", NULL};
  bfam_domain_get_subdomains(domain, BFAM_DOMAIN_OR, tags1,
    domain->numSubdomains, matchedSubdomains, &numMatchedSubdomains);
  BFAM_ABORT_IF(numMatchedSubdomains != 3, "Error matching tags1: %jd",
      (intmax_t) numMatchedSubdomains);
  BFAM_ABORT_IF(matchedSubdomains[0]->id > matchedSubdomains[1]->id,
      "Error matching tags1 order");
  BFAM_ABORT_IF(matchedSubdomains[1]->id > matchedSubdomains[2]->id,
      "Error matching tags1 order");

  const char *tags2[] = {"a11", "testing 1 2 3", NULL};
  bfam_domain_get_subdomains(domain, BFAM_DOMAIN_AND, tags2,
    domain->numSubdomains, matchedSubdomains, &numMatchedSubdomains);
  BFAM_ABORT_IF(numMatchedSubdomains != 1, "Error matching tags2: %jd",
      (intmax_t) numMatchedSubdomains);

  bfam_critbit0_tree_t ctags1 = {0};
  bfam_critbit0_insert(&ctags1, "a11");
  bfam_critbit0_insert(&ctags1, "a10");
  bfam_critbit0_insert(&ctags1, "b");
  bfam_domain_get_subdomains_critbit(domain, BFAM_DOMAIN_OR, &ctags1,
    domain->numSubdomains, matchedSubdomains, &numMatchedSubdomains);
  BFAM_ABORT_IF(numMatchedSubdomains != 3, "Error matching ctags1: %jd",
      (intmax_t) numMatchedSubdomains);
  BFAM_ABORT_IF(matchedSubdomains[0]->id > matchedSubdomains[1]->id,
      "Error matching ctags1 order");
  BFAM_ABORT_IF(matchedSubdomains[1]->id > matchedSubdomains[2]->id,
      "Error matching ctags1 order");
  bfam_critbit0_clear(&ctags1);

  bfam_critbit0_tree_t ctags2 = {0};
  bfam_critbit0_insert(&ctags2, "a11");
  bfam_critbit0_insert(&ctags2, "testing 1 2 3");
  bfam_domain_get_subdomains_critbit(domain, BFAM_DOMAIN_AND, &ctags2,
    domain->numSubdomains, matchedSubdomains, &numMatchedSubdomains);
  BFAM_ABORT_IF(numMatchedSubdomains != 1, "Error matching ctags2: %jd",
      (intmax_t) numMatchedSubdomains);
  bfam_critbit0_clear(&ctags2);

  bfam_free(matchedSubdomains);
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
