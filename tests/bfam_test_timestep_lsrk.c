#include <bfam.h>
#include <bfam_subdomain_dummy.h>

void
test_lsrk(bfam_ts_lsrk_method_t type,int N)
{
  bfam_ts_lsrk_t* ts;

  bfam_domain_t* dom = bfam_domain_new(NULL);
  bfam_subdomain_dummy_t* subDom = NULL;

  subDom = bfam_subdomain_dummy_new("a",N);
  bfam_subdomain_add_tag((bfam_subdomain_t*)subDom,"yes");
  bfam_subdomain_add_tag((bfam_subdomain_t*)subDom,"other");
  bfam_subdomain_add_tag((bfam_subdomain_t*)subDom,"dummy");
  bfam_domain_add_subdomain(dom, (bfam_subdomain_t*) subDom);

  subDom = bfam_subdomain_dummy_new("b",N);
  bfam_subdomain_add_tag((bfam_subdomain_t*)subDom,"no");
  bfam_subdomain_add_tag((bfam_subdomain_t*)subDom,"other");
  bfam_domain_add_subdomain(dom, (bfam_subdomain_t*) subDom);

  subDom = bfam_subdomain_dummy_new("c",N);
  bfam_subdomain_add_tag((bfam_subdomain_t*)subDom,"dummy");
  bfam_subdomain_add_tag((bfam_subdomain_t*)subDom,"yes");
  bfam_subdomain_add_tag((bfam_subdomain_t*)subDom,"random");
  bfam_domain_add_subdomain(dom, (bfam_subdomain_t*) subDom);

  subDom = bfam_subdomain_dummy_new("d",N);
  bfam_subdomain_add_tag((bfam_subdomain_t*)subDom,"other");
  bfam_subdomain_add_tag((bfam_subdomain_t*)subDom,"yes");
  bfam_subdomain_add_tag((bfam_subdomain_t*)subDom,"random");
  bfam_domain_add_subdomain(dom, (bfam_subdomain_t*) subDom);

  ts = bfam_ts_lsrk_new(dom,type);

  // free
  bfam_ts_lsrk_free(ts);
  bfam_free(ts);

  bfam_domain_free(dom);
  bfam_free(dom);
}

int
main (int argc, char *argv[])
{
  test_lsrk(BFAM_TS_LSRK_KC54,4);
  test_lsrk(BFAM_TS_LSRK_W33 ,3);
  test_lsrk(BFAM_TS_LSRK_HEUN,2);
  test_lsrk(BFAM_TS_LSRK_FE  ,1);

  return EXIT_SUCCESS;
}
