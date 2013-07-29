#include <bfam.h>
#include <bfam_subdomain_dummy.h>

void
test_lsrk(bfam_ts_lsrk_method_t type,int N,MPI_Comm mpicomm, int commtag)
{
  bfam_ts_lsrk_t* ts;

  bfam_domain_t* dom = bfam_domain_new(mpicomm);
  bfam_subdomain_dummy_t* subDom = NULL;

  subDom = bfam_subdomain_dummy_new(1,"a",N);
  bfam_subdomain_add_tag((bfam_subdomain_t*)subDom,"yes");
  bfam_subdomain_add_tag((bfam_subdomain_t*)subDom,"other");
  bfam_subdomain_add_tag((bfam_subdomain_t*)subDom,"dummy");
  bfam_domain_add_subdomain(dom, (bfam_subdomain_t*) subDom);

  subDom = bfam_subdomain_dummy_new(2,"b",N);
  bfam_subdomain_add_tag((bfam_subdomain_t*)subDom,"no");
  bfam_subdomain_add_tag((bfam_subdomain_t*)subDom,"other");
  bfam_domain_add_subdomain(dom, (bfam_subdomain_t*) subDom);

  subDom = bfam_subdomain_dummy_new(3,"c",N);
  bfam_subdomain_add_tag((bfam_subdomain_t*)subDom,"dummy");
  bfam_subdomain_add_tag((bfam_subdomain_t*)subDom,"yes");
  bfam_subdomain_add_tag((bfam_subdomain_t*)subDom,"random");
  bfam_domain_add_subdomain(dom, (bfam_subdomain_t*) subDom);

  subDom = bfam_subdomain_dummy_new(4,"d",N);
  bfam_subdomain_add_tag((bfam_subdomain_t*)subDom,"other");
  bfam_subdomain_add_tag((bfam_subdomain_t*)subDom,"yes");
  bfam_subdomain_add_tag((bfam_subdomain_t*)subDom,"random");
  bfam_domain_add_subdomain(dom, (bfam_subdomain_t*) subDom);

  const char *tags[] = {NULL};
  bfam_communicator_t *comm =
    bfam_communicator_new(dom,BFAM_DOMAIN_AND,tags,mpicomm,commtag);


  ts = bfam_ts_lsrk_new(dom,comm,type);

  // free
  bfam_ts_lsrk_free(ts);
  bfam_free(ts);

  bfam_communicator_free(comm);
  bfam_free(comm);

  bfam_domain_free(dom);
  bfam_free(dom);
}

int
main (int argc, char *argv[])
{
  BFAM_MPI_CHECK(MPI_Init(&argc,&argv));
  test_lsrk(BFAM_TS_LSRK_KC54,4,MPI_COMM_WORLD,10);
  test_lsrk(BFAM_TS_LSRK_W33 ,3,MPI_COMM_WORLD,10);
  test_lsrk(BFAM_TS_LSRK_HEUN,2,MPI_COMM_WORLD,10);
  test_lsrk(BFAM_TS_LSRK_FE  ,1,MPI_COMM_WORLD,10);
  BFAM_MPI_CHECK(MPI_Finalize());

  return EXIT_SUCCESS;
}
