#include <bfam.h>

int
main (int argc, char *argv[])
{
  int rank, size;

  /* startup MPI */
  BFAM_MPI_CHECK(MPI_Init(&argc,&argv));
  BFAM_MPI_CHECK(MPI_Comm_size(MPI_COMM_WORLD, &size));
  BFAM_MPI_CHECK(MPI_Comm_rank(MPI_COMM_WORLD, &rank));
  bfam_log_init(rank,stdout,BFAM_LL_VERBOSE);

  /* setup the domain*/
  bfam_domain_t domain;
  bfam_domain_init(&domain,MPI_COMM_WORLD);


  /* clean up */
  bfam_domain_free(&domain);
  BFAM_MPI_CHECK(MPI_Finalize());
}
