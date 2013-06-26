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

  /* set up the block information */
  bfam_long_real_t x0[4] = {0,0.25,0,0.25};
  bfam_long_real_t y0[4] = {0,0,0.5,0.25};

  bfam_long_real_t x1[4] = {0.25,1,0.25,0.5};
  bfam_long_real_t y1[4] = {0,0,0.25,0.5};

  bfam_long_real_t x2[4] = {0.25,0.5,0,0};
  bfam_long_real_t y2[4] = {0.25,0.5,0.5,1};

  /* figure out which part of the domain I handle */
  if(size < 3)
  {
  }
  else
  {
  }

  /* clean up */
  bfam_domain_free(&domain);
  BFAM_MPI_CHECK(MPI_Finalize());
}
