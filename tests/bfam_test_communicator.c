#include <bfam.h>

int
main (int argc, char *argv[])
{
  int rank, size;

  BFAM_MPI_CHECK(MPI_Init(&argc,&argv));
  BFAM_MPI_CHECK(MPI_Comm_size(MPI_COMM_WORLD, &size));
  BFAM_MPI_CHECK(MPI_Comm_rank(MPI_COMM_WORLD, &rank));


  BFAM_MPI_CHECK(MPI_Finalize());
}
