#include <bfam.h>

int
main (int argc, char *argv[])
{
  int rank, size;
  int source, count;
  int buffer;
  const int val = 9;

  BFAM_CHECK_MPI(MPI_Init(&argc,&argv));
  BFAM_CHECK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &size));
  BFAM_CHECK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &rank));
  source=0;
  if(rank == source)
  {
    buffer=val;
  }

  count=1;
  BFAM_CHECK_MPI(MPI_Bcast(&buffer, count, MPI_INT, source, MPI_COMM_WORLD));

  BFAM_CHECK_MPI(MPI_Finalize());

  if(buffer==val)
  {
    return EXIT_SUCCESS;
  }
  else
  {
    return EXIT_FAILURE;
  }
}
