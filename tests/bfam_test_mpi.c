#include <bfam.h>

int
main (int argc, char *argv[])
{
  int rank, size;
  int source, count;
  int buffer;
  const int val = 9;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  source=0;
  if(rank == source)
  {
    buffer=val;
  }

  count=1;
  MPI_Bcast(&buffer, count, MPI_INT, source, MPI_COMM_WORLD);

  MPI_Finalize();

  if(buffer==val)
  {
    return EXIT_SUCCESS;
  }
  else
  {
    return EXIT_FAILURE;
  }
}
