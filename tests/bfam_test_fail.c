#include <bfam.h>

int main(int argc, char *argv[])
{
  int rank;

  BFAM_MPI_CHECK(MPI_Init(&argc, &argv));
  BFAM_MPI_CHECK(MPI_Comm_rank(MPI_COMM_WORLD, &rank));

  bfam_log_init(rank, stdout, BFAM_LL_DEFAULT);

  BFAM_INFO("My rank is %3d", rank);

  BFAM_MPI_CHECK(MPI_Finalize());
  return EXIT_FAILURE;
}
