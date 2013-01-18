#include <bfam.h>

int
main (int argc, char *argv[])
{
  int rank;

  BFAM_CHECK_MPI(MPI_Init(&argc,&argv));
  BFAM_CHECK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &rank));

  bfam_log_init(rank, stdout, BFAM_LL_DEFAULT);

  BFAM_ROOT_INFO("BFAM Version: %s", bfam_version_get());


  BFAM_ROOT_TRACE("TRACE");
  BFAM_ROOT_LDEBUG("LDEBUG");
  BFAM_ROOT_VERBOSE("VERBOSE");
  BFAM_ROOT_INFO("INFO");
  BFAM_ROOT_WARNING("WARNING");
  BFAM_ROOT_LERROR("LERROR");

  BFAM_TRACE("TRACE");
  BFAM_LDEBUG("LDEBUG");
  BFAM_VERBOSE("VERBOSE");
  BFAM_INFO("INFO");
  BFAM_WARNING("WARNING");
  BFAM_LERROR("LERROR");

  BFAM_ROOT_TRACE(   "TRACE   %09d", 1);
  BFAM_ROOT_LDEBUG(  "LDEBUG  %09d", 2);
  BFAM_ROOT_VERBOSE( "VERBOSE %09d", 3);
  BFAM_ROOT_INFO(    "INFO    %09d", 4);
  BFAM_ROOT_WARNING( "WARNING %09d", 7);
  BFAM_ROOT_LERROR(  "LERROR  %09d", 8);

  BFAM_TRACE(   "TRACE   %09d", 1);
  BFAM_LDEBUG(  "LDEBUG  %09d", 2);
  BFAM_VERBOSE( "VERBOSE %09d", 3);
  BFAM_INFO(    "INFO    %09d", 4);
  BFAM_WARNING( "WARNING %09d", 7);
  BFAM_LERROR(  "LERROR  %09d", 8);

  BFAM_CHECK_MPI(MPI_Finalize());

  return EXIT_SUCCESS;
}
