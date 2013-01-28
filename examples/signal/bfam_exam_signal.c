#include <bfam.h>

// static int
// divide_by_zero()
// {
//   int a = 1;
//   int b = 0;
//   return a / b;
// }

static void
cause_segfault()
{
  int * p = (int*)0x12345678;
  *p = 0;
}

// static void stack_overflow();
//
// static void
// stack_overflow()
// {
//   int foo[1000];
//   (void)foo;
//   stack_overflow();
// }
//
// /* break out with ctrl+c to test SIGINT handling */
// static void
// infinite_loop()
// {
//   while(1) {};
// }
//
// static void
// illegal_instruction()
// {
//   /* I couldn't find an easy way to cause this one, so I'm cheating */
//   raise(SIGILL);
// }

static void
cause_calamity()
{
  /* uncomment one of the following error conditions to cause a calamity of 
   your choosing! */

  // (void)divide_by_zero();
  cause_segfault();
  // assert(false);
  // infinite_loop();
  // illegal_instruction();
  // stack_overflow();
}

int
main (int argc, char *argv[])
{
  int rank;

  BFAM_MPI_CHECK(MPI_Init(&argc,&argv));
  BFAM_MPI_CHECK(MPI_Comm_rank(MPI_COMM_WORLD, &rank));

  bfam_log_init(rank, stdout, BFAM_LL_DEFAULT);
  bfam_signal_handler_set();

  cause_calamity();

  BFAM_WARNING("Nothing bad happend!");

  BFAM_MPI_CHECK(MPI_Finalize());

  return EXIT_SUCCESS;
}

