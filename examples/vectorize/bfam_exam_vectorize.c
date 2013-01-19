#include <bfam.h>

#define SIZE (1L << 16)

void test1(double *a, double *b)
{
  int i;

  for (i = 0; i < SIZE; i++)
  {
    a[i] += b[i];
  }
}

void test4(double * restrict a, double * restrict b)
{
  int i;

  BFAM_ASSUME_ALIGNED(a, 32);
  BFAM_ASSUME_ALIGNED(b, 32);

  for (i = 0; i < SIZE; i++)
  {
    a[i] += b[i];
  }
}

int main (int argc, char *argv[])
{
  double *a = bfam_malloc_aligned(SIZE*sizeof(double));
  double *b = bfam_malloc_aligned(SIZE*sizeof(double));

  test4(a, b);

  BFAM_INFO("a[0] = %g", a[0]);

  bfam_free_aligned(a);
  bfam_free_aligned(b);
  return EXIT_SUCCESS;
}
