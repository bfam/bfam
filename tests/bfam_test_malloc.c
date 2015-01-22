#include <bfam.h>

int main(int argc, char *argv[])
{
  char *r1 = bfam_malloc_aligned(1024);
  char *r2 = bfam_malloc_aligned(3432);
  char *r3 = bfam_malloc_aligned(1293);
  char *r4 = bfam_malloc_aligned(12);
  char *r5 = bfam_malloc_aligned(1);
  char *r6 = bfam_malloc_aligned(0);
  char *r7 = bfam_malloc_aligned(1209);
  char *r8 = bfam_malloc_aligned(9999);

  bfam_free_aligned(r1);
  bfam_free_aligned(r2);
  bfam_free_aligned(r3);
  bfam_free_aligned(r4);
  bfam_free_aligned(r5);
  bfam_free_aligned(r6);
  bfam_free_aligned(r7);
  bfam_free_aligned(r8);

  return EXIT_SUCCESS;
}
