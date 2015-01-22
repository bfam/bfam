#include <bfam.h>

int main(int argc, char *argv[])
{

#define LIST_OF_TESTS                                                          \
  X(8, double, 1024)                                                           \
  X(16, double, 1024)                                                          \
  X(32, double, 1024)                                                          \
  X(64, double, 1024)                                                          \
  X(128, double, 1024)                                                         \
  X(8, char, 1024)                                                             \
  X(16, char, 1024)                                                            \
  X(32, char, 1024)                                                            \
  X(64, char, 1024)                                                            \
  X(128, char, 1024)

#define X(n, type, size)                                                       \
  {                                                                            \
    BFAM_INFO("Testing " #n " : " #type " : " #size);                          \
    BFAM_ALIGN(n) type tmp[size];                                              \
    BFAM_ABORT_IF_NOT(BFAM_IS_ALIGNED(tmp, n),                                 \
                      "Memory " #n " " #type " " #size " is not aligned: %p",  \
                      tmp);                                                    \
  }
  LIST_OF_TESTS
#undef X

  return EXIT_SUCCESS;
}
