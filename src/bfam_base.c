#include <bfam_base.h>

void
bfam_abort_verbose(const char *file, int line, const char *note)
{
  printf("Abort: [%s:%d] %s\n", file, line, note);
  bfam_abort();
}

void
bfam_abort()
{
  abort();
}
