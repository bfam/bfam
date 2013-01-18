#include <bfam_base.h>
#include <bfam_log.h>

void
bfam_abort_verbose(const char *file, int line, const char *note)
{
  BFAM_LERROR("Abort: [%s:%d] %s\n", file, line, note);
  bfam_abort();
}

void
bfam_abort()
{
  abort();
}
