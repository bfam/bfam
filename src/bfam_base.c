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

void *
bfam_malloc(size_t size)
{
  void *r;

  r = malloc(size);

  if (!r && !size)
  {
    r = malloc(1);
  }

  if(!r)
  {
    BFAM_LERROR("Memory allocation of %lu bytes failed", (unsigned long)size);
    BFAM_ABORT("Failed malloc");
  }

#ifdef BFAM_DEBUG
  memset(r, 0xA3, size);
#endif

  return r;
}

void *
bfam_calloc(size_t nmemb, size_t size)
{
  void *r;

  r = calloc(nmemb, size);

  if (!r && (!nmemb || !size))
  {
    r = calloc(1, 1);
  }

  if(!r)
  {
    BFAM_LERROR("Memory allocation of %lu elements of size %lu bytes failed",
        (unsigned long)nmemb, (unsigned long)size);
    BFAM_ABORT("Failed calloc");
  }

  return r;
}

void *
bfam_realloc(void *ptr, size_t size)
{
  void *r;

  r = realloc(ptr, size);

  if (!r && !size)
  {
    r = realloc(ptr, 1);
  }

  if (!r)
  {
    BFAM_LERROR("Memory reallocation of size %lu bytes failed",
        (unsigned long)size);
    BFAM_ABORT("Failed realloc");
  }

  return r;
}


void
bfam_free(void *ptr)
{
  free(ptr);
}
