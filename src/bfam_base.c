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

static size_t
bfam_page_size()
{
  long page_size = sysconf(_SC_PAGE_SIZE);

  return (size_t)page_size;
}

#if defined(__APPLE__)

static size_t
bfam_cache_line_size()
{
  size_t line_size = 0;
  size_t sizeof_line_size = sizeof(line_size);
  sysctlbyname("hw.cachelinesize", &line_size, &sizeof_line_size, 0, 0);
  return line_size;
}

#elif defined(_WIN32)

static size_t
bfam_cache_line_size() {
  size_t line_size = 0;
  DWORD buffer_size = 0;
  DWORD i = 0;
  SYSTEM_LOGICAL_PROCESSOR_INFORMATION * buffer = 0;

  GetLogicalProcessorInformation(0, &buffer_size);
  buffer = (SYSTEM_LOGICAL_PROCESSOR_INFORMATION *)bfam_malloc(buffer_size);
  GetLogicalProcessorInformation(&buffer[0], &buffer_size);

  for (i=0; i != buffer_size/sizeof(SYSTEM_LOGICAL_PROCESSOR_INFORMATION); ++i)
  {
    if (buffer[i].Relationship == RelationCache && buffer[i].Cache.Level == 1)
    {
      line_size = buffer[i].Cache.LineSize;
      break;
    }
  }

  bfam_free(buffer);
  return line_size;
}

#elif defined(__linux__)

static size_t
bfam_cache_line_size()
{
  long line_size = sysconf(_SC_LEVEL1_DCACHE_LINESIZE);

  return (size_t)line_size;
}

#else
#error Unrecognized platform for cache line size
#endif

static void *
bfam_malloc_aligned_cache_line(size_t size, size_t line,
   size_t line_size, size_t page_size)
{
  void *r;
  intptr_t a;

  BFAM_ASSERT(page_size >= 1);
  BFAM_ASSERT(page_size >= line*line_size);

  r = bfam_malloc(size +  sizeof(intptr_t) + page_size);

  a = (((intptr_t)r + sizeof(intptr_t) + page_size - line*line_size - 1)/
      page_size)*page_size + line*line_size;

  ((intptr_t *)a)[-1] = (intptr_t)r;

  r = (void *)a;

  return r;
}

void *
bfam_malloc_aligned(size_t size)
{
  void * r;
  static size_t line_no = 0;

  const size_t line_size = bfam_cache_line_size();
  const size_t page_size = bfam_page_size();
  const size_t line_count = page_size/line_size;

  r = bfam_malloc_aligned_cache_line(size, line_no, line_size, page_size);

  BFAM_ABORT_IF_NOT(BFAM_IS_ALIGNED(r,16), "Memory not 16 bit aligned");
  BFAM_ABORT_IF_NOT(BFAM_IS_ALIGNED(r,32), "Memory not 32 bit aligned");
  BFAM_ABORT_IF_NOT(BFAM_IS_ALIGNED(r,64), "Memory not 64 bit aligned");

  line_no = (line_no + 1) % line_count;

#ifdef BFAM_DEBUG
  memset(r, 0xA3, size);
#endif

  return r;
}

void
bfam_free_aligned(void *ptr)
{
  ptr = (void *) ((intptr_t *) ptr)[-1];
  BFAM_ASSERT(ptr != NULL);
  free(ptr);
}
