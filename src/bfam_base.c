#include <bfam_base.h>
#include <bfam_log.h>

void
bfam_abort_verbose(const char *file, int line, ...)
{
  va_list ap;
  const char * fmt;
  char tmpStr[256];
  va_start(ap,line);
  fmt = va_arg(ap, const char *);
  sprintf(tmpStr,"Abort: [%s:%d] %s\n", file, line, fmt);
  BFAM_LERROR(tmpStr,ap);
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

/*
 * This signal handler code is a modified version of the one presented in
 *
 *     http://spin.atomicobject.com/2013/01/13/exceptions-stack-traces-c/
 *
 * by Job Vranish.  The code can be found here
 *
 *     https://gist.github.com/4441299
 */

#define BFAM_MAX_STACK_FRAMES 1024
static void *stack_traces[BFAM_MAX_STACK_FRAMES];

static void bfam_posix_print_stack_trace()
{
  int i, trace_size = 0;
  char **messages = (char **)NULL;

  trace_size = backtrace(stack_traces, BFAM_MAX_STACK_FRAMES);
  messages = backtrace_symbols(stack_traces, trace_size);
  BFAM_SYS_ERROR_CHECK(messages == NULL, "backtrace_symbols");

  for(i = 0; i < trace_size; ++i)
  {
    BFAM_LERROR("%s", messages[i]);
  }

  if(messages)
  {
    free(messages);
  }
}

void bfam_posix_signal_handler(int sig, siginfo_t *siginfo, void *context)
{
  (void)context;
  switch(sig)
  {
    case SIGSEGV:
      BFAM_LERROR("Caught SIGSEGV: Segmentation Fault");
      break;
    case SIGINT:
      BFAM_LERROR(
          "Caught SIGINT: Interactive attention signal, (usually ctrl+c)");
      break;
    case SIGFPE:
      switch(siginfo->si_code)
      {
        case FPE_INTDIV:
          BFAM_LERROR("Caught SIGFPE: (integer divide by zero)");
          break;
        case FPE_INTOVF:
          BFAM_LERROR("Caught SIGFPE: (integer overflow)");
          break;
        case FPE_FLTDIV:
          BFAM_LERROR("Caught SIGFPE: (floating-point divide by zero)");
          break;
        case FPE_FLTOVF:
          BFAM_LERROR("Caught SIGFPE: (floating-point overflow)");
          break;
        case FPE_FLTUND:
          BFAM_LERROR("Caught SIGFPE: (floating-point underflow)");
          break;
        case FPE_FLTRES:
          BFAM_LERROR("Caught SIGFPE: (floating-point inexact result)");
          break;
        case FPE_FLTINV:
          BFAM_LERROR("Caught SIGFPE: (floating-point invalid operation)");
          break;
        case FPE_FLTSUB:
          BFAM_LERROR("Caught SIGFPE: (subscript out of range)");
          break;
        default:
          BFAM_LERROR("Caught SIGFPE: Arithmetic Exception");
          break;
      }
    case SIGILL:
      switch(siginfo->si_code)
      {
        case ILL_ILLOPC:
          BFAM_LERROR("Caught SIGILL: (illegal opcode)");
          break;
        case ILL_ILLOPN:
          BFAM_LERROR("Caught SIGILL: (illegal operand)");
          break;
        case ILL_ILLADR:
          BFAM_LERROR("Caught SIGILL: (illegal addressing mode)");
          break;
        case ILL_ILLTRP:
          BFAM_LERROR("Caught SIGILL: (illegal trap)");
          break;
        case ILL_PRVOPC:
          BFAM_LERROR("Caught SIGILL: (privileged opcode)");
          break;
        case ILL_PRVREG:
          BFAM_LERROR("Caught SIGILL: (privileged register)");
          break;
        case ILL_COPROC:
          BFAM_LERROR("Caught SIGILL: (coprocessor error)");
          break;
        case ILL_BADSTK:
          BFAM_LERROR("Caught SIGILL: (internal stack error)");
          break;
        default:
          BFAM_LERROR("Caught SIGILL: Illegal Instruction");
          break;
      }
      break;
    case SIGTERM:
      BFAM_LERROR(
          "Caught SIGTERM: a termination request was sent to the program");
      break;
    case SIGABRT:
      BFAM_LERROR("Caught SIGABRT: usually caused by an abort() or assert()");
      break;
    default:
      break;
  }
  bfam_posix_print_stack_trace();
  _Exit(1);
}

static uint8_t alternate_stack[SIGSTKSZ];
void bfam_signal_handler_set()
{
  /* setup alternate stack */
  {
    stack_t ss;
    /* malloc is usually used here, I'm not 100% sure my static allocation
       is valid but it seems to work just fine. */
    ss.ss_sp = (void*)alternate_stack;
    ss.ss_size = SIGSTKSZ;
    ss.ss_flags = 0;

    if(sigaltstack(&ss, NULL) != 0) { err(1, "sigaltstack"); }
  }

  /* register our signal handlers */
  {
    struct sigaction sig_action;
    sig_action.sa_sigaction = bfam_posix_signal_handler;
    sigemptyset(&sig_action.sa_mask);

#ifdef __APPLE__
    /* for some reason we backtrace() doesn't work on osx
       when we use an alternate stack */
    sig_action.sa_flags = SA_SIGINFO;
#else
    sig_action.sa_flags = SA_SIGINFO | SA_ONSTACK;
#endif

    if(sigaction(SIGSEGV, &sig_action, NULL) != 0) { err(1, "sigaction"); }
    if(sigaction(SIGFPE,  &sig_action, NULL) != 0) { err(1, "sigaction"); }
    if(sigaction(SIGINT,  &sig_action, NULL) != 0) { err(1, "sigaction"); }
    if(sigaction(SIGILL,  &sig_action, NULL) != 0) { err(1, "sigaction"); }
    if(sigaction(SIGTERM, &sig_action, NULL) != 0) { err(1, "sigaction"); }
    if(sigaction(SIGABRT, &sig_action, NULL) != 0) { err(1, "sigaction"); }
  }
}
