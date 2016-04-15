#include <err.h>
#include <execinfo.h>
#include <fcntl.h>
#include <signal.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sysexits.h>

#include "bfam.h"

// {{{ base

#ifdef __GNUC__
#define BFAM_UNUSED_VAR __attribute__((unused))
#else
#define BFAM_UNUSED_VAR
#endif

#define BFAM_PXEST_DIMENSION BFAM_DGX_DIMENSION
#define DIM BFAM_PXEST_DIMENSION

#define BFAM_APPROX_EQ(x, y, K, abs, eps, min)                                 \
  ((abs)((x) - (y)) < (min) + (K) * (eps)*BFAM_MAX((abs)((x)), (abs)((y))))

#define BFAM_LITTLE_ENDIAN 0
#define BFAM_BIG_ENDIAN 1

static void bfam_abort() { abort(); }

void bfam_abort_verbose(const char *file, int line, ...)
{
  va_list ap;
  const char *fmt;
  char note[BFAM_BUFSIZ];

  va_start(ap, line);
  fmt = va_arg(ap, const char *);
  vsnprintf(note, BFAM_BUFSIZ, fmt, ap);
  va_end(ap);
  BFAM_LERROR("Abort: [%s:%d] %s", file, line, note);
  bfam_abort();
}

void *bfam_malloc(size_t size)
{
  void *r;

  r = malloc(size);

  if (!r && !size)
  {
    r = malloc(1);
  }

  if (!r)
  {
    BFAM_LERROR("Memory allocation of %lu bytes failed", (unsigned long)size);
    BFAM_ABORT("Failed malloc");
  }

#ifdef BFAM_DEBUG
  memset(r, 0xA3, size);
#endif

  return r;
}

void *bfam_calloc(size_t nmemb, size_t size)
{
  void *r;

  r = calloc(nmemb, size);

  if (!r && (!nmemb || !size))
  {
    r = calloc(1, 1);
  }

  if (!r)
  {
    BFAM_LERROR("Memory allocation of %lu elements of size %lu bytes failed",
                (unsigned long)nmemb, (unsigned long)size);
    BFAM_ABORT("Failed calloc");
  }

  return r;
}

void *bfam_realloc(void *ptr, size_t size)
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

void bfam_free(void *ptr) { free(ptr); }

static size_t bfam_page_size()
{
  long page_size = sysconf(_SC_PAGE_SIZE);

  return (size_t)page_size;
}

#if defined(__APPLE__)
#include <sys/sysctl.h>
#include <sys/types.h>

static size_t bfam_cache_line_size()
{
  size_t line_size = 0;
  size_t sizeof_line_size = sizeof(line_size);
  sysctlbyname("hw.cachelinesize", &line_size, &sizeof_line_size, 0, 0);
  return line_size;
}

#elif defined(_WIN32)

static size_t bfam_cache_line_size()
{
  size_t line_size = 0;
  DWORD buffer_size = 0;
  DWORD i = 0;
  SYSTEM_LOGICAL_PROCESSOR_INFORMATION *buffer = 0;

  GetLogicalProcessorInformation(0, &buffer_size);
  buffer = (SYSTEM_LOGICAL_PROCESSOR_INFORMATION *)bfam_malloc(buffer_size);
  GetLogicalProcessorInformation(&buffer[0], &buffer_size);

  for (i = 0; i != buffer_size / sizeof(SYSTEM_LOGICAL_PROCESSOR_INFORMATION);
       ++i)
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

static size_t bfam_cache_line_size()
{
  long line_size = sysconf(_SC_LEVEL1_DCACHE_LINESIZE);

  return (size_t)line_size;
}

#else
#error Unrecognized platform for cache line size
#endif

static void *bfam_malloc_aligned_cache_line(size_t size, size_t line,
                                            size_t line_size, size_t page_size)
{
  void *r;
  intptr_t a;

  BFAM_ASSERT(page_size >= 1);
  BFAM_ASSERT(page_size >= line * line_size);

  r = bfam_malloc(size + sizeof(intptr_t) + page_size);

  a = (((intptr_t)r + sizeof(intptr_t) + page_size - line * line_size - 1) /
       page_size) *
          page_size +
      line * line_size;

  ((intptr_t *)a)[-1] = (intptr_t)r;

  r = (void *)a;

  return r;
}

void *bfam_malloc_aligned(size_t size)
{
  void *r;
  static size_t line_no = 0;

  const size_t line_size = bfam_cache_line_size();
  const size_t page_size = bfam_page_size();
  const size_t line_count = page_size / line_size;

  r = bfam_malloc_aligned_cache_line(size, line_no, line_size, page_size);

  BFAM_ABORT_IF_NOT(BFAM_IS_ALIGNED(r, 16), "Memory not 16 bit aligned");
  BFAM_ABORT_IF_NOT(BFAM_IS_ALIGNED(r, 32), "Memory not 32 bit aligned");
  BFAM_ABORT_IF_NOT(BFAM_IS_ALIGNED(r, 64), "Memory not 64 bit aligned");

  line_no = (line_no + 1) % line_count;

#ifdef BFAM_DEBUG
  memset(r, 0xA3, size);
#endif

  return r;
}

void bfam_free_aligned(void *ptr)
{
  ptr = (void *)((intptr_t *)ptr)[-1];
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

  for (i = 0; i < trace_size; ++i)
  {
    BFAM_LERROR("%s", messages[i]);
  }

  if (messages)
  {
    free(messages);
  }
}

static void bfam_posix_signal_handler(int sig, siginfo_t *siginfo,
                                      void *context)
{
  (void)context;
  switch (sig)
  {
  case SIGSEGV:
    BFAM_LERROR("Caught SIGSEGV: Segmentation Fault");
    break;
  case SIGINT:
    BFAM_LERROR(
        "Caught SIGINT: Interactive attention signal, (usually ctrl+c)");
    break;
  case SIGFPE:
    switch (siginfo->si_code)
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
    switch (siginfo->si_code)
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
    ss.ss_sp = (void *)alternate_stack;
    ss.ss_size = SIGSTKSZ;
    ss.ss_flags = 0;

    if (sigaltstack(&ss, NULL) != 0)
    {
      err(1, "sigaltstack");
    }
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

    if (sigaction(SIGSEGV, &sig_action, NULL) != 0)
    {
      err(1, "sigaction");
    }
    if (sigaction(SIGFPE, &sig_action, NULL) != 0)
    {
      err(1, "sigaction");
    }
    if (sigaction(SIGINT, &sig_action, NULL) != 0)
    {
      err(1, "sigaction");
    }
    if (sigaction(SIGILL, &sig_action, NULL) != 0)
    {
      err(1, "sigaction");
    }
    if (sigaction(SIGTERM, &sig_action, NULL) != 0)
    {
      err(1, "sigaction");
    }
    if (sigaction(SIGABRT, &sig_action, NULL) != 0)
    {
      err(1, "sigaction");
    }
  }
}

/*
 *  The follow code snippet is from:
 *
 *    http://www.ibm.com/developerworks/aix/library/au-endianc/
 */
static int bfam_endian()
{
  int i = 1;
  char *p = (char *)&i;

  if (p[0] == 1)
    return BFAM_LITTLE_ENDIAN;
  else
    return BFAM_BIG_ENDIAN;
}

// }}}

// {{{ critbit

typedef struct
{
  void *child[2];
  uint32_t byte;
  uint8_t otherbits;
} bfam_critbit0_node_t;

/** Membership testing.
 *
 * The following function takes a tree, \a t, and a \c NULL terminated string,
 * \a u,  and returns non-zero iff \a u in \a t.
 *
 * \param [in] t tree
 * \param [in] u possible member
 * \returns non-zero iff \a u in \a t
 */
static int bfam_critbit0_contains(bfam_critbit0_tree_t *t, const char *u)
{
  const uint8_t *ubytes = (void *)u;
  const size_t ulen = strlen(u);
  uint8_t *p = t->root;

  if (!p)
    return 0;

  while (1 & (intptr_t)p)
  {
    bfam_critbit0_node_t *q = (void *)(p - 1);

    uint8_t c = 0;
    if (q->byte < ulen)
      c = ubytes[q->byte];
    const int direction = (1 + (q->otherbits | c)) >> 8;

    p = q->child[direction];
  }

  return 0 == strcmp(u, (const char *)p);
}

/** Inserting into the tree.
 *
 * It takes a tree, \a t, and possibly mutates it such that a \c NULL
 * terminated string, \a u, is a member on exit.
 *
 * \param [in,out] t tree
 * \param [in] u possible member
 * \returns:
 *   $\cases{ 0 &if {\rm out of memory} \cr
 *            1 &if {\it u} {\rm was already a member} \cr
 *           2 &if {\it t} {\rm was mutated successfully}}$.
 */
static int bfam_critbit0_insert(bfam_critbit0_tree_t *t, const char *u)
{
  const uint8_t *const ubytes = (void *)u;
  const size_t ulen = strlen(u);
  uint8_t *p = t->root;

  if (!p)
  {
    char *x;
    int a = posix_memalign((void **)&x, sizeof(void *), ulen + 1);
    if (a)
      return 0;
    memcpy(x, u, ulen + 1);
    t->root = x;
    return 2;
  }

  while (1 & (intptr_t)p)
  {
    bfam_critbit0_node_t *q = (void *)(p - 1);

    uint8_t c = 0;
    if (q->byte < ulen)
      c = ubytes[q->byte];
    const int direction = (1 + (q->otherbits | c)) >> 8;

    p = q->child[direction];
  }

  uint32_t newbyte;
  uint32_t newotherbits;

  for (newbyte = 0; newbyte < ulen; ++newbyte)
  {
    if (p[newbyte] != ubytes[newbyte])
    {
      newotherbits = p[newbyte] ^ ubytes[newbyte];
      goto different_byte_found;
    }
  }

  if (p[newbyte] != 0)
  {
    newotherbits = p[newbyte];
    goto different_byte_found;
  }
  return 1;

different_byte_found:

  newotherbits |= newotherbits >> 1;
  newotherbits |= newotherbits >> 2;
  newotherbits |= newotherbits >> 4;
  newotherbits = (newotherbits & ~(newotherbits >> 1)) ^ 255;
  uint8_t c = p[newbyte];
  int newdirection = (1 + (newotherbits | c)) >> 8;

  bfam_critbit0_node_t *newnode;

  if (posix_memalign((void **)&newnode, sizeof(void *),
                     sizeof(bfam_critbit0_node_t)))
  {
    return 0;
  }

  char *x;
  if (posix_memalign((void **)&x, sizeof(void *), ulen + 1))
  {
    free(newnode);
    return 0;
  }
  memcpy(x, ubytes, ulen + 1);

  newnode->byte = newbyte;
  newnode->otherbits = (int8_t)newotherbits;
  newnode->child[1 - newdirection] = x;

  void **wherep = &t->root;
  for (;;)
  {
    uint8_t *p = *wherep;
    if (!(1 & (intptr_t)p))
      break;
    bfam_critbit0_node_t *q = (void *)(p - 1);
    if (q->byte > newbyte)
      break;
    if (q->byte == newbyte && q->otherbits > newotherbits)
      break;
    uint8_t c = 0;
    if (q->byte < ulen)
      c = ubytes[q->byte];
    const int direction = (1 + (q->otherbits | c)) >> 8;
    wherep = q->child + direction;
  }

  newnode->child[newdirection] = *wherep;
  *wherep = (void *)(1 + (char *)newnode);

  return 2;
}

/** Deleting elements.
 *
 * This function takes a tree, \a t, and a \c NULL terminated string,
 * \a u, and possibly mutates the tree such that $u \notin t$.
 *
 * \param [in,out] t tree
 * \param [in] u possible member to remove from \a t
 * \returns It returns 1 if the tree was mutated, 0 otherwise.
 */
static int bfam_critbit0_delete(bfam_critbit0_tree_t *t, const char *u)
{
  const uint8_t *ubytes = (void *)u;
  const size_t ulen = strlen(u);
  uint8_t *p = t->root;
  void **wherep = &t->root;
  void **whereq = 0;
  bfam_critbit0_node_t *q = 0;
  int direction = 0;

  if (!p)
    return 0;

  while (1 & (intptr_t)p)
  {
    whereq = wherep;
    q = (void *)(p - 1);
    uint8_t c = 0;
    if (q->byte < ulen)
      c = ubytes[q->byte];
    direction = (1 + (q->otherbits | c)) >> 8;
    wherep = q->child + direction;
    p = *wherep;
  }

  if (0 != strcmp(u, (const char *)p))
    return 0;
  free(p);

  if (!whereq)
  {
    t->root = 0;
    return 1;
  }

  *whereq = q->child[1 - direction];
  free(q);

  return 1;
}

static void traverse(void *top)
{
  uint8_t *p = top;

  if (1 & (intptr_t)p)
  {
    bfam_critbit0_node_t *q = (void *)(p - 1);
    traverse(q->child[0]);
    traverse(q->child[1]);
    free(q);
  }
  else
  {
    free(p);
  }
}

/** Clearing a tree.
 *
 * Clearing a tree (freeing all members) brings us our first code for walking
 * the whole tree rather than just tracing a path through it.
 *
 * So, the \c critbit0_clear function takes a tree, \a t, and frees every
 * member of it, mutating the tree such that it is empty on exit.
 *
 * \param [in,out] t tree
 */
static void bfam_critbit0_clear(bfam_critbit0_tree_t *t)
{
  if (t->root)
    traverse(t->root);
  t->root = NULL;
}

static int allprefixed_traverse(uint8_t *top,
                                int (*handle)(const char *, void *), void *arg)
{

  if (1 & (intptr_t)top)
  {
    bfam_critbit0_node_t *q = (void *)(top - 1);
    for (int direction = 0; direction < 2; ++direction)
      switch (allprefixed_traverse(q->child[direction], handle, arg))
      {
      case 1:
        break;
      case 0:
        return 0;
      default:
        return -1;
      }
    return 1;
  }

  return handle((const char *)top, arg);
}

/** Fetching elements with a given prefix.
 *
 * One of the operations which crit-bit trees can perform efficiently that hash
 * tables cannot is the extraction of the subset of elements with a given
 * prefix.
 *
 * The following function takes a tree, \a t, and a \c NULL terminated string,
 * \a prefix. Let $S \subseteq t$ where $x \in S$ iff \a prefix is a prefix of
 * \c x, then $\forall x : S.$ \a handle is called with arguments \c x and
 * \c arg.
 * \returns:
 *   $\cases{ 0 &if {\it handle} {\rm returned 0} \cr
 *            1 &if {\rm successful} \cr
 *            2 &if {\it handle} {\rm returned a value} $\notin [0,1]$}$
 * \note (Note that, if |handle| returns 0, the iteration is aborted)
 */
static int bfam_critbit0_allprefixed(bfam_critbit0_tree_t *t,
                                     const char *prefix,
                                     int (*handle)(const char *, void *),
                                     void *arg)
{
  const uint8_t *ubytes = (void *)prefix;
  const size_t ulen = strlen(prefix);
  uint8_t *p = t->root;
  uint8_t *top = p;

  if (!p)
    return 1;

  while (1 & (intptr_t)p)
  {
    bfam_critbit0_node_t *q = (void *)(p - 1);
    uint8_t c = 0;
    if (q->byte < ulen)
      c = ubytes[q->byte];
    const int direction = (1 + (q->otherbits | c)) >> 8;
    p = q->child[direction];
    if (q->byte < ulen)
      top = p;
  }

  for (size_t i = 0; i < ulen; ++i)
  {
    if (p[i] != ubytes[i])
      return 1;
  }

  return allprefixed_traverse(top, handle, arg);
}

// }}}

// {{{ dictionary

#define BFAM_KEYVALUE_SPLIT '\255'
#define BFAM_PTR_STR_LEN BFAM_BUFSIZ

void bfam_dictionary_init(bfam_dictionary_t *d)
{
  d->num_entries = 0;
  d->t.root = NULL;
}

void bfam_dictionary_clear(bfam_dictionary_t *d)
{
  bfam_critbit0_clear(&(d->t));
  d->num_entries = 0;
}

static int bfam_dictionary_contains_check(const char *value, void *arg)
{
  (*(int *)arg)++;
  return 1;
}

int bfam_dictionary_contains(bfam_dictionary_t *d, const char *key)
{
  const size_t keylen = strlen(key);

  char *u = bfam_malloc(sizeof(char) * (keylen + 2));

  memcpy(u, key, keylen);
  u[keylen] = BFAM_KEYVALUE_SPLIT;
  u[keylen + 1] = '\0';
  int found = 0;
  bfam_critbit0_allprefixed(&(d->t), u, &bfam_dictionary_contains_check,
                            &found);
  bfam_free(u);
  return found;
}

/** Inserting key and value pair into a dictionary
 *
 * It takes a dictionary, \a d, and possibly mutates it such that a \c NULL
 * terminated strings, \a key and \a value, is a member on exit.
 *
 * \param [in,out] d dictionary
 * \param [in] key possible key
 * \param [in] val possible value
 * \returns:
 *   $\cases{ 0 &if {\rm out of memory} \cr
 *            1 &if {\it key} {\rm was already a member} \cr
 *            2 &if {\it d} {\rm was mutated successfully}}$.
 */
static int bfam_dictionary_insert(bfam_dictionary_t *d, const char *key,
                                  const char *val)
{
  const size_t keylen = strlen(key);
  const size_t vallen = strlen(val);

  char *keyval = bfam_malloc(sizeof(char) * (keylen + vallen + 2));

  memcpy(keyval, key, keylen);

  if (bfam_dictionary_contains(d, key))
  {
    free(keyval);
    return 1;
  }

  keyval[keylen] = BFAM_KEYVALUE_SPLIT;
  memcpy(&keyval[keylen + 1], val, vallen);
  keyval[keylen + vallen + 1] = '\0';

  int rval = bfam_critbit0_insert(&(d->t), keyval);

  if (rval == 2)
    ++d->num_entries;

  bfam_free(keyval);
  return rval;
}

int bfam_dictionary_insert_ptr(bfam_dictionary_t *d, const char *key,
                               const void *val_ptr)
{
  char val_str[BFAM_PTR_STR_LEN + 1];
  snprintf(val_str, BFAM_PTR_STR_LEN + 1, "%p", val_ptr);
  return bfam_dictionary_insert(d, key, val_str);
}

int bfam_dictionary_insert_int(bfam_dictionary_t *d, const char *key,
                               const int val)
{
  char val_str[BFAM_BUFSIZ];
  snprintf(val_str, BFAM_BUFSIZ, "%d", val);
  return bfam_dictionary_insert(d, key, val_str);
}

int bfam_dictionary_insert_locidx(bfam_dictionary_t *d, const char *key,
                                  const bfam_locidx_t val)
{
  char val_str[BFAM_BUFSIZ];
  snprintf(val_str, BFAM_BUFSIZ, "%" BFAM_LOCIDX_PRId, val);
  return bfam_dictionary_insert(d, key, val_str);
}

static int bfam_dictionary_get_value_handle(const char *keyval, void *arg)
{
  char *key = (char *)((void **)arg)[0];
  char **val = (char **)((void **)arg)[1];
  const size_t keylen = strlen(key);

  *val = (char *)&keyval[keylen];

  return 1;
}

/** Return a value given a key
 *
 * It takes a dictionary, \a d, returns a pointer to the value associated with
 * a \c NULL terminated \a key.
 *
 * \param [in] d dictionary
 * \param [in] key possible key
 * \returns:
 *   $\cases{ \c NULL & if {\it key} {\rm is not a member} \cr
 *          {\rm pointer to value} & if {\it key} {\rm is a member}$
 */
static char *bfam_dictionary_get_value(bfam_dictionary_t *d, const char *key)
{
  if (!bfam_dictionary_contains(d, key))
    return NULL;

  const size_t keylen = strlen(key);

  char *u = bfam_malloc(sizeof(char) * (keylen + 2));

  memcpy(u, key, keylen);
  u[keylen] = BFAM_KEYVALUE_SPLIT;
  u[keylen + 1] = '\0';

  char *value = NULL;
  void *arg[2];
  arg[0] = u;
  arg[1] = &value;

  bfam_critbit0_allprefixed(&(d->t), u, &bfam_dictionary_get_value_handle, arg);

  bfam_free(u);

  return value;
}

/** Delete key and value pair into a dictionary
 *
 * It takes a dictionary, \a d, and possibly mutates it such that a \c NULL
 * terminated string, \a key  is not member on exit.
 *
 * \param [in,out] d dictionary
 * \param [in] key possible key
 * \returns:
 *   $\cases{ 0 &if {\rm key not found} \cr
 *            1 &if {\it key deleted} }$.
 */
static int bfam_dictionary_delete(bfam_dictionary_t *d, const char *key)
{
  const size_t keylen = strlen(key);
  char *val = bfam_dictionary_get_value(d, key);
  if (val == NULL)
    return 0;
  const char *keyval = val - keylen - 1;
  return bfam_critbit0_delete(&(d->t), keyval);
}

void *bfam_dictionary_get_value_ptr(bfam_dictionary_t *d, const char *key)
{
  char *val_str = bfam_dictionary_get_value(d, key);
  if (val_str == NULL)
    return NULL;
  void *val_ptr = NULL;
  sscanf(val_str, "%p", &val_ptr);
  return val_ptr;
}

int bfam_dictionary_get_value_locidx(bfam_dictionary_t *d, const char *key,
                                     bfam_locidx_t *val)
{
  char *val_str = bfam_dictionary_get_value(d, key);
  if (val_str == NULL)
    return 0;
  int n = sscanf(val_str, "%" BFAM_LOCIDX_SCNd, val);
  return n;
}

typedef struct
{
  int (*handle)(const char *, const char *, void *);
  void *arg;
} bfam_dict_allprex;

static int bfam_dictionary_allprefixed_usercall(const char *keyval, void *arg)
{
  char *key = (char *)keyval;
  char *split = strchr(key, BFAM_KEYVALUE_SPLIT);
  *split = '\0';
  bfam_dict_allprex *s_arg = (bfam_dict_allprex *)arg;
  s_arg->handle(key, split + 1, s_arg->arg);
  *split = BFAM_KEYVALUE_SPLIT;
  return 1;
}

int bfam_dictionary_allprefixed(bfam_dictionary_t *d, const char *prefix,
                                int (*handle)(const char *, const char *,
                                              void *),
                                void *arg)
{
  bfam_dict_allprex args = {0, 0};
  args.handle = handle;
  args.arg = arg;
  return bfam_critbit0_allprefixed(
      &(d->t), prefix, &bfam_dictionary_allprefixed_usercall, &args);
}

typedef struct
{
  int (*handle)(const char *, void *, void *);
  void *arg;
} bfam_dict_allprex_ptr;

static int bfam_dictionary_allprefixed_usercall_ptr(const char *keyval,
                                                    void *arg)
{
  char *key = (char *)keyval;
  char *split = strchr(key, BFAM_KEYVALUE_SPLIT);
  *split = '\0';

  void *val_ptr = NULL;
  sscanf(split + 1, "%p", &val_ptr);

  bfam_dict_allprex_ptr *s_arg = (bfam_dict_allprex_ptr *)arg;
  s_arg->handle(key, val_ptr, s_arg->arg);

  *split = BFAM_KEYVALUE_SPLIT;
  return 1;
}

int bfam_dictionary_allprefixed_ptr(bfam_dictionary_t *d, const char *prefix,
                                    int (*handle)(const char *, void *, void *),
                                    void *arg)
{
  bfam_dict_allprex_ptr args = {0, 0};
  args.handle = handle;
  args.arg = arg;
  return bfam_critbit0_allprefixed(
      &(d->t), prefix, &bfam_dictionary_allprefixed_usercall_ptr, &args);
}

// }}}

// {{{ logging

static const char *const bfam_log_level_names[] = {
    "ALWAYS", "TRACE", "DEBUG", " VERB", " INFO", " WARN", "ERROR", "SILENT",
};

static const int bfam_log_root_rank = 0;
static int bfam_log_rank = 0;
static FILE *bfam_log_stream = NULL;
static int bfam_log_threshold = BFAM_LL_ALWAYS;

void bfam_log_init(int rank, FILE *stream, int threshold)
{

  bfam_log_rank = rank;
  if (stream == NULL)
  {
    bfam_log_stream = stdout;
  }
  else
  {
    bfam_log_stream = stream;
  }

  if (threshold == BFAM_LL_DEFAULT)
  {
    bfam_log_threshold = BFAM_LL_INFO;
  }
  else
  {
    BFAM_ABORT_IF_NOT(threshold <= BFAM_LL_SILENT &&
                          threshold >= BFAM_LL_ALWAYS,
                      "Invalid logging threshold");
    bfam_log_threshold = threshold;
  }
}

void bfam_log_printf(const char *file, int line, int category, int level, ...)
{
  FILE *stream = (bfam_log_stream) ? bfam_log_stream : stdout;
  va_list ap;
  const char *fmt;

  BFAM_ASSERT(level <= BFAM_LL_SILENT && level >= BFAM_LL_ALWAYS);

  if (bfam_log_threshold > level)
    return;

  if (category == BFAM_LC_ROOT && bfam_log_rank != bfam_log_root_rank)
    return;

  if (category == BFAM_LC_ALL)
  {
    fprintf(stream, "[%3d] ", bfam_log_rank);
  }
  else
  {
    fprintf(stream, "[   ] ");
  }
  fprintf(stream, "%s: ", bfam_log_level_names[level]);

  if (level == BFAM_LL_TRACE)
  {

    fprintf(stream, "%s:%3d ", file, line);
  }

  va_start(ap, level);
  fmt = va_arg(ap, const char *);
  vfprintf(stream, fmt, ap);
  fprintf(stream, "\n");
  va_end(ap);
}

// }}}

// {{{ utilities

static void bfam_util_strcsl(char *str, const char **list)
{
  str[0] = '\0';

  if (list != NULL)
  {
    for (size_t l = 0; list[l]; ++l)
    {
      strcat(str, list[l]);
      if (list[l + 1])
        strcat(str, ",");
    }
  }
}

static void bfam_util_mtranspose(size_t m, size_t n,
                                 bfam_long_real_t *restrict A, size_t lda,
                                 bfam_long_real_t *restrict B, size_t ldb)
{
  for (size_t i = 0; i < m; ++i)
    for (size_t j = 0; j < n; ++j)
      B[i * ldb + j] = A[j * lda + i];
}

static void bfam_util_mTmmult(size_t m, size_t n, size_t k,
                              bfam_long_real_t *restrict A, size_t lda,
                              bfam_long_real_t *restrict B, size_t ldb,
                              bfam_long_real_t *restrict C, size_t ldc)
{
  for (size_t i = 0; i < m; ++i)
    for (size_t j = 0; j < n; ++j)
      for (size_t p = 0; p < k; ++p)
        C[j * ldc + i] += A[i * lda + p] * B[j * ldb + p];
}

static void bfam_util_mmmult(size_t m, size_t n, size_t k,
                             bfam_long_real_t *restrict A, size_t lda,
                             bfam_long_real_t *restrict B, size_t ldb,
                             bfam_long_real_t *restrict C, size_t ldc)
{
  for (size_t i = 0; i < m; ++i)
    for (size_t j = 0; j < n; ++j)
      for (size_t p = 0; p < k; ++p)
        C[j * ldc + i] += A[p * lda + i] * B[j * ldb + p];
}

static void bfam_util_lu_factor(size_t n, bfam_long_real_t *restrict A,
                                size_t *restrict p, size_t *restrict q)
{
  /*
   * Algorithm 3.4.2 (Gaussian Elimination with Complete Pivoting) from p. 118
   * of Matrix Computations, Third Edition, by Golub and van Loan.
   */

  /*
   * Initialize pivots
   */
  for (size_t k = 0; k < n; ++k)
  {
    p[k] = q[k] = k;
  }

  for (size_t k = 0; k < n - 1; ++k)
  {
    size_t mu = k, lambda = k;
    bfam_long_real_t a_max = 0.0;

    for (size_t j = k; j < n; ++j)
    {
      for (size_t i = k; i < n; ++i)
      {
        const bfam_long_real_t a_abs = BFAM_LONG_REAL_ABS(A[i + n * j]);
        if (a_abs > a_max)
        {
          a_max = a_abs;
          mu = i;
          lambda = j;
        }
      }
    }

    /*
     * Swap rows
     */
    const size_t ptmp = p[k];
    p[k] = p[mu];
    p[mu] = ptmp;
    for (size_t j = 0; j < n; ++j)
    {
      const bfam_long_real_t rtmp = A[k + n * j];
      A[k + n * j] = A[mu + n * j];
      A[mu + n * j] = rtmp;
    }

    /*
     * Swap columns
     */
    const size_t qtmp = q[k];
    q[k] = q[lambda];
    q[lambda] = qtmp;
    for (size_t i = 0; i < n; ++i)
    {
      const bfam_long_real_t rtmp = A[i + n * k];
      A[i + n * k] = A[i + n * lambda];
      A[i + n * lambda] = rtmp;
    }

    if (BFAM_LONG_REAL_ABS(A[k + n * k]) >
        (BFAM_LONG_REAL_EPS * BFAM_LONG_REAL_EPS))
    {
      for (size_t i = k + 1; i < n; ++i)
      {
        A[i + n * k] = A[i + n * k] / A[k + n * k];
      }

      for (size_t i = k + 1; i < n; ++i)
      {
        for (size_t j = k + 1; j < n; ++j)
        {
          A[i + n * j] = A[i + n * j] - A[i + n * k] * A[k + n * j];
        }
      }
    }
  }
}

static void bfam_util_lu_solve(size_t n, bfam_long_real_t *restrict LU,
                               bfam_long_real_t *restrict x, size_t *restrict p,
                               size_t *restrict q,
                               bfam_long_real_t *restrict work)
{
  /*
   * Compute $Pb$
   */
  for (size_t k = 0; k < n; ++k)
  {
    work[k] = x[p[k]];
  }

  /*
   * Algorithm 3.1.3 (Forward Substitution: Column Version) from p. 90 of
   * Matrix Computations, Third Edition, by Golub and van Loan.
   *
   * Note: here we have L is unit lower triangular.
   *
   * Solve $Ly=Pb$.
   */
  for (size_t j = 0; j < n - 1; ++j)
  {
    /*
     * work[j] = work[j] / LU[j + n * j];
     */
    for (size_t i = j + 1; i < n; ++i)
    {
      work[i] = work[i] - work[j] * LU[i + n * j];
    }
  }
  /*
   * work[n - 1] = work[n - 1] / LU[n - 1 + n * (n - 1)];
   */

  /*
   * Algorithm 3.1.4 (Back Substitution: Column Version) from p. 90 of
   * Matrix Computations, Third Edition, by Golub and van Loan.
   *
   * Solve $Uw=y$.
   */
  for (size_t j = n - 1; j > 0; --j)
  {
    work[j] = work[j] / LU[j + n * j];
    for (size_t i = 0; i < j; ++i)
    {
      work[i] = work[i] - work[j] * LU[i + n * j];
    }
  }
  work[0] = work[0] / LU[0 + n * 0];

  /*
   * Compute $Qw$
   */
  for (size_t k = 0; k < n; ++k)
  {
    x[q[k]] = work[k];
  }
}

static void bfam_util_backslash(size_t m, size_t n,
                                bfam_long_real_t *restrict A,
                                bfam_long_real_t *restrict B,
                                bfam_long_real_t *restrict C)
{
  bfam_long_real_t *LU = bfam_malloc_aligned(m * m * sizeof(bfam_long_real_t));
  bfam_long_real_t *work = bfam_malloc_aligned(m * sizeof(bfam_long_real_t));

  size_t *p = bfam_malloc_aligned(m * sizeof(size_t));
  size_t *q = bfam_malloc_aligned(m * sizeof(size_t));

  memcpy(LU, A, m * m * sizeof(bfam_long_real_t));
  memcpy(C, B, m * n * sizeof(bfam_long_real_t));

  bfam_util_lu_factor(m, LU, p, q);

  for (size_t j = 0; j < n; ++j)
    bfam_util_lu_solve(m, LU, C + j * m, p, q, work);

  bfam_free_aligned(q);
  bfam_free_aligned(p);
  bfam_free_aligned(work);
  bfam_free_aligned(LU);
}

static void bfam_util_forwardslash(size_t m, size_t n,
                                   bfam_long_real_t *restrict A,
                                   bfam_long_real_t *restrict B,
                                   bfam_long_real_t *restrict C)
{
  bfam_long_real_t *AT = bfam_malloc_aligned(n * m * sizeof(bfam_long_real_t));
  bfam_long_real_t *BT = bfam_malloc_aligned(n * n * sizeof(bfam_long_real_t));
  bfam_long_real_t *CT = bfam_malloc_aligned(n * m * sizeof(bfam_long_real_t));

  bfam_util_mtranspose(m, n, A, m, AT, n);
  bfam_util_mtranspose(n, n, B, n, BT, n);

  bfam_util_backslash(n, m, BT, AT, CT);

  bfam_util_mtranspose(n, m, CT, n, C, m);

  bfam_free_aligned(CT);
  bfam_free_aligned(BT);
  bfam_free_aligned(AT);
}

/*
 * Integer power routine from:
 *   http://stackoverflow.com/questions/101439/the-most-efficient-way-to-implement-an-integer-based-power-function-powint-int
 */
static int bfam_ipow(int base, int exp)
{
  BFAM_ASSERT(exp >= 0);

  int result = 1;
  while (exp)
  {
    if (exp & 1)
      result *= base;
    exp >>= 1;
    base *= base;
  }

  return result;
}

/*
 * This uses a posix compilent solution from:
 *   https://www.securecoding.cert.org/confluence/display/c/FIO19-C.+Do+not+use+fseek%28%29+and+ftell%28%29+to+compute+the+size+of+a+regular+file
 */
static size_t bfam_util_file_size(const char *filename)
{
  struct stat stbuf;
  int fd, err;

  fd = open(filename, O_RDONLY);
  if (fd == -1)
  {
    BFAM_LERROR("Error opening file: %s", filename);
    BFAM_ABORT("Failed open");
  }

  if ((fstat(fd, &stbuf) != 0) || (!S_ISREG(stbuf.st_mode)))
  {
    BFAM_LERROR("Error determining size of the file: %s", filename);
    BFAM_ABORT("Failed fstat");
  }

  BFAM_ASSERT(stbuf.st_size >= 0);

  err = close(fd);
  if (err)
  {
    BFAM_LERROR("Error closing file: %s", filename);
    BFAM_ABORT("Failed close");
  }

  return (size_t)stbuf.st_size;
}

char *bfam_util_read_file(const char *filename, size_t *len)
{
  size_t readsize;
  char *buffer;
  size_t filesize = bfam_util_file_size(filename);
  FILE *stream = fopen(filename, "r");
  if (!stream)
  {
    BFAM_LERROR("Error opening the file: %s", filename);
    BFAM_ABORT("Failed fopen");
  }

  buffer = bfam_malloc(filesize + 1);
  readsize = fread(buffer, sizeof(char), filesize, stream);
  buffer[filesize] = '\0';

  if (readsize != filesize)
  {
    BFAM_LERROR("Error determining reading the file: %s", filename);
    BFAM_ABORT("Failed fread");
  }

  if (fclose(stream))
  {
    BFAM_LERROR("Error closing the file: %s", filename);
    BFAM_ABORT("Failed fclose");
  }

  if (len)
    *len = filesize;

  return buffer;
}
/** Return a value given a key assuming value is a \c int.
 *
 * It takes a dictionary, \a d, returns an \c int associated with a \c NULL
 * terminated \a key.
 *
 * \param [in]  d dictionary
 * \param [in]  key possible key
 * \param [out] val value if function returned \c 1
 *
 * \returns:
 *   $\cases{ \c 0 & if {\it key} {\rm is not a member} \cr
 *            \c 1 & if {\it key} {\rm is a member}$
 */
static int bfam_dictionary_get_value_int(bfam_dictionary_t *d, const char *key,
                                         int *val)
{
  char *val_str = bfam_dictionary_get_value(d, key);
  if (val_str == NULL)
    return 0;
  int n = sscanf(val_str, "%d", val);
  return n;
}

int bfam_util_get_host_rank(MPI_Comm comm)
{
  int host_rank, rank, size, length;
  char name[MPI_MAX_PROCESSOR_NAME];

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Get_processor_name(name, &length);

  int *host_ranks = bfam_malloc(size * sizeof(int));
  char *allnames = bfam_malloc(size * MPI_MAX_PROCESSOR_NAME * sizeof(char));

  BFAM_MPI_CHECK(MPI_Allgather(name, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, allnames,
                               MPI_MAX_PROCESSOR_NAME, MPI_CHAR, comm));

  bfam_dictionary_t hosts;
  bfam_dictionary_init(&hosts);

  for (int h = 0; h < size; ++h)
  {
    const char *key = allnames + h * MPI_MAX_PROCESSOR_NAME;
    int hr = 0;
    if (bfam_dictionary_get_value_int(&hosts, key, &hr))
    {
      ++hr;
      if (1 != bfam_dictionary_delete(&hosts, key))
        BFAM_ABORT("get host rank dictionary delete fail");
    }
    if (2 != bfam_dictionary_insert_int(&hosts, key, hr))
      BFAM_ABORT("get host rank dictionary insert fail");

    host_ranks[h] = hr;
  }

  BFAM_ROOT_VERBOSE("-----------Host Ranks----------");
  for (int h = 0; h < size; ++h)
  {
    BFAM_ROOT_VERBOSE("  host_rank[%3d] = %2d", h, host_ranks[h]);
  }
  BFAM_ROOT_VERBOSE("-------------------------------");

  host_rank = host_ranks[rank];

  bfam_dictionary_clear(&hosts);
  bfam_free(allnames);
  bfam_free(host_ranks);

  return host_rank;
}

// }}}

// {{{ gopt

struct opt_spec_s
{
  int key;
  int flags;
  const char *shorts;
  const char *const *longs;
};
typedef struct opt_spec_s opt_spec_t;

struct opt_s
{
  int key;
  const char *arg;
};
typedef struct opt_s opt_t;

void *bfam_gopt_sort(int *argc, const char **argv, const void *opt_specs)
{
  void *opts;
  {{{const char *const *arg_p = argv + 1;
  size_t opt_count = 1;
  for (; *arg_p; ++arg_p)
  {
    if ('-' == (*arg_p)[0] && (*arg_p)[1])
    {
      if ('-' == (*arg_p)[1])
      {
        if ((*arg_p)[2])
          ++opt_count;
        else
          break;
      }
      else
      {
        const opt_spec_t *opt_spec_p = opt_specs;
        for (; opt_spec_p->key; ++opt_spec_p)
        {
          if (strchr(opt_spec_p->shorts, (*arg_p)[1]))
          {
            opt_count +=
                opt_spec_p->flags & BFAM_GOPT_ARG ? 1 : strlen((*arg_p) + 1);
            break;
          }
        }
      }
    }
  }
  opts = bfam_malloc(opt_count * sizeof(opt_t));
}
}
}
{
  const char **arg_p = argv + 1;
  const char **next_operand = arg_p;
  opt_t *next_option = opts;

  if (!opts)
  {
    perror(argv[0]);
    exit(EX_OSERR);
  }
  for (; *arg_p; ++arg_p)
    if ('-' == (*arg_p)[0] && (*arg_p)[1])
      if ('-' == (*arg_p)[1])
        if ((*arg_p)[2])
        {
          {
            {
              const opt_spec_t *opt_spec_p = opt_specs;
              const char *const *longs = opt_spec_p->longs;
              next_option->key = 0;
              while (*longs)
              {
                const char *option_cp = (*arg_p) + 2;
                const char *name_cp = *longs;
                while (*option_cp && *option_cp == *name_cp)
                {
                  ++option_cp;
                  ++name_cp;
                }
                if ('=' == *option_cp || !*option_cp)
                {
                  if (*name_cp)
                  {
                    if (next_option->key)
                    {
                      BFAM_LERROR(
                          "%s: --%.*s: abbreviated option is ambiguous\n",
                          argv[0], (int)(option_cp - ((*arg_p) + 2)),
                          (*arg_p) + 2);
                      bfam_free(opts);
                      exit(EX_USAGE);
                    }
                    next_option->key = opt_spec_p->key;
                  }
                  else
                  {
                    next_option->key = opt_spec_p->key;
                    goto found_long;
                  }
                }
                if (!*++longs)
                {
                  ++opt_spec_p;
                  if (opt_spec_p->key)
                    longs = opt_spec_p->longs;
                }
              }
              if (!next_option->key)
              {
                BFAM_LERROR("%s: --%.*s: unknown option\n", argv[0],
                            (int)strcspn((*arg_p) + 2, "="), (*arg_p) + 2);
                bfam_free(opts);
                exit(EX_USAGE);
              }
              for (opt_spec_p = opt_specs; opt_spec_p->key != next_option->key;
                   ++opt_spec_p)
                ;
            found_long:

              if (!(opt_spec_p->flags & BFAM_GOPT_REPEAT))
              {
                const opt_t *opt_p = opts;
                for (; opt_p != next_option; ++opt_p)
                  if (opt_p->key == opt_spec_p->key)
                  {
                    BFAM_LERROR("%s: --%.*s: option may not be repeated (in "
                                "any long or short form)\n",
                                argv[0], (int)strcspn((*arg_p) + 2, "="),
                                (*arg_p) + 2);
                    bfam_free(opts);
                    exit(EX_USAGE);
                  }
              }
              if (opt_spec_p->flags & BFAM_GOPT_ARG)
              {
                next_option->arg = strchr((*arg_p) + 2, '=') + 1;
                if ((char *)0 + 1 == next_option->arg)
                {
                  ++arg_p;
                  if (!*arg_p || ('-' == (*arg_p)[0] && (*arg_p)[1]))
                  {
                    BFAM_LERROR(
                        "%s: --%s: option requires an option argument\n",
                        argv[0], (*(arg_p - 1)) + 2);
                    bfam_free(opts);
                    exit(EX_USAGE);
                  }
                  next_option->arg = *arg_p;
                }
              }
              else
              {
                if (strchr((*arg_p) + 2, '='))
                {
                  BFAM_LERROR(
                      "%s: --%.*s: option may not take an option argument\n",
                      argv[0], (int)strcspn((*arg_p) + 2, "="), (*arg_p) + 2);
                  bfam_free(opts);
                  exit(EX_USAGE);
                }
                next_option->arg = NULL;
              }
              ++next_option;
            }
          }
        }
        else
        {
          for (++arg_p; *arg_p; ++arg_p)
            *next_operand++ = *arg_p;
          break;
        }
      else
      {
        {
          {
            const char *short_opt = (*arg_p) + 1;
            for (; *short_opt; ++short_opt)
            {
              const opt_spec_t *opt_spec_p = opt_specs;

              for (; opt_spec_p->key; ++opt_spec_p)
                if (strchr(opt_spec_p->shorts, *short_opt))
                {
                  if (!(opt_spec_p->flags & BFAM_GOPT_REPEAT))
                  {
                    const opt_t *opt_p = opts;
                    for (; opt_p != next_option; ++opt_p)
                      if (opt_p->key == opt_spec_p->key)
                      {
                        BFAM_LERROR("%s: -%c: option may not be repeated (in "
                                    "any long or short form)\n",
                                    argv[0], *short_opt);
                        bfam_free(opts);
                        exit(EX_USAGE);
                      }
                  }
                  next_option->key = opt_spec_p->key;

                  if (opt_spec_p->flags & BFAM_GOPT_ARG)
                  {
                    if (short_opt[1])
                      next_option->arg = short_opt + 1;

                    else
                    {
                      ++arg_p;
                      if (!*arg_p || ('-' == (*arg_p)[0] && (*arg_p)[1]))
                      {
                        BFAM_LERROR(
                            "%s: -%c: option requires an option argument\n",
                            argv[0], *short_opt);
                        bfam_free(opts);
                        exit(EX_USAGE);
                      }
                      next_option->arg = *arg_p;
                    }
                    ++next_option;
                    goto break_2;
                  }
                  next_option->arg = NULL;
                  ++next_option;
                  goto continue_2;
                }
              BFAM_LERROR("%s: -%c: unknown option\n", argv[0], *short_opt);
              bfam_free(opts);
              exit(EX_USAGE);
            continue_2:
              BFAM_NOOP();
            }
          break_2:
            BFAM_NOOP();
          }
        }
      }
    else
      *next_operand++ = *arg_p;

  next_option->key = 0;
  *next_operand = NULL;
  *argc = (int)(next_operand - argv);
}
return opts;
}

size_t bfam_gopt(const void *vptr_opts, int key)
{
  const opt_t *opts = vptr_opts;
  size_t count = 0;
  for (; opts->key; ++opts)
    count += opts->key == key;

  return count;
}

void bfam_gopt_free(void *vptr_opts) { bfam_free(vptr_opts); }

// }}}

#ifdef BFAM_USE_LUA
// {{{ lua

#define BFAM_LUA_MAX_COMMAND_LEN 4096
#define BFAM_LUA_EVALEXP_VAR "XXX_bfam_evalexpr_XXX"

int bfam_lua_global_function_call(lua_State *L, const char *name,
                                  const char *sig, ...)
{
  va_list vl;
  int num_arg = 0;
  int num_res = 0;

  va_start(vl, sig);

  char buf[BFAM_LUA_MAX_COMMAND_LEN];
  /* Assign the Lua expression to a Lua global variable. */
  snprintf(buf, BFAM_LUA_MAX_COMMAND_LEN, BFAM_LUA_EVALEXP_VAR "=%s", name);
  if (!luaL_dostring(L, buf))
  {
    /* Get the value of the global varibable */
    lua_getglobal(L, BFAM_LUA_EVALEXP_VAR);
  }
  else
  {
    BFAM_ROOT_WARNING("function `%s' not found in lua file", name);
    return 1;
  }

  if (!lua_isfunction(L, -1))
  {
    BFAM_ROOT_WARNING("function `%s' not found in lua file", name);
    lua_pop(L, 1);
    return 1;
  }

  for (num_arg = 0; sig[num_arg] && sig[num_arg] != '>'; num_arg++)
  {
    luaL_checkstack(L, 1, "too many arguments");

    switch (sig[num_arg])
    {
    case 'd':
      lua_pushnumber(L, va_arg(vl, double));
      break;
    case 'l':
      lua_pushnumber(L, (double)va_arg(vl, bfam_long_real_t));
      break;
    case 'r':
      lua_pushnumber(L, (double)va_arg(vl, bfam_real_t));
      break;
    case 'i':
      lua_pushinteger(L, va_arg(vl, int));
      break;
    case 's':
      lua_pushstring(L, va_arg(vl, char *));
      break;
    case '>':
      break;
    default:
      BFAM_ABORT("function '%s' invalid input argument (%c)", name,
                 sig[num_arg]);
    }
  }

  BFAM_ABORT_IF_NOT(sig[num_arg] == '>', "arguments for '%s' does not contain "
                                         " a '>' character",
                    name);

  num_res = (int)strlen(sig) - num_arg - 1;

  BFAM_ABORT_IF_NOT(lua_pcall(L, num_arg, num_res, 0) == 0,
                    "error running function %s: %s", name, lua_tostring(L, -1));

  for (int n = 0; n < num_res; n++)
  {
    switch (sig[num_arg + 1 + n])
    {
    case 'd':
      BFAM_ABORT_IF_NOT(lua_isnumber(L, n - num_res),
                        "for '%s' return %d expected number got '%s'", name, n,
                        lua_tostring(L, n - num_res));
      *va_arg(vl, double *) = (double)lua_tonumber(L, n - num_res);
      break;
    case 'l':
      BFAM_ABORT_IF_NOT(lua_isnumber(L, n - num_res),
                        "for '%s' return %d expected number got '%s'", name, n,
                        lua_tostring(L, n - num_res));
      *va_arg(vl, bfam_long_real_t *) =
          (bfam_long_real_t)lua_tonumber(L, n - num_res);
      break;
    case 'r':
      BFAM_ABORT_IF_NOT(lua_isnumber(L, n - num_res),
                        "for '%s' return %d expected number got '%s'", name, n,
                        lua_tostring(L, n - num_res));
      *va_arg(vl, bfam_real_t *) = (bfam_real_t)lua_tonumber(L, n - num_res);
      break;
    case 'i':
      BFAM_ABORT_IF_NOT(lua_isnumber(L, n - num_res),
                        "for '%s' return %d expected number got '%s'", name, n,
                        lua_tostring(L, n - num_res));
      *va_arg(vl, int *) = (int)lua_tointeger(L, n - num_res);
      break;
    case 's':
      BFAM_ABORT_IF_NOT(lua_isstring(L, n - num_res),
                        "for '%s' return %d expected string got '%s'", name, n,
                        lua_tostring(L, n - num_res));
      *va_arg(vl, const char **) = lua_tostring(L, n - num_res);
      break;
    default:
      BFAM_ABORT("function '%s' invalid output argument (%c)", name,
                 sig[num_arg]);
    }
  }

  lua_pop(L, num_res);

  va_end(vl);
  return 0;
}

int bfam_lua_expr_boolean(lua_State *L, const char *expr, int def)
{
  int r = def;
  char buf[BFAM_LUA_MAX_COMMAND_LEN];
  snprintf(buf, BFAM_LUA_MAX_COMMAND_LEN, BFAM_LUA_EVALEXP_VAR "=%s", expr);
  if (!luaL_dostring(L, buf))
  {
    /* Get the value of the global varibable */
    lua_getglobal(L, BFAM_LUA_EVALEXP_VAR);
    if (lua_isboolean(L, -1))
      r = lua_toboolean(L, -1);
    /* remove lua_getglobal value */
    lua_pop(L, 1);
  }
  return r;
}

lua_Integer bfam_lua_expr_integer(lua_State *L, const char *expr,
                                  lua_Integer def)
{
  lua_Integer r = def;
  char buf[BFAM_LUA_MAX_COMMAND_LEN];
  /* Assign the Lua expression to a Lua global variable. */
  snprintf(buf, BFAM_LUA_MAX_COMMAND_LEN, BFAM_LUA_EVALEXP_VAR "=%s", expr);
  if (!luaL_dostring(L, buf))
  {
    /* Get the value of the global varibable */
    lua_getglobal(L, BFAM_LUA_EVALEXP_VAR);
    if (lua_isnumber(L, -1))
      r = lua_tointeger(L, -1);
    /* remove lua_getglobal value */
    lua_pop(L, 1);
  }
  return r;
}

lua_Number bfam_lua_expr_number(lua_State *L, const char *expr, lua_Number def)
{
  lua_Number r = def;
  char buf[BFAM_LUA_MAX_COMMAND_LEN];
  /* Assign the Lua expression to a Lua global variable. */
  snprintf(buf, BFAM_LUA_MAX_COMMAND_LEN, BFAM_LUA_EVALEXP_VAR "=%s", expr);
  if (!luaL_dostring(L, buf))
  {
    /* Get the value of the global varibable */
    lua_getglobal(L, BFAM_LUA_EVALEXP_VAR);
    if (lua_isnumber(L, -1))
      r = lua_tonumber(L, -1);
    /* remove lua_getglobal value */
    lua_pop(L, 1);
  }
  return r;
}

char *bfam_lua_expr_string(lua_State *L, const char *expr, const char *def)
{
  const char *r = def;
  size_t len = strlen(r);
  char buf[BFAM_LUA_MAX_COMMAND_LEN];
  /* Assign the Lua expression to a Lua global variable. */
  snprintf(buf, BFAM_LUA_MAX_COMMAND_LEN, BFAM_LUA_EVALEXP_VAR "=%s", expr);
  if (!luaL_dostring(L, buf))
  {
    /* Get the value of the global varibable */
    lua_getglobal(L, BFAM_LUA_EVALEXP_VAR);
    if (lua_isstring(L, -1))
      r = lua_tolstring(L, -1, &len);

    /* remove lua_getglobal value */
    lua_pop(L, 1);
  }

  char *s = bfam_malloc((len + 1) * sizeof(char));
  strncpy(s, r, len + 1);
  s[len] = '\0';

  return s;
}

// }}}
#endif

// {{{ subdomain

/*
 * Compare function which sorts a bfam_subdomain_face_map_entry_t array
 * in sending order.
 */

typedef struct bfam_subdomain_face_map_entry
{
  bfam_locidx_t np; /* Neighbor's processor number */
  bfam_locidx_t ns; /* Neighbor's subdomain id */
  bfam_locidx_t nk; /* Neighbor's element number */
  int8_t nf;        /* Neighbor's face number */
  int8_t nh;        /* Neighbor's hanging number */

  bfam_locidx_t s; /* Local subdomain id */
  bfam_locidx_t k; /* Local element number */
  int8_t f;        /* Local face number */
  int8_t h;        /* Local hanging number */
  int8_t o;        /* Local orientation */

  bfam_locidx_t id; /* Interface id */

  bfam_locidx_t gi; /* Index variable */
  bfam_locidx_t i;  /* Index variable */
} bfam_subdomain_face_map_entry_t;

static int bfam_subdomain_face_send_cmp(const void *a, const void *b)
{
  const bfam_subdomain_face_map_entry_t *la = a;
  const bfam_subdomain_face_map_entry_t *lb = b;

  if (la->np < lb->np)
    return -1;
  else if (la->np > lb->np)
    return 1;
  else if (la->id < lb->id)
    return -1;
  else if (la->id > lb->id)
    return 1;
  else if (la->ns < lb->ns)
    return -1;
  else if (la->ns > lb->ns)
    return 1;
  else if (la->s < lb->s)
    return -1;
  else if (la->s > lb->s)
    return 1;
  else if (la->nk < lb->nk)
    return -1;
  else if (la->nk > lb->nk)
    return 1;
  else if (la->nf < lb->nf)
    return -1;
  else if (la->nf > lb->nf)
    return 1;
  else if (la->nh < lb->nh)
    return -1;
  else if (la->nh > lb->nh)
    return 1;
  else
    return 0;

  BFAM_ABORT("We should never reach here.");
}

/*
 * Compare function which sorts a bfam_subdomain_face_map_entry_t array
 * in receiving order.
 */
static int bfam_subdomain_face_recv_cmp(const void *a, const void *b)
{
  const bfam_subdomain_face_map_entry_t *la = a;
  const bfam_subdomain_face_map_entry_t *lb = b;

  if (la->np < lb->np)
    return -1;
  else if (la->np > lb->np)
    return 1;
  else if (la->id < lb->id)
    return -1;
  else if (la->id > lb->id)
    return 1;
  else if (la->s < lb->s)
    return -1;
  else if (la->s > lb->s)
    return 1;
  else if (la->ns < lb->ns)
    return -1;
  else if (la->ns > lb->ns)
    return 1;
  else if (la->k < lb->k)
    return -1;
  else if (la->k > lb->k)
    return 1;
  else if (la->f < lb->f)
    return -1;
  else if (la->f > lb->f)
    return 1;
  else if (la->h < lb->h)
    return -1;
  else if (la->h > lb->h)
    return 1;
  else
    return 0;

  BFAM_ABORT("We should never reach here.");
}

/** initializes a subdomain glue data
 *
 * There is no new function for subdomains glue data since these are really just
 * a base class and a concrete grid and physics type should be defined
 *
 * \param [in,out] thisGlue pointer to the subdomain glue data
 * \param [in]     id of this side
 * \param [in]     sort id of this side
 * \param [in]     mpirank for this glue subdomain data
 * \param [in]     pointer to the minus side subdomain (can be \c NULL);
 */
static void bfam_subdomain_glue_init(bfam_subdomain_glue_data_t *glue,
                                     const bfam_locidx_t rank,
                                     const bfam_locidx_t id,
                                     const bfam_locidx_t id_s,
                                     bfam_subdomain_t *sub_m)
{
  glue->rank = rank;
  glue->sub_m = sub_m;
  glue->id = id;
  glue->id_s = id_s;
  bfam_dictionary_init(&glue->fields);
}

/** free up the memory allocated by the subdomain
 *
 * \param [in,out] thisSubdomain subdomain to clean up
 */
static void bfam_subdomain_free(bfam_subdomain_t *thisSubdomain)
{
  thisSubdomain->id = -1;
  thisSubdomain->uid = -1;
  bfam_free(thisSubdomain->name);
  bfam_critbit0_clear(&thisSubdomain->tags);
  bfam_dictionary_clear(&thisSubdomain->fields);
  bfam_dictionary_clear(&thisSubdomain->fields_face);

  thisSubdomain->tags.root = NULL;

  thisSubdomain->vtk_write_vtu_piece = NULL;
  thisSubdomain->field_add = NULL;

  thisSubdomain->glue_comm_info = NULL;

  if (thisSubdomain->glue_m)
  {
    bfam_dictionary_clear(&thisSubdomain->glue_m->fields);
    thisSubdomain->glue_m->rank = -1;
    thisSubdomain->glue_m->sub_m = NULL;
    thisSubdomain->glue_m->id = -1;
    thisSubdomain->glue_m->id_s = -1;
  }
  if (thisSubdomain->glue_p)
  {
    bfam_dictionary_clear(&thisSubdomain->glue_p->fields);
    thisSubdomain->glue_p->rank = -1;
    thisSubdomain->glue_p->sub_m = NULL;
    thisSubdomain->glue_p->id = -1;
    thisSubdomain->glue_p->id_s = -1;
  }
}

/** initializes a subdomain
 *
 * There is no new function for subdomains since these are really just a base
 * class and a concrete grid and physics type should be defined
 *
 * \param [in,out] thisSubdomain pointer to the subdomain
 * \param [in]     id   Unique id number for this subdomain
 * \param [in]     uid  user id number for this subdomain
 * \param [in]     name Name of this subdomain
 */
static void bfam_subdomain_init(bfam_subdomain_t *thisSubdomain,
                                bfam_locidx_t id, bfam_locidx_t uid,
                                const char *name)
{
  thisSubdomain->id = id;
  thisSubdomain->uid = uid;
  const size_t len = strlen(name);
  thisSubdomain->name = bfam_malloc((len + 1) * sizeof(char));
  strncpy(thisSubdomain->name, name, len + 1);

  bfam_dictionary_init(&thisSubdomain->fields);
  bfam_dictionary_init(&thisSubdomain->fields_face);

  thisSubdomain->tags.root = NULL;

  thisSubdomain->free = bfam_subdomain_free;

  thisSubdomain->vtk_write_vtu_piece = NULL;

  thisSubdomain->field_add = NULL;

  thisSubdomain->glue_comm_info = NULL;

  thisSubdomain->glue_m = NULL;
  thisSubdomain->glue_p = NULL;

  thisSubdomain->user_data = NULL;
}

void bfam_subdomain_add_tag(bfam_subdomain_t *thisSubdomain, const char *tag)
{
  BFAM_LDEBUG("subdomain %s: adding tag %s", thisSubdomain->name, tag);
  int r = bfam_critbit0_insert(&thisSubdomain->tags, tag);

  BFAM_ABORT_IF(!r, "Out of memory when adding tag: %s to subdomain %s", tag,
                thisSubdomain->name);
}

static int tag_prefix_match(const char *tag, void *args) { return 0; }

int bfam_subdomain_has_tag(bfam_subdomain_t *thisSubdomain, const char *tag)
{
  const size_t tlen = strlen(tag);
  if (tag[tlen - 1] == '*')
  {
    char newtag[BFAM_BUFSIZ];
    strncpy(newtag, tag, tlen - 1);
    newtag[tlen - 1] = '\0';
    return 0 == bfam_critbit0_allprefixed(&thisSubdomain->tags, newtag,
                                          tag_prefix_match, NULL);
  }
  else
    return bfam_critbit0_contains(&thisSubdomain->tags, tag);
}

/** Add a field to the subdomain
 *
 * \param [in,out] thisSubdomain subdomain to search for the tag
 * \param [in]     name          name of the field to add to the subdomain
 *
 * \returns:
 *   $\cases{ 0 &if {\rm out of memory} \cr
 *            1 &if {\it name} {\rm was already a field} \cr
 *            2 &if {\it name} {\rm was added successfully}}$.
 */
static int bfam_subdomain_field_add(bfam_subdomain_t *thisSubdomain,
                                    const char *name)
{
  if (thisSubdomain->field_add)
  {
    BFAM_VERBOSE("subdomain %s: adding field %s", thisSubdomain->name, name);
    return thisSubdomain->field_add(thisSubdomain, name);
  }
  else
  {
    BFAM_VERBOSE("subdomain %s cannot add field %s (no field_add)",
                 thisSubdomain->name, name);
    return 0;
  }
}

// }}}

// {{{ domain

#define BFAM_DEFAULT_SUBDOMAIN_ROOT 0

/* Domain based functions */

void bfam_domain_init(bfam_domain_t *thisDomain, MPI_Comm domComm)
{
  const bfam_locidx_t sizeSubdomains = 16;
  thisDomain->comm = domComm; // Perhaps we should duplicate it?
  thisDomain->num_subdomains = 0;
  thisDomain->sizeSubdomains = sizeSubdomains;
  thisDomain->subdomains =
      bfam_malloc(sizeSubdomains * sizeof(bfam_subdomain_t *));
  bfam_dictionary_init(&(thisDomain->name2num));
}

/** Clean up domain
 *
 * frees any mememory allocated by the domain and calls free command on all
 * subdomains
 *
 * \param [in,out] domain domain to clean up
 */
static void bfam_domain_free(bfam_domain_t *thisDomain)
{
  for (bfam_locidx_t i = 0; i < thisDomain->num_subdomains; i++)
  {
    thisDomain->subdomains[i]->free(thisDomain->subdomains[i]);
    bfam_free(thisDomain->subdomains[i]);
  }
  thisDomain->comm = MPI_COMM_NULL;
  thisDomain->num_subdomains = 0;
  thisDomain->sizeSubdomains = 0;
  bfam_free(thisDomain->subdomains);
  thisDomain->subdomains = NULL;
  bfam_dictionary_clear(&thisDomain->name2num);
}

static int bfam_domain_compare_subdomain_by_id(const void *a, const void *b)
{
  const bfam_subdomain_t *subA = *(bfam_subdomain_t * const *)a;
  const bfam_subdomain_t *subB = *(bfam_subdomain_t * const *)b;

  int rval;

  if (subA->id < subB->id)
    rval = -1;
  else if (subA->id > subB->id)
    rval = 1;
  else
    rval = 0;

  return rval;
}

void bfam_domain_get_subdomains(bfam_domain_t *thisDomain,
                                bfam_domain_match_t matchType,
                                const char **tags, bfam_locidx_t numEntries,
                                bfam_subdomain_t **subdomains,
                                bfam_locidx_t *num_subdomains)
{
  BFAM_ASSERT(subdomains != NULL);
  BFAM_ASSERT(num_subdomains != NULL);

  if (numEntries <= 0)
    return;

  *num_subdomains = 0;
  for (bfam_locidx_t d = 0; d < thisDomain->num_subdomains; ++d)
  {
    bfam_subdomain_t *subdomain = thisDomain->subdomains[d];
    int matched = 0;
    switch (matchType)
    {
    case BFAM_DOMAIN_OR:
      matched = 0;
      for (size_t t = 0; !matched && tags[t]; ++t)
      {
        int hasTag = bfam_subdomain_has_tag(subdomain, tags[t]);
        matched = hasTag || matched;
      }
      break;
    case BFAM_DOMAIN_AND:
      matched = 1;
      for (size_t t = 0; matched && tags[t]; ++t)
        matched = matched && bfam_subdomain_has_tag(subdomain, tags[t]);
      break;
    default:
      BFAM_ABORT("Unsupported Match Type");
    }

    if (matched)
    {
      subdomains[*num_subdomains] = subdomain;
      ++(*num_subdomains);
    }

    if (*num_subdomains == numEntries)
      break;
  }

  /*
   * Sort subdomains by id number
   */
  qsort(subdomains, *num_subdomains, sizeof(bfam_subdomain_t *),
        bfam_domain_compare_subdomain_by_id);

  return;
}

/** Add subdomain
 *
 * \param [in,out] thisDomain domain to add subdomain to
 * \param [in]     newSubdomain subdomain to add to the domain
 *
 * \return subdomain id
 */
static bfam_locidx_t bfam_domain_add_subdomain(bfam_domain_t *thisDomain,
                                               bfam_subdomain_t *newSubdomain)
{
  bfam_locidx_t sub_id;
  // double size
  if (thisDomain->num_subdomains == thisDomain->sizeSubdomains)
  {
    BFAM_VERBOSE("Doubling domain size");
    thisDomain->sizeSubdomains = 2 * thisDomain->sizeSubdomains;
    thisDomain->subdomains =
        bfam_realloc(thisDomain->subdomains,
                     thisDomain->sizeSubdomains * sizeof(bfam_subdomain_t *));
  }

  // create the key value pair
  BFAM_VERBOSE("adding subdomain %3" BFAM_LOCIDX_PRId " with name %s",
               thisDomain->num_subdomains, newSubdomain->name);

  // check if it's already there
  if (1 == bfam_dictionary_insert_locidx(&(thisDomain->name2num),
                                         newSubdomain->name,
                                         thisDomain->num_subdomains))
  {
    BFAM_ABORT("domain already contains subdomain \"%s\"", newSubdomain->name);
  }

  // add block
  sub_id = thisDomain->num_subdomains;
  thisDomain->subdomains[sub_id] = newSubdomain;
  thisDomain->num_subdomains++;
  return sub_id;
}

bfam_subdomain_t *bfam_domain_get_subdomain_by_num(bfam_domain_t *thisDomain,
                                                   bfam_locidx_t id)
{
  BFAM_ABORT_IF(id >= thisDomain->num_subdomains || id < 0,
                "Bad subdomain id: %jd", (intmax_t)id);
  return thisDomain->subdomains[id];
}

/** Add fields to subdomains matching the tags passed in.
 *
 * \param [in,out]  thisDomain    domain to search for subdomains in
 * \param [in]      matchType     type of match, \c BFAM_DOMAIN_OR will
 *                                match subdomains with any of the tags
 *                                and \c BFAM_DOMAIN_AND will match subdomains
 *                                with all of the tags.
 * \param [in]      mtags         \c NULL terminated array of the tags to match
 * \param [in]      tags          tags to add to the subdomains
 *
 */
static void bfam_domain_add_tags(bfam_domain_t *thisDomain,
                                 bfam_domain_match_t match, const char **mtags,
                                 const char **tags)
{
  bfam_subdomain_t **subdomains =
      bfam_malloc(thisDomain->num_subdomains * sizeof(bfam_subdomain_t **));

  bfam_locidx_t num_subdomains = 0;

  bfam_domain_get_subdomains(thisDomain, match, mtags,
                             thisDomain->num_subdomains, subdomains,
                             &num_subdomains);

  for (bfam_locidx_t s = 0; s < num_subdomains; ++s)
    for (size_t t = 0; tags[t]; ++t)
      bfam_subdomain_add_tag(subdomains[s], tags[t]);

  bfam_free(subdomains);
}

void bfam_domain_add_tag(bfam_domain_t *thisDomain, bfam_domain_match_t match,
                         const char **mtags, const char *tag)
{
  const char *tags[] = {tag, NULL};
  bfam_domain_add_tags(thisDomain, match, mtags, tags);
}

void bfam_domain_add_fields(bfam_domain_t *thisDomain,
                            bfam_domain_match_t match, const char **tags,
                            const char **fields)
{
  bfam_subdomain_t **subdomains =
      bfam_malloc(thisDomain->num_subdomains * sizeof(bfam_subdomain_t **));

  bfam_locidx_t num_subdomains = 0;

  bfam_domain_get_subdomains(thisDomain, match, tags,
                             thisDomain->num_subdomains, subdomains,
                             &num_subdomains);

  for (bfam_locidx_t s = 0; s < num_subdomains; ++s)
    for (size_t f = 0; fields[f]; ++f)
      bfam_subdomain_field_add(subdomains[s], fields[f]);

  bfam_free(subdomains);
}

// }}}

// {{{ communicator

#define BFAM_COMM_NUM_SRT (6)

typedef struct bfam_communicator_map_entry
{
  bfam_locidx_t np; /* Neighbor's processor number */
  bfam_locidx_t s[BFAM_COMM_NUM_SRT];
  /* sort param 0: typically Neighbor's subdomain id */
  /* sort param 1: typically My subdomain id */
  /* sort param 2: typically local ID (sbp) */
  /* sort param 3: typically fce ID (sbp) */
  bfam_locidx_t rank;          /*my rank*/
  bfam_subdomain_t *subdomain; /*pointer to subdomain*/

  bfam_locidx_t orig_order; /* original order */
} bfam_communicator_map_entry_t;

static int bfam_communicator_send_compare(const void *a, const void *b)
{
  const bfam_communicator_map_entry_t *mapA =
      (const bfam_communicator_map_entry_t *)a;
  const bfam_communicator_map_entry_t *mapB =
      (const bfam_communicator_map_entry_t *)b;

  if (mapA->np < mapB->np)
    return -1;
  if (mapA->np > mapB->np)
    return 1;
  for (int i = 0; i < BFAM_COMM_NUM_SRT; i++)
  {
    if (mapA->s[i] < mapB->s[i])
      return -1;
    if (mapA->s[i] > mapB->s[i])
      return 1;
  }
  BFAM_ABORT_IF(0 == 0, "Should not be same map!");

  return 0;
}

static int bfam_communicator_recv_compare(const void *a, const void *b)
{
  const bfam_communicator_map_entry_t *mapA =
      (const bfam_communicator_map_entry_t *)a;
  const bfam_communicator_map_entry_t *mapB =
      (const bfam_communicator_map_entry_t *)b;

  if (mapA->np < mapB->np)
    return -1;
  if (mapA->np > mapB->np)
    return 1;
  for (int i = 0; i < BFAM_COMM_NUM_SRT / 2; i++)
  {
    if (mapA->s[2 * i + 1] < mapB->s[2 * i + 1])
      return -1;
    if (mapA->s[2 * i + 1] > mapB->s[2 * i + 1])
      return 1;
    if (mapA->s[2 * i] < mapB->s[2 * i])
      return -1;
    if (mapA->s[2 * i] > mapB->s[2 * i])
      return 1;
  }
  BFAM_ABORT_IF(0 == 0, "Should not be same map!");

  return 0;
}

/** initializes a communicator
 *
 * \param [in,out] communicator pointer to the communicator
 * \param [in]     domain       domain to output to communicate
 * \param [in]     match        type of match, \c BFAM_DOMAIN_OR will
 *                              match subdomains with any of the tags
 *                              and \c BFAM_DOMAIN_AND will match subdomains
 *                              with all of the tags.
 * \param [in]     tags         \c NULL terminated array of the tags to match
 *                              glue grids doing the communication
 * \param [in]     comm         MPI communicator
 * \param [in]     tag          user specified communicator tag
 * \param [in]     userdata     user custom data to pass through
 */
static void bfam_communicator_init(bfam_communicator_t *communicator,
                                   bfam_domain_t *domain,
                                   bfam_domain_match_t match, const char **tags,
                                   MPI_Comm comm, int tag, void *user_args)
{
  BFAM_LDEBUG("Communicator Init");
  communicator->comm = comm;
  communicator->tag = tag;

  /* get the subdomains */
  bfam_subdomain_t *subdomains[domain->num_subdomains];
  bfam_domain_get_subdomains(domain, match, tags, domain->num_subdomains,
                             subdomains, &communicator->num_subs);

  communicator->sub_data =
      bfam_malloc(communicator->num_subs * sizeof(bfam_comm_subdata_t));

  size_t send_sz = 0;
  size_t recv_sz = 0;
  communicator->num_procs = 0;
  communicator->user_args = user_args;

  /* figure out the info for everyone */
  bfam_communicator_map_entry_t map[communicator->num_subs];
  bfam_critbit0_tree_t procs = {0};
  char procStr[BFAM_BUFSIZ];
  for (int s = 0; s < communicator->num_subs; s++)
  {
    map[s].subdomain = subdomains[s];

    BFAM_ABORT_IF(subdomains[s]->glue_comm_info == NULL,
                  "glue_comm_info not initialized for sudbomain %s",
                  subdomains[s]->name);

    for (int i = 0; i < BFAM_COMM_NUM_SRT; i++)
      map[s].s[i] = 0;
    subdomains[s]->glue_comm_info(
        subdomains[s], &map[s].np, map[s].s, BFAM_COMM_NUM_SRT,
        &communicator->sub_data[s].send_sz, &communicator->sub_data[s].recv_sz,
        communicator->user_args);

    BFAM_MPI_CHECK(MPI_Comm_rank(comm, &map[s].rank));

    send_sz += communicator->sub_data[s].send_sz;
    recv_sz += communicator->sub_data[s].recv_sz;

    communicator->sub_data[s].subdomain = subdomains[s];

    communicator->sub_data[s].send_buf = NULL;
    communicator->sub_data[s].recv_buf = NULL;

    map[s].orig_order = s;

    snprintf(procStr, BFAM_BUFSIZ, "%jd", (intmax_t)map[s].np);
    if (!bfam_critbit0_contains(&procs, procStr))
    {
      communicator->num_procs++;
      bfam_critbit0_insert(&procs, procStr);
    }
  }
  bfam_critbit0_clear(&procs);

  /* allocate everything now */
  communicator->send_sz = send_sz;
  communicator->send_buf = bfam_malloc_aligned(communicator->send_sz);

  communicator->recv_sz = recv_sz;
  communicator->recv_buf = bfam_malloc_aligned(communicator->recv_sz);

  communicator->proc_data =
      bfam_malloc(communicator->num_procs * sizeof(bfam_comm_procdata_t));

  communicator->send_request =
      bfam_malloc(2 * communicator->num_procs * sizeof(MPI_Request));
  communicator->recv_request =
      communicator->send_request + communicator->num_procs;

  communicator->send_status =
      bfam_malloc(2 * communicator->num_procs * sizeof(MPI_Status));
  communicator->recv_status =
      communicator->send_status + communicator->num_procs;

  for (int i = 0; communicator->num_procs > i; i++)
  {
    communicator->proc_data[i].send_sz = 0;
    communicator->proc_data[i].send_buf = NULL;

    communicator->proc_data[i].recv_sz = 0;
    communicator->proc_data[i].recv_buf = NULL;

    communicator->send_request[i] = MPI_REQUEST_NULL;
    communicator->recv_request[i] = MPI_REQUEST_NULL;
  }

  /* sort for send  and fill struct */
  qsort((void *)map, communicator->num_subs,
        sizeof(bfam_communicator_map_entry_t), bfam_communicator_send_compare);
  char *send_buf_ptr = communicator->send_buf;
  size_t send_offset = 0;
  bfam_locidx_t np = -1;   /* global proc ID */
  bfam_locidx_t proc = -1; /* local storage proc ID */
  for (int s = 0; s < communicator->num_subs; s++)
  {
    bfam_locidx_t t = map[s].orig_order;
    communicator->sub_data[t].send_buf = send_buf_ptr;
    communicator->sub_data[t].send_offset = send_offset;

    if (map[s].np != np)
    {
      np = map[s].np;
      proc++;
      communicator->proc_data[proc].send_buf = send_buf_ptr;
      communicator->proc_data[proc].rank = np;
    }
    BFAM_ABORT_IF_NOT(communicator->proc_data[proc].rank == np,
                      "problem with local proc ID");
    communicator->proc_data[proc].send_sz += communicator->sub_data[t].send_sz;

    send_buf_ptr += communicator->sub_data[t].send_sz;
    send_offset += communicator->sub_data[t].send_sz;
  }

  /* sort for recv and fill struct*/
  qsort((void *)map, communicator->num_subs,
        sizeof(bfam_communicator_map_entry_t), bfam_communicator_recv_compare);
  char *recv_buf_ptr = communicator->recv_buf;
  size_t recv_offset = 0;
  np = -1;   /* global proc ID */
  proc = -1; /* local storage proc ID */
  for (int s = 0; s < communicator->num_subs; s++)
  {
    bfam_locidx_t t = map[s].orig_order;
    communicator->sub_data[t].recv_buf = recv_buf_ptr;
    communicator->sub_data[t].recv_offset = recv_offset;
    if (map[s].np != np)
    {
      np = map[s].np;
      proc++;
      communicator->proc_data[proc].recv_buf = recv_buf_ptr;
    }
    BFAM_ABORT_IF_NOT(communicator->proc_data[proc].rank == np,
                      "problem with local proc ID");
    communicator->proc_data[proc].recv_sz += communicator->sub_data[t].recv_sz;

    recv_buf_ptr += communicator->sub_data[t].recv_sz;
    recv_offset += communicator->sub_data[t].recv_sz;
  }
}

bfam_communicator_t *bfam_communicator_new(bfam_domain_t *domain,
                                           bfam_domain_match_t match,
                                           const char **tags, MPI_Comm comm,
                                           int tag, void *user_args)
{
  bfam_communicator_t *newCommunicator =
      bfam_malloc(sizeof(bfam_communicator_t));

  bfam_communicator_init(newCommunicator, domain, match, tags, comm, tag,
                         user_args);
  return newCommunicator;
}

void bfam_communicator_free(bfam_communicator_t *communicator)
{
  BFAM_LDEBUG("Communicator Free");

  /* Just make sure there are no pending requests */
  BFAM_MPI_CHECK(MPI_Waitall(2 * communicator->num_procs,
                             communicator->send_request,
                             communicator->send_status));

  bfam_free_aligned(communicator->send_buf);
  bfam_free_aligned(communicator->recv_buf);
  bfam_free(communicator->sub_data);
  bfam_free(communicator->proc_data);
  bfam_free(communicator->send_request);
  bfam_free(communicator->send_status);
}

// }}}

// {{{ jacobi

/*
 * This function computes the normalization of the Jacobi polynomial
 * $\left\{h_N^{(\alpha,\beta)}\right\}^{-\frac12}$ where [see @Szego39
 * (4.3.3)]
 *
 * $$
 *  \begin{aligned}
 *    h_N^{(\alpha,\beta)} &=
 *      \int_{-1}^{+1} (1-x)^{\alpha} (1+x)^{\beta}
 *    \left\{P_N^{(\alpha,\beta)} (x)\right\}^2 \, dx \\
 *    &=
 *    \frac{2^{\alpha+\beta+1}}{2N+\alpha+\beta+1}
 *    \frac{\Gamma(N+\alpha+1)\Gamma(N+\beta+1)}
 *         {\Gamma(N+\alpha+\beta+1)\Gamma(N+1)}.
 *  \end{aligned}
 * $$
 */
static bfam_long_real_t bfam_jacobi_h_inv_sqrt(bfam_long_real_t alpha,
                                               bfam_long_real_t beta, int N)
{
  BFAM_ASSERT(N >= 0);
  BFAM_ASSERT(alpha >= BFAM_LONG_REAL(-1.0));
  BFAM_ASSERT(beta >= BFAM_LONG_REAL(-1.0));
  BFAM_ASSERT(!(BFAM_LONG_REAL_APPROX_EQ(alpha, BFAM_LONG_REAL(-0.5), 10) &&
                BFAM_LONG_REAL_APPROX_EQ(beta, BFAM_LONG_REAL(-0.5), 10)));

  bfam_long_real_t lgn = -(alpha + beta + 1) * BFAM_LONG_REAL_LOG(2) -
                         BFAM_LONG_REAL_LGAMMA(N + alpha + 1) -
                         BFAM_LONG_REAL_LGAMMA(N + beta + 1) +
                         BFAM_LONG_REAL_LOG(2 * N + alpha + beta + 1) +
                         BFAM_LONG_REAL_LGAMMA(N + 1) +
                         BFAM_LONG_REAL_LGAMMA(N + alpha + beta + 1);
  return BFAM_LONG_REAL_SQRT(BFAM_LONG_REAL_EXP(lgn));
}

/*
 * This function evaluates the orthonormal polynomial $p_N(x)$ associated with
 * the Jacobi polynomial where $p_N(x) =
 * \left\{h_N^{(\alpha,\beta)}\right\}^{-\frac12} P_N^{(\alpha,\beta)} (x)$.
 *
 * The Jacobi polynomials are a set of polynomial functions for $\alpha > -1$,
 * $\beta > -1$, and $N$ a non-negative integer.  The functions are defined
 * on $[-1, +1]$ and orthogonal with respect to the weight function
 * $$w(x)=(1-x)^\alpha(1+x)^\beta.$$  Here we use the same normalization as
 * @Szego39, i.e., $P_N^{(\alpha,\beta)}(1) = {n+\alpha\choose n}$.  Thus we
 * have
 * $$
 *   \int_{-1}^{+1} p_N(x) p_M(x) w(x) \, dx = \delta_{NM}.
 * $$
 *
 * The three term recurrence relation arrived at by rearranging @Szego39
 * [(4.5.1)] is
 * $$
 *   P_n^{(\alpha,\beta)}(x) = (ax-b) P_{n-1}^{(\alpha,\beta)}(x)
 *                             - c P_{n-2}^{(\alpha,\beta)}
 * $$
 * where
 * $$
 * \begin{aligned}
 *   a &= \frac{(2n + \alpha + \beta -1)(2n + \alpha + \beta)}
 *             {2n (n + \alpha + \beta)} \\
 *   b &= \frac{(\beta^2 - \alpha^2)(2n + \alpha + \beta - 1)}
 *             {2n(n + \alpha + \beta)(2n + \alpha + \beta - 2)} \\
 *   c &= \frac{(n + \alpha - 1)(n + \beta - 1)(2n + \alpha + \beta)}
 *             {n(n + \alpha + \beta)(2n + \alpha + \beta - 2)}
 * \end{aligned}
 * $$
 * with $P_0^{(\alpha,\beta)}(x) = 1$ and
 * $P_1^{(\alpha,\beta)}(x) =  \frac12(\alpha + \beta + 2)x
 *                           + \frac12(\alpha - \beta)$.
 */
static void bfam_jacobi_p(bfam_long_real_t alpha, bfam_long_real_t beta, int N,
                          size_t nx, bfam_long_real_t *x, bfam_long_real_t *P)
{
  BFAM_ASSERT(N >= 0);
  BFAM_ASSERT(alpha >= BFAM_LONG_REAL(-1.0));
  BFAM_ASSERT(beta >= BFAM_LONG_REAL(-1.0));
  BFAM_ASSERT(!(BFAM_LONG_REAL_APPROX_EQ(alpha, BFAM_LONG_REAL(-0.5), 10) &&
                BFAM_LONG_REAL_APPROX_EQ(beta, BFAM_LONG_REAL(-0.5), 10)));

  for (size_t i = 0; i < nx; ++i)
  {
    bfam_long_real_t P_n_2;
    bfam_long_real_t P_n_1 = 1;
    bfam_long_real_t P_n_0 =
        ((alpha + beta + 2) / 2) * x[i] + (alpha - beta) / 2;
    if (N == 0)
    {
      P[i] = P_n_1;
    }
    else if (N == 1)
    {
      P[i] = P_n_0;
    }
    else
    {
      for (int n = 2; n < N + 1; ++n)
      {
        bfam_long_real_t a = (2 * n + alpha + beta - 1) *
                             (2 * n + alpha + beta) /
                             (2 * n * (n + alpha + beta));
        bfam_long_real_t b =
            (beta * beta - alpha * alpha) * (2 * n + alpha + beta - 1) /
            (2 * n * (n + alpha + beta) * (2 * n + alpha + beta - 2));
        bfam_long_real_t c =
            (n + alpha - 1) * (n + beta - 1) * (2 * n + alpha + beta) /
            (n * (n + alpha + beta) * (2 * n + alpha + beta - 2));

        P_n_2 = P_n_1;
        P_n_1 = P_n_0;
        P_n_0 = (a * x[i] - b) * P_n_1 - c * P_n_2;
      }
      P[i] = P_n_0;
    }
  }

  /*
   * Normalize the Jacobi polynomials
   */
  bfam_long_real_t h_inv_sqrt = bfam_jacobi_h_inv_sqrt(alpha, beta, N);
  for (size_t i = 0; i < nx; ++i)
    P[i] *= h_inv_sqrt;

  return;
}

/*
 * This function evaluates the derivative of the orthonormal polynomial
 * $p_N(x)$ associated with the Jacobi polynomial where
 * $p_N(x) = \left\{h_N^{(\alpha,\beta)}\right\}^{-\frac12}
 * P_N^{(\alpha,\beta)} (x)$.
 *
 * For the evaluation of the derivative we use the identity
 * $$
 *   \frac{d}{dx} P_N^{(\alpha,\beta)} (x) =
 *    \frac{N+\alpha+\beta+1}{2} P_{N-1}^{(\alpha+1,\beta+1)} (x)
 * $$
 * along with
 * $$
 *   h_N^{(\alpha,\beta)} =
 *     \frac{N+\alpha+\beta+1}{4N} h_{N-1}^{(\alpha+1,\beta+1)}
 * $$
 * to get
 * $$
 * \begin{aligned}
 *   \frac{d}{dx} p_N^{(\alpha,\beta)} (x)
 *     &= \left\{h_N^{(\alpha,\beta)}\right\}^{-\frac12}
 *        \frac{d}{dx} P_N^{(\alpha,\beta)} (x) \\
 *     &= \left\{h_N^{(\alpha,\beta)}\right\}^{-\frac12}
 *        \frac{N+\alpha+\beta+1}{2}
 *        P_{N-1}^{(\alpha+1,\beta+1)} (x) \\
 *     &= \left(\frac{4N}{N+\alpha+\beta+1}\right)^{\frac12}
 *        \frac{N+\alpha+\beta+1}{2}
 *        \left\{h_{N-1}^{(\alpha+1,\beta+1)}\right\}^{-\frac12}
 *        P_{N-1}^{(\alpha+1,\beta+1)} (x) \\
 *     &= \left(N(N+\alpha+\beta+1)\right)^{\frac12}
 *        p_{N-1}^{(\alpha+1,\beta+1)} (x).
 * \end{aligned}
 * $$
 */
static void bfam_grad_jacobi_p(bfam_long_real_t alpha, bfam_long_real_t beta,
                               int N, size_t nx, bfam_long_real_t *x,
                               bfam_long_real_t *dP)
{
  BFAM_ASSERT(N >= 0);
  BFAM_ASSERT(alpha >= BFAM_LONG_REAL(-1.0));
  BFAM_ASSERT(beta >= BFAM_LONG_REAL(-1.0));

  if (N == 0)
  {
    for (size_t i = 0; i < nx; ++i)
    {
      dP[i] = BFAM_LONG_REAL(0.0);
    }
  }
  else
  {
    bfam_jacobi_p(alpha + 1, beta + 1, N - 1, nx, x, dP);
    bfam_long_real_t scale = BFAM_LONG_REAL_SQRT(N * (N + alpha + beta + 1));
    for (size_t i = 0; i < nx; ++i)
    {
      dP[i] *= scale;
    }
  }

  return;
}

static void bfam_jacobi_gauss_quadrature_half(bfam_long_real_t alpha,
                                              bfam_long_real_t beta, int N,
                                              int half,
                                              bfam_long_real_t *restrict x,
                                              bfam_long_real_t *restrict w)
{
  BFAM_ASSERT(N >= 0);
  BFAM_ASSERT(alpha >= BFAM_LONG_REAL(-1.0));
  BFAM_ASSERT(beta >= BFAM_LONG_REAL(-1.0));
  BFAM_ASSERT(!(BFAM_LONG_REAL_APPROX_EQ(alpha, BFAM_LONG_REAL(-0.5), 10) &&
                BFAM_LONG_REAL_APPROX_EQ(beta, BFAM_LONG_REAL(-0.5), 10)));

  const int MAX_ITERATIONS = 200;

  int nk = (half) ? (((N + 1) % 2) ? (N + 1) / 2 + 1
                                   : (N + 1) / 2) /* ceil((N + 1)/2) */
                  : (N + 1) / 2;                  /* floor((N + 1)/2) */

  if (nk == 0)
    return;

  bfam_long_real_t tworho = 2 * (N + 1) + alpha + beta + 1;
  bfam_long_real_t *tmp;

  bfam_long_real_t *restrict theta0 =
      bfam_malloc_aligned(nk * sizeof(bfam_long_real_t));
  bfam_long_real_t *restrict theta1 =
      bfam_malloc_aligned(nk * sizeof(bfam_long_real_t));
  bfam_long_real_t *restrict p0 =
      bfam_malloc_aligned(nk * sizeof(bfam_long_real_t));
  bfam_long_real_t *restrict dp0 =
      bfam_malloc_aligned(nk * sizeof(bfam_long_real_t));

  BFAM_ASSUME_ALIGNED(theta0, 32);
  BFAM_ASSUME_ALIGNED(theta1, 32);
  BFAM_ASSUME_ALIGNED(p0, 32);
  BFAM_ASSUME_ALIGNED(dp0, 32);

  /*
   * Use Gatteschi and Pittaluga's approximation for the roots of the Jacobi
   * polynomials as an initial guess.  See equation (3.19) of
   * Nicholas Hale and Alex Townsend ``Fast and Accurate Computation of
   * Gauss–Legendre and Gauss–Jacobi Quadrature Nodes and Weights'' SIAM J.
   * SCI. COMPUT. Vol. 35, No. 2, pp. A652–A674.
   */
  for (int k = nk; k > 0; --k)
  {
    int khat = (half) ? nk - k : k - 1;

    bfam_long_real_t phik =
        (2 * k + alpha - BFAM_LONG_REAL(0.5)) * BFAM_LONG_REAL_PI / tworho;

    theta1[khat] = phik +
                   1 / (tworho * tworho) *
                       ((BFAM_LONG_REAL(0.25) - alpha * alpha) * 1 /
                            BFAM_LONG_REAL_TAN(BFAM_LONG_REAL(0.5) * phik) -
                        (BFAM_LONG_REAL(0.25) - beta * beta) *
                            BFAM_LONG_REAL_TAN(BFAM_LONG_REAL(0.5) * phik));
  }

  /*
   * Use Newton's method for finding the roots of the Jacobi polynomial.
   */
  int converged = 0;
  for (int i = 0; i < MAX_ITERATIONS; ++i)
  {
    tmp = theta0;
    theta0 = theta1;
    theta1 = tmp;

    for (int k = 0; k < nk; ++k)
    {
      x[k] = BFAM_LONG_REAL_COS(theta0[k]);
    }

    bfam_jacobi_p(alpha, beta, N + 1, nk, x, p0);
    bfam_grad_jacobi_p(alpha, beta, N + 1, nk, x, dp0);

    for (int k = 0; k < nk; ++k)
    {
      theta1[k] = theta0[k] - p0[k] / (-BFAM_LONG_REAL_SIN(theta0[k]) * dp0[k]);
    }

    int diff = 0;
    for (int k = 0; k < nk; ++k)
    {
      diff += !BFAM_LONG_REAL_APPROX_EQ(theta0[k], theta1[k], 10);
    }
    if (!diff)
    {
      converged = 1;
      break;
    }
  }

  BFAM_ABORT_IF(
      !converged,
      "Newton's method does not converge when computing Jacobi Gauss points");
  /*
   * Nodes
   */
  for (int k = 0; k < nk; ++k)
  {
    x[k] = BFAM_LONG_REAL_COS(theta1[k]);
  }

  /*
   * Weights
   */
  bfam_grad_jacobi_p(alpha, beta, N + 1, nk, x, dp0);

  for (int k = 0; k < nk; ++k)
  {
    bfam_long_real_t sint = BFAM_LONG_REAL_SIN(theta1[k]);
    w[k] = tworho / (sint * sint * dp0[k] * dp0[k]);
  }

  bfam_free_aligned(theta0);
  bfam_free_aligned(theta1);
  bfam_free_aligned(p0);
  bfam_free_aligned(dp0);

  return;
}

static void bfam_jacobi_gauss_quadrature(bfam_long_real_t alpha,
                                         bfam_long_real_t beta, int N,
                                         bfam_long_real_t *restrict x,
                                         bfam_long_real_t *restrict w)
{
  int nk_floor = (N + 1) / 2; /* floor((N + 1)/2) */

  bfam_jacobi_gauss_quadrature_half(alpha, beta, N, 1, x + nk_floor,
                                    w + nk_floor);
  bfam_jacobi_gauss_quadrature_half(beta, alpha, N, 0, x, w);

  for (int k = 0; k < nk_floor; ++k)
    x[k] *= -1;

  return;
}

static void bfam_jacobi_gauss_lobatto_quadrature(bfam_long_real_t alpha,
                                                 bfam_long_real_t beta, int N,
                                                 bfam_long_real_t *restrict x,
                                                 bfam_long_real_t *restrict w)
{
  BFAM_ASSERT(N >= 1);

  x[0] = -1;
  x[N] = 1;

  if (N > 1)
  {
    bfam_jacobi_gauss_quadrature(alpha + 1, beta + 1, N - 2, x + 1, w + 1);
  }

  bfam_jacobi_p(alpha, beta, N, N + 1, x, w);
  bfam_long_real_t fac =
      (2 * N + alpha + beta + 1) / (N * (N + alpha + beta + 1));
  for (int k = 0; k < N + 1; ++k)
  {
    w[k] = fac / (w[k] * w[k]);
  }

  w[0] *= (1 + beta);
  w[N] *= (1 + alpha);

  return;
}

static void bfam_jacobi_p_vandermonde(bfam_long_real_t alpha,
                                      bfam_long_real_t beta, int N, size_t nx,
                                      bfam_long_real_t *x, bfam_long_real_t *V)
{
  for (int j = 0; j <= N; ++j)
    bfam_jacobi_p(alpha, beta, j, nx, x, V + j * nx);

  return;
}

static void bfam_grad_jacobi_p_vandermonde(bfam_long_real_t alpha,
                                           bfam_long_real_t beta, int N,
                                           size_t nx, bfam_long_real_t *x,
                                           bfam_long_real_t *V)
{
  for (int j = 0; j <= N; ++j)
    bfam_grad_jacobi_p(alpha, beta, j, nx, x, V + j * nx);

  return;
}

static void bfam_jacobi_p_interpolation(bfam_long_real_t alpha,
                                        bfam_long_real_t beta, int N, size_t nx,
                                        bfam_long_real_t *x,
                                        bfam_long_real_t *V,
                                        bfam_long_real_t *I)
{

  bfam_long_real_t *Vx =
      bfam_malloc_aligned((N + 1) * nx * sizeof(bfam_long_real_t));

  bfam_jacobi_p_vandermonde(alpha, beta, N, nx, x, Vx);

  bfam_util_forwardslash(nx, N + 1, Vx, V, I);

  bfam_free_aligned(Vx);

  return;
}

static void bfam_jacobi_p_differentiation(bfam_long_real_t alpha,
                                          bfam_long_real_t beta, int N,
                                          size_t nx, bfam_long_real_t *x,
                                          bfam_long_real_t *V,
                                          bfam_long_real_t *D)
{

  bfam_long_real_t *Vx =
      bfam_malloc_aligned((N + 1) * nx * sizeof(bfam_long_real_t));

  bfam_grad_jacobi_p_vandermonde(alpha, beta, N, nx, x, Vx);

  bfam_util_forwardslash(nx, N + 1, Vx, V, D);

  bfam_free_aligned(Vx);

  return;
}

static void bfam_jacobi_p_mass(bfam_long_real_t alpha, bfam_long_real_t beta,
                               int N, bfam_long_real_t *V, bfam_long_real_t *M)
{
  bfam_long_real_t *I =
      bfam_malloc_aligned((N + 1) * (N + 1) * sizeof(bfam_long_real_t));

  bfam_long_real_t *invV =
      bfam_malloc_aligned((N + 1) * (N + 1) * sizeof(bfam_long_real_t));

  bfam_long_real_t *invVT =
      bfam_malloc_aligned((N + 1) * (N + 1) * sizeof(bfam_long_real_t));

  for (int i = 0; i < (N + 1) * (N + 1); ++i)
    I[i] = 0;

  for (int i = 0; i < (N + 1) * (N + 1); ++i)
    M[i] = 0;

  for (int i = 0; i <= N; ++i)
    I[(N + 1) * i + i] = 1;

  bfam_util_backslash(N + 1, N + 1, V, I, invV);

  bfam_util_mtranspose(N + 1, N + 1, invV, N + 1, invVT, N + 1);

  bfam_util_mmmult(N + 1, N + 1, N + 1, invVT, N + 1, invV, N + 1, M, N + 1);

  bfam_free_aligned(I);
  bfam_free_aligned(invV);
  bfam_free_aligned(invVT);

  return;
}

// }}}

// {{{ kron

/** $y += (I \otimes A) x$
 *
 * \param [in]  N number of rows and columns of $A$
 * \param [in]  A column first representation of $A$
 * \param [in]  x vector $x$
 * \param [in,out] y vector $y$
 */
#define BFAM_KRON_IXA_PE(N, A, x, y)                                           \
  do                                                                           \
  {                                                                            \
    for (int bfam_kron_k = 0; bfam_kron_k < (N); ++bfam_kron_k)                \
      for (int bfam_kron_i = 0; bfam_kron_i < (N); ++bfam_kron_i)              \
        for (int bfam_kron_j = 0; bfam_kron_j < (N); ++bfam_kron_j)            \
          (y)[(N)*bfam_kron_k + bfam_kron_j] +=                                \
              (A)[(N)*bfam_kron_i + bfam_kron_j] *                             \
              (x)[(N)*bfam_kron_k + bfam_kron_i];                              \
  } while (0)

/** $y += (I \otimes A^T) x$.
 *
 * \param [in]  N number of rows and columns of $A$
 * \param [in]  A column first representation of $A$
 * \param [in]  x vector $x$
 * \param [out] y vector $y$
 */
#define BFAM_KRON_IXAT_PE(N, A, x, y)                                          \
  do                                                                           \
  {                                                                            \
    for (int bfam_kron_k = 0; bfam_kron_k < (N); ++bfam_kron_k)                \
      for (int bfam_kron_j = 0; bfam_kron_j < (N); ++bfam_kron_j)              \
        for (int bfam_kron_i = 0; bfam_kron_i < (N); ++bfam_kron_i)            \
          (y)[(N)*bfam_kron_k + bfam_kron_j] +=                                \
              (A)[(N)*bfam_kron_j + bfam_kron_i] *                             \
              (x)[(N)*bfam_kron_k + bfam_kron_i];                              \
  } while (0)

/** $y += (A \otimes I) x$.
 *
 * \param [in]  N number of rows and columns of $A$
 * \param [in]  A column first representation of $A$
 * \param [in]  x vector $x$
 * \param [out] y vector $y$
 */
#define BFAM_KRON_AXI_PE(N, A, x, y)                                           \
  do                                                                           \
  {                                                                            \
    for (int bfam_kron_i = 0; bfam_kron_i < (N); ++bfam_kron_i)                \
      for (int bfam_kron_k = 0; bfam_kron_k < (N); ++bfam_kron_k)              \
        for (int bfam_kron_j = 0; bfam_kron_j < (N); ++bfam_kron_j)            \
          (y)[(N)*bfam_kron_k + bfam_kron_j] +=                                \
              (A)[(N)*bfam_kron_i + bfam_kron_k] *                             \
              (x)[(N)*bfam_kron_i + bfam_kron_j];                              \
  } while (0)

/** $y += (A^T \otimes I) x$.
 *
 * \param [in]  N number of rows and columns of $A$
 * \param [in]  A column first representation of $A$
 * \param [in]  x vector $x$
 * \param [out] y vector $y$
 */
#define BFAM_KRON_ATXI_PE(N, A, x, y)                                          \
  do                                                                           \
  {                                                                            \
    for (int bfam_kron_k = 0; bfam_kron_k < (N); ++bfam_kron_k)                \
      for (int bfam_kron_i = 0; bfam_kron_i < (N); ++bfam_kron_i)              \
        for (int bfam_kron_j = 0; bfam_kron_j < (N); ++bfam_kron_j)            \
          (y)[(N)*bfam_kron_k + bfam_kron_j] +=                                \
              (A)[(N)*bfam_kron_k + bfam_kron_i] *                             \
              (x)[(N)*bfam_kron_i + bfam_kron_j];                              \
  } while (0)

/** $y = (I \otimes A) x$.
 *
 * \param [in]  N number of rows and columns of $A$
 * \param [in]  A column first representation of $A$
 * \param [in]  x vector $x$
 * \param [out] y vector $y$
 */
#define BFAM_KRON_IXA(N, A, x, y)                                              \
  do                                                                           \
  {                                                                            \
    for (int bfam_kron_n = 0; bfam_kron_n < (N) * (N); ++bfam_kron_n)          \
      (y)[bfam_kron_n] = 0;                                                    \
    BFAM_KRON_IXA_PE(N, A, x, y);                                              \
  } while (0)

/** $y = (I \otimes A^T) x$.
 *
 * \param [in]  N number of rows and columns of $A$
 * \param [in]  A column first representation of $A$
 * \param [in]  x vector $x$
 * \param [out] y vector $y$
 */
#define BFAM_KRON_IXAT(N, A, x, y)                                             \
  do                                                                           \
  {                                                                            \
    for (int bfam_kron_n = 0; bfam_kron_n < (N) * (N); ++bfam_kron_n)          \
      (y)[bfam_kron_n] = 0;                                                    \
    BFAM_KRON_IXAT_PE(N, A, x, y);                                             \
  } while (0)

/** $y = (A \otimes I) x$.
 *
 * \param [in]  N number of rows and columns of $A$
 * \param [in]  A column first representation of $A$
 * \param [in]  x vector $x$
 * \param [out] y vector $y$
 */
#define BFAM_KRON_AXI(N, A, x, y)                                              \
  do                                                                           \
  {                                                                            \
    for (int bfam_kron_n = 0; bfam_kron_n < (N) * (N); ++bfam_kron_n)          \
      (y)[bfam_kron_n] = 0;                                                    \
    BFAM_KRON_AXI_PE(N, A, x, y);                                              \
  } while (0)

/** $y = (A^T \otimes I) x$.
 *
 * \param [in]  N number of rows and columns of $A$
 * \param [in]  A column first representation of $A$
 * \param [in]  x vector $x$
 * \param [out] y vector $y$
 */
#define BFAM_KRON_ATXI(N, A, x, y)                                             \
  do                                                                           \
  {                                                                            \
    for (int bfam_kron_n = 0; bfam_kron_n < (N) * (N); ++bfam_kron_n)          \
      (y)[bfam_kron_n] = 0;                                                    \
    BFAM_KRON_ATXI_PE(N, A, x, y);                                             \
  } while (0)

/** $y = a \dot\times x$
 *
 * \param [in]  N number of elements of $a$, $x$, $y$
 * \param [in]  a vector $a$
 * \param [in]  x vector $x$
 * \param [out] y vector $y$
 */
#define BFAM_DOT_AX(N, a, x, y)                                                \
  do                                                                           \
  {                                                                            \
    for (int bfam_dot_n = 0; bfam_dot_n < (N); ++bfam_dot_n)                   \
      (y)[bfam_dot_n] = (a)[bfam_dot_n] * (x)[bfam_dot_n];                     \
  } while (0)

/** $y = a \dot\times b \dot\times x$
 *
 * \param [in]  N number of elements of $a$, $b$, $x$, $y$
 * \param [in]  a vector $a$
 * \param [in]  b vector $b$
 * \param [in]  x vector $x$
 * \param [out] y vector $y$
 */
#define BFAM_DOT_ABX(N, a, b, x, y)                                            \
  do                                                                           \
  {                                                                            \
    for (int bfam_dot_n = 0; bfam_dot_n < (N); ++bfam_dot_n)                   \
      (y)[bfam_dot_n] = (a)[bfam_dot_n] * (b)[bfam_dot_n] * (x)[bfam_dot_n];   \
  } while (0)

/** $y += a \dot\times x$
 *
 * \param [in]  N number of elements of $a$, $x$, $y$
 * \param [in]  a vector $a$
 * \param [in]  x vector $x$
 * \param [out] y vector $y$
 */
#define BFAM_DOT_AX_PE(N, a, x, y)                                             \
  do                                                                           \
  {                                                                            \
    for (int bfam_dot_n = 0; bfam_dot_n < (N); ++bfam_dot_n)                   \
      (y)[bfam_dot_n] += (a)[bfam_dot_n] * (x)[bfam_dot_n];                    \
  } while (0)

/** $y -= a \dot\times x$
 *
 * \param [in]  N number of elements of $a$, $x$, $y$
 * \param [in]  a vector $a$
 * \param [in]  x vector $x$
 * \param [out] y vector $y$
 */
#define BFAM_DOT_AX_ME(N, a, x, y)                                             \
  do                                                                           \
  {                                                                            \
    for (int bfam_dot_n = 0; bfam_dot_n < (N); ++bfam_dot_n)                   \
      (y)[bfam_dot_n] -= (a)[bfam_dot_n] * (x)[bfam_dot_n];                    \
  } while (0)

/** $y -= 2*a \dot\times x$
 *
 * \param [in]  N number of elements of $a$, $x$, $y$
 * \param [in]  a vector $a$
 * \param [in]  x vector $x$
 * \param [out] y vector $y$
 */
#define BFAM_DOT_2AX_ME(N, a, x, y)                                            \
  do                                                                           \
  {                                                                            \
    for (int bfam_dot_n = 0; bfam_dot_n < (N); ++bfam_dot_n)                   \
      (y)[bfam_dot_n] -= 2 * (a)[bfam_dot_n] * (x)[bfam_dot_n];                \
  } while (0)

/** $y += a \dot\times b \dot\times x$
 *
 * \param [in]  N number of elements of $a$, $b$, $x$, $y$
 * \param [in]  a vector $a$
 * \param [in]  b vector $b$
 * \param [in]  x vector $x$
 * \param [out] y vector $y$
 */
#define BFAM_DOT_ABX_PE(N, a, b, x, y)                                         \
  do                                                                           \
  {                                                                            \
    for (int bfam_dot_n = 0; bfam_dot_n < (N); ++bfam_dot_n)                   \
      (y)[bfam_dot_n] += (a)[bfam_dot_n] * (b)[bfam_dot_n] * (x)[bfam_dot_n];  \
  } while (0)

/** $y += (A \otimes B) \dot\times C \dot\times D \dot\times x$
 *
 * \param [in]  N number of elements of $a$ and $b$
 * \param [in]  a vector $a$
 * \param [in]  b vector $b$
 * \param [in]  c vector $c$
 * \param [in]  d vector $d$
 * \param [in]  x vector $x$
 * \param [out] y vector $y$
 */
#define BFAM_KRON_AB_DOT_CD_PE(N, a, b, c, d, x, y)                            \
  do                                                                           \
  {                                                                            \
    for (int bfam_kron_i = 0; bfam_kron_i < (N); ++bfam_kron_i)                \
      for (int bfam_kron_j = 0; bfam_kron_j < (N); ++bfam_kron_j)              \
        (y)[(N)*bfam_kron_i + bfam_kron_j] +=                                  \
            (a)[bfam_kron_i] * (b)[bfam_kron_j] *                              \
            (c)[(N)*bfam_kron_i + bfam_kron_j] *                               \
            (d)[(N)*bfam_kron_i + bfam_kron_j] *                               \
            (x)[(N)*bfam_kron_i + bfam_kron_j];                                \
  } while (0)

/** $y = (A \otimes B) \dot\times x$
 *
 * \param [in]  N number of elements of $a$ and $b$
 * \param [in]  a vector $a$
 * \param [in]  b vector $b$
 * \param [in]  x vector $x$
 * \param [out] y vector $y$
 */
#define BFAM_KRON_AB_DOT_C(N, a, b, c, x, y)                                   \
  do                                                                           \
  {                                                                            \
    for (int bfam_kron_i = 0; bfam_kron_i < (N); ++bfam_kron_i)                \
      for (int bfam_kron_j = 0; bfam_kron_j < (N); ++bfam_kron_j)              \
        (y)[(N)*bfam_kron_i + bfam_kron_j] =                                   \
            (a)[bfam_kron_i] * (b)[bfam_kron_j] *                              \
            (c)[(N)*bfam_kron_i + bfam_kron_j] *                               \
            (x)[(N)*bfam_kron_i + bfam_kron_j];                                \
  } while (0)

/** $z = x +  a(B \otimes C) \dot\times D \dot\times y$
 *
 * \param [in]  N number of elements of $b$ and $c$
 * \param [in]  a scalara
 * \param [in]  a vector $a$
 * \param [in]  b vector $b$
 * \param [in]  c vector $c$
 * \param [in]  d vector $d$
 * \param [in]  x vector $x$
 * \param [in]  y vector $y$
 * \param [out] z vector $z$
 */
#define BFAM_KRON_A_BC_DOT_D_PE(N, a, b, c, d, x, y, z)                        \
  do                                                                           \
  {                                                                            \
    for (int bfam_kron_i = 0; bfam_kron_i < (N); ++bfam_kron_i)                \
      for (int bfam_kron_j = 0; bfam_kron_j < (N); ++bfam_kron_j)              \
        (z)[(N)*bfam_kron_i + bfam_kron_j] =                                   \
            (x)[(N)*bfam_kron_i + bfam_kron_j] +                               \
            a * (b)[bfam_kron_i] * (c)[bfam_kron_j] *                          \
                (d)[(N)*bfam_kron_i + bfam_kron_j] *                           \
                (y)[(N)*bfam_kron_i + bfam_kron_j];                            \
  } while (0)

/** $y = (A \otimes I \otimes I) x$.
 *
 * \param [in]  N number of rows and columns of $A$
 * \param [in]  A column first representation of $A$
 * \param [in]  x vector $x$
 * \param [out] y vector $y$
 */
#define BFAM_KRON_AXIXI(N, A, x, y)                                            \
  do                                                                           \
  {                                                                            \
    for (int bfam_kron_n = 0; bfam_kron_n < (N) * (N) * (N); ++bfam_kron_n)    \
      (y)[bfam_kron_n] = 0;                                                    \
    BFAM_KRON_AXIXI_PE(N, A, x, y);                                            \
  } while (0)

/** $y = (A^T \otimes I \otimes I) x$.
 *
 * \param [in]  N number of rows and columns of $A$
 * \param [in]  A column first representation of $A$
 * \param [in]  x vector $x$
 * \param [out] y vector $y$
 */
#define BFAM_KRON_ATXIXI(N, A, x, y)                                           \
  do                                                                           \
  {                                                                            \
    for (int bfam_kron_n = 0; bfam_kron_n < (N) * (N) * (N); ++bfam_kron_n)    \
      (y)[bfam_kron_n] = 0;                                                    \
    BFAM_KRON_ATXIXI_PE(N, A, x, y);                                           \
  } while (0)

/** $y = (I \otimes A \otimes I) x$.
 *
 * \param [in]  N number of rows and columns of $A$
 * \param [in]  A column first representation of $A$
 * \param [in]  x vector $x$
 * \param [out] y vector $y$
 */
#define BFAM_KRON_IXAXI(N, A, x, y)                                            \
  do                                                                           \
  {                                                                            \
    for (int bfam_kron_n = 0; bfam_kron_n < (N) * (N) * (N); ++bfam_kron_n)    \
      (y)[bfam_kron_n] = 0;                                                    \
    BFAM_KRON_IXAXI_PE(N, A, x, y);                                            \
  } while (0)

/** $y = (I \otimes A^T \otimes I) x$.
 *
 * \param [in]  N number of rows and columns of $A$
 * \param [in]  A column first representation of $A$
 * \param [in]  x vector $x$
 * \param [out] y vector $y$
 */
#define BFAM_KRON_IXATXI(N, A, x, y)                                           \
  do                                                                           \
  {                                                                            \
    for (int bfam_kron_n = 0; bfam_kron_n < (N) * (N) * (N); ++bfam_kron_n)    \
      (y)[bfam_kron_n] = 0;                                                    \
    BFAM_KRON_IXATXI_PE(N, A, x, y);                                           \
  } while (0)

/** $y = (I \otimes I \otimes A) x$.
 *
 * \param [in]  N number of rows and columns of $A$
 * \param [in]  A column first representation of $A$
 * \param [in]  x vector $x$
 * \param [out] y vector $y$
 */
#define BFAM_KRON_IXIXA(N, A, x, y)                                            \
  do                                                                           \
  {                                                                            \
    for (int bfam_kron_n = 0; bfam_kron_n < (N) * (N) * (N); ++bfam_kron_n)    \
      (y)[bfam_kron_n] = 0;                                                    \
    BFAM_KRON_IXIXA_PE(N, A, x, y);                                            \
  } while (0)

/** $y = (I \otimes I \otimes A^T) x$.
 *
 * \param [in]  N number of rows and columns of $A$
 * \param [in]  A column first representation of $A$
 * \param [in]  x vector $x$
 * \param [out] y vector $y$
 */
#define BFAM_KRON_IXIXAT(N, A, x, y)                                           \
  do                                                                           \
  {                                                                            \
    for (int bfam_kron_n = 0; bfam_kron_n < (N) * (N) * (N); ++bfam_kron_n)    \
      (y)[bfam_kron_n] = 0;                                                    \
    BFAM_KRON_IXIXAT_PE(N, A, x, y);                                           \
  } while (0)

/** $y += (A \otimes I \otimes I) x$.
 *
 * \param [in]  N number of rows and columns of $A$
 * \param [in]  A column first representation of $A$
 * \param [in]  x vector $x$
 * \param [out] y vector $y$
 */
#define BFAM_KRON_AXIXI_PE(N, A, x, y)                                         \
  do                                                                           \
  {                                                                            \
    for (int bfam_kron_i = 0; bfam_kron_i < (N); ++bfam_kron_i)                \
      for (int bfam_kron_l = 0; bfam_kron_l < (N); ++bfam_kron_l)              \
        for (int bfam_kron_k = 0; bfam_kron_k < (N); ++bfam_kron_k)            \
          for (int bfam_kron_j = 0; bfam_kron_j < (N); ++bfam_kron_j)          \
            (y)[(N) * ((N)*bfam_kron_l + bfam_kron_k) + bfam_kron_j] +=        \
                (A)[(N)*bfam_kron_i + bfam_kron_l] *                           \
                (x)[(N) * ((N)*bfam_kron_i + bfam_kron_k) + bfam_kron_j];      \
  } while (0)

/** $y += (A^T \otimes I \otimes I) x$.
 *
 * \param [in]  N number of rows and columns of $A$
 * \param [in]  A column first representation of $A$
 * \param [in]  x vector $x$
 * \param [out] y vector $y$
 */
#define BFAM_KRON_ATXIXI_PE(N, A, x, y)                                        \
  do                                                                           \
  {                                                                            \
    for (int bfam_kron_i = 0; bfam_kron_i < (N); ++bfam_kron_i)                \
      for (int bfam_kron_l = 0; bfam_kron_l < (N); ++bfam_kron_l)              \
        for (int bfam_kron_k = 0; bfam_kron_k < (N); ++bfam_kron_k)            \
          for (int bfam_kron_j = 0; bfam_kron_j < (N); ++bfam_kron_j)          \
            (y)[(N) * ((N)*bfam_kron_l + bfam_kron_k) + bfam_kron_j] +=        \
                (A)[(N)*bfam_kron_l + bfam_kron_i] *                           \
                (x)[(N) * ((N)*bfam_kron_i + bfam_kron_k) + bfam_kron_j];      \
  } while (0)

/** $y += (I \otimes A \otimes I) x$.
 *
 * \param [in]  N number of rows and columns of $A$
 * \param [in]  A column first representation of $A$
 * \param [in]  x vector $x$
 * \param [out] y vector $y$
 */
#define BFAM_KRON_IXAXI_PE(N, A, x, y)                                         \
  do                                                                           \
  {                                                                            \
    for (int bfam_kron_i = 0; bfam_kron_i < (N); ++bfam_kron_i)                \
      for (int bfam_kron_l = 0; bfam_kron_l < (N); ++bfam_kron_l)              \
        for (int bfam_kron_k = 0; bfam_kron_k < (N); ++bfam_kron_k)            \
          for (int bfam_kron_j = 0; bfam_kron_j < (N); ++bfam_kron_j)          \
            (y)[(N) * ((N)*bfam_kron_l + bfam_kron_k) + bfam_kron_j] +=        \
                (A)[(N)*bfam_kron_i + bfam_kron_k] *                           \
                (x)[(N) * ((N)*bfam_kron_l + bfam_kron_i) + bfam_kron_j];      \
  } while (0)

/** $y += (I \otimes A^T \otimes I) x$.
 *
 * \param [in]  N number of rows and columns of $A$
 * \param [in]  A column first representation of $A$
 * \param [in]  x vector $x$
 * \param [out] y vector $y$
 */
#define BFAM_KRON_IXATXI_PE(N, A, x, y)                                        \
  do                                                                           \
  {                                                                            \
    for (int bfam_kron_i = 0; bfam_kron_i < (N); ++bfam_kron_i)                \
      for (int bfam_kron_l = 0; bfam_kron_l < (N); ++bfam_kron_l)              \
        for (int bfam_kron_k = 0; bfam_kron_k < (N); ++bfam_kron_k)            \
          for (int bfam_kron_j = 0; bfam_kron_j < (N); ++bfam_kron_j)          \
            (y)[(N) * ((N)*bfam_kron_l + bfam_kron_k) + bfam_kron_j] +=        \
                (A)[(N)*bfam_kron_k + bfam_kron_i] *                           \
                (x)[(N) * ((N)*bfam_kron_l + bfam_kron_i) + bfam_kron_j];      \
  } while (0)

/** $y += (I \otimes I \otimes A) x$.
 *
 * \param [in]  N number of rows and columns of $A$
 * \param [in]  A column first representation of $A$
 * \param [in]  x vector $x$
 * \param [out] y vector $y$
 */
#define BFAM_KRON_IXIXA_PE(N, A, x, y)                                         \
  do                                                                           \
  {                                                                            \
    for (int bfam_kron_i = 0; bfam_kron_i < (N); ++bfam_kron_i)                \
      for (int bfam_kron_l = 0; bfam_kron_l < (N); ++bfam_kron_l)              \
        for (int bfam_kron_k = 0; bfam_kron_k < (N); ++bfam_kron_k)            \
          for (int bfam_kron_j = 0; bfam_kron_j < (N); ++bfam_kron_j)          \
            (y)[(N) * ((N)*bfam_kron_l + bfam_kron_k) + bfam_kron_j] +=        \
                (A)[(N)*bfam_kron_i + bfam_kron_j] *                           \
                (x)[(N) * ((N)*bfam_kron_l + bfam_kron_k) + bfam_kron_i];      \
  } while (0)

/** $y += (I \otimes I \otimes A^T) x$.
 *
 * \param [in]  N number of rows and columns of $A$
 * \param [in]  A column first representation of $A$
 * \param [in]  x vector $x$
 * \param [out] y vector $y$
 */
#define BFAM_KRON_IXIXAT_PE(N, A, x, y)                                        \
  do                                                                           \
  {                                                                            \
    for (int bfam_kron_i = 0; bfam_kron_i < (N); ++bfam_kron_i)                \
      for (int bfam_kron_l = 0; bfam_kron_l < (N); ++bfam_kron_l)              \
        for (int bfam_kron_k = 0; bfam_kron_k < (N); ++bfam_kron_k)            \
          for (int bfam_kron_j = 0; bfam_kron_j < (N); ++bfam_kron_j)          \
            (y)[(N) * ((N)*bfam_kron_l + bfam_kron_k) + bfam_kron_j] +=        \
                (A)[(N)*bfam_kron_j + bfam_kron_i] *                           \
                (x)[(N) * ((N)*bfam_kron_l + bfam_kron_k) + bfam_kron_i];      \
  } while (0)

/* // }}} */

// {{{ domain pxest

/* conversions between p4est and standard orientation */
static const int8_t bfam_p8est_FToF_code[6][6] = {
    {0, 1, 1, 0, 0, 1}, {2, 0, 0, 1, 1, 0}, {2, 0, 0, 1, 1, 0},
    {0, 2, 2, 0, 0, 1}, {0, 2, 2, 0, 0, 1}, {2, 0, 0, 2, 2, 0}};

static const int8_t bfam_p8est_code_to_perm[3][4] = {
    {1, 2, 5, 6}, {0, 3, 4, 7}, {0, 4, 3, 7}};

#if DIM == 3
static const int8_t bfam_p8est_perm_to_order[8][4] = {
    {0, 1, 2, 3}, {0, 2, 1, 3}, {1, 0, 3, 2}, {1, 3, 0, 2},
    {2, 0, 3, 1}, {2, 3, 0, 1}, {3, 1, 2, 0}, {3, 2, 1, 0}};

static const int8_t bfam_p8est_perm_to_order_inv[8][4] = {
    {0, 1, 2, 3}, {0, 2, 1, 3}, {1, 0, 3, 2}, {2, 0, 3, 1},
    {1, 3, 0, 2}, {2, 3, 0, 1}, {3, 1, 2, 0}, {3, 2, 1, 0}};
#endif

/** BFAM_P8EST_ORIENTATION
 * f    my face
 * nf   neighbor face
 * o    p8est orientation code
 */
#define BFAM_P8EST_ORIENTATION(f, nf, o)                                       \
  (bfam_p8est_code_to_perm[bfam_p8est_FToF_code[f][nf]][o])

/*
 * orientation 0:
 *   2---3     2---3
 *   |   | --> |   |
 *   0---1     0---1
 *   same:
 *   (a,b) --> (a,b)
 *
 * orientation 1:
 *   2---3     1---3
 *   |   | --> |   |
 *   0---1     0---2
 *   switch indices:
 *   (a,b) --> (b,a)
 *
 * orientation 2:
 *   2---3     3---2
 *   |   | --> |   |
 *   0---1     1---0
 *   reverse first index:
 *   (a,b) --> (Na-a,b)
 *
 * orientation 3:
 *   2---3     0---2
 *   |   | --> |   |
 *   0---1     1---3
 *   reverse first index and switch:
 *   (a,b) --> (b,Na-a)
 *
 * orientation 4:
 *   2---3     3---1
 *   |   | --> |   |
 *   0---1     2---0
 *   reverse second index and switch:
 *   (a,b) --> (Nb-b,a)
 *
 * orientation 5:
 *   2---3     0---1
 *   |   | --> |   |
 *   0---1     2---3
 *   reverse second index:
 *   (a,b) --> (a,Nb-b)
 *
 * orientation 6:
 *   2---3     2---0
 *   |   | --> |   |
 *   0---1     3---1
 *   reverse both and switch:
 *   (a,b) --> (Nb-b,Na-a)
 *
 * orientation 7:
 *   2---3     1---0
 *   |   | --> |   |
 *   0---1     3---2
 *   reverse both:
 *   (a,b) --> (Na-a,Nb-b)
 */

#if DIM == 2
#define BFAM_PXEST_RELATIVE_ORIENTATIONS 2
#define BDIM 1
#define BFAM_PXEST_CONNECT P4EST_CONNECT_FULL
#define BFAM_PXEST_ORIENTATION(n, nf, o) (o)
#define BFAM_PXEST_BOHTONH(bo, h) (((bo) + (h)) % (2))
#define BFAM_PXEST_BOHTONH_INV(bo, h) BFAM_PXEST_BOHTONH(bo, h)
#define BFAM_D3_AP(A1, A2) (A1)
#elif DIM == 3
#define BFAM_PXEST_RELATIVE_ORIENTATIONS 4
#define BFAM_PXEST_ORIENTATION(n, nf, o) BFAM_P8EST_ORIENTATION(n, nf, o)
#define BDIM 2
#define BFAM_PXEST_CONNECT P8EST_CONNECT_FULL
#define BFAM_PXEST_BOHTONH(bo, h) (bfam_p8est_perm_to_order[(bo)][(h)])
#define BFAM_PXEST_BOHTONH_INV(bo, h) (bfam_p8est_perm_to_order_inv[(bo)][(h)])
#define BFAM_D3_AP(A1, A2) (A1 A2)
#else
#error "Bad Dimension"
#endif

/** Callback for p4est user data init function.
 */
static void bfam_domain_pxest_init_callback(p4est_t *p4est,
                                            p4est_topidx_t which_tree,
                                            p4est_quadrant_t *quadrant)
{

  bfam_pxest_user_data_t *ud = quadrant->p.user_data;

  ud->flags = 0;
  ud->N = -1;
  ud->Nold = -1;
  ud->subd_id = -1;
  ud->elem_id = -1;
  ud->root_id = BFAM_DEFAULT_SUBDOMAIN_ROOT;
  ud->glue_id[0] = -1;
  ud->glue_id[1] = -1;
  ud->glue_id[2] = -1;
  ud->glue_id[3] = -1;
#if DIM == 3
  ud->glue_id[4] = -1;
  ud->glue_id[5] = -1;
#endif
}

/** initializes a domain
 *
 * \param [in,out] domain          pointer to the pxest managed domain
 * \param [in]     domComm         pointer to the communicator for the domain
 * \param [in]     conn            pointer to the pxest connectivity for the
 *                                 domain
 * \param [in]     min_quadrants   Minimum initial quadrants per processor.
 *                                 Makes the refinement pattern
 *                                 mpisize-specific.
 * \param [in]     min_level       The forest is refined at least to this level.
 *                                 May be negative or 0, then it has no effect.
 * \param [in]     fill_uniform    If true, fill the forest with a uniform mesh
 *                                 instead of the coarsest possible one.
 *                                 The latter is partition-specific so that
 *                                 is usually not a good idea.
 */
static void bfam_domain_pxest_init_ext(bfam_domain_pxest_t *domain,
                                       MPI_Comm domComm,
                                       p4est_connectivity_t *conn,
                                       p4est_locidx_t min_quadrants,
                                       int min_level, int fill_uniform)
{
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wmissing-field-initializers"
#endif
  bfam_pxest_user_data_t default_user_data = {0};
#if defined(__clang__)
#pragma clang diagnostic pop
#endif
  bfam_domain_init(&domain->base, domComm);

  domain->conn = conn;
  domain->pxest =
      p4est_new_ext(domComm, conn, min_quadrants, min_level, fill_uniform,
                    sizeof(bfam_pxest_user_data_t),
                    bfam_domain_pxest_init_callback, &default_user_data);
  domain->N2N = bfam_malloc(sizeof(bfam_dictionary_t));
  bfam_dictionary_init(domain->N2N);

  domain->dgx_ops = bfam_malloc(sizeof(bfam_dictionary_t));
  bfam_dictionary_init(domain->dgx_ops);
}

/* Domain managed by pxest based functions */

bfam_domain_pxest_t *bfam_domain_pxest_new_ext(MPI_Comm domComm,
                                               p4est_connectivity_t *conn,
                                               p4est_locidx_t min_quadrants,
                                               int min_level, int fill_uniform)
{
  bfam_domain_pxest_t *newDomain = bfam_malloc(sizeof(bfam_domain_pxest_t));
  bfam_domain_pxest_init_ext(newDomain, domComm, conn, min_quadrants, min_level,
                             fill_uniform);
  return newDomain;
}

/** free the dgx_ops data put in the dgx_ops dictionary
 */
static int bfam_subdomain_dgx_clear_dgx_ops_dict(const char *key, void *val,
                                                 void *args)
{
  bfam_free_aligned(val);
  return 1;
}

/** free the interpolator data put in the interpolation dictionary
 */

/*
 * -+-            -+-
 *  |  Coasen ->   |
 * -+-             |
 *  |  <- Refine   |
 * -+-            -+-
 */
typedef struct
{
  int N_src;
  int N_dst;
  bfam_locidx_t num_prj;
  bfam_real_t **prj; /* array of projection operators;
                      * no refinement (NULL is no change in order)
                      * coarsen from bottom
                      * coarsen from top
                      * refine  to   bottom
                      * refine  to   top
                      */
  bfam_real_t *
      *mass_prj; /* array of exact mass projection operators:
                  * (i.e., the mass matrix for the destination is multiplied
                  * times the projection: M * Pr * q)
                  * no refinement (i.e., just mass)
                  * coarsen from bottom
                  * coarsen from top
                  * refine  to   bottom
                  * refine  to   top
                  */
  bfam_real_t *
      *wi_mass_prj; /* array of LGL mass inv exact mass projection operators:
                  * (i.e., the LGL diag mass and exact mass matrix for the
                  * destination
                  * are multiplied times the projection: diag(wi) * M * Pr * q)
                  * no refinement (i.e., just mass)
                  * coarsen from bottom
                  * coarsen from top
                  * refine  to   bottom
                  * refine  to   top
                  */
} bfam_subdomain_dgx_interpolator_t;
static int bfam_subdomain_dgx_clear_interpolation_dict(const char *key,
                                                       void *val, void *args)
{
  bfam_subdomain_dgx_interpolator_t *interp =
      (bfam_subdomain_dgx_interpolator_t *)val;
  for (bfam_locidx_t k = 0; k < interp->num_prj; k++)
    if (interp->prj[k])
      bfam_free_aligned(interp->prj[k]);
  bfam_free(interp->prj);
  for (bfam_locidx_t k = 0; k < interp->num_prj; k++)
    if (interp->mass_prj[k])
      bfam_free_aligned(interp->mass_prj[k]);
  bfam_free(interp->mass_prj);
  for (bfam_locidx_t k = 0; k < interp->num_prj; k++)
    if (interp->wi_mass_prj[k])
      bfam_free_aligned(interp->wi_mass_prj[k]);
  bfam_free(interp->wi_mass_prj);
  bfam_free(val);
  return 1;
}

void bfam_domain_pxest_free(bfam_domain_pxest_t *domain)
{
  /* Memory we don't manage */
  domain->conn = NULL;

  /* Memory we do manage */
  if (domain->pxest)
    p4est_destroy(domain->pxest);
  domain->pxest = NULL;

  if (domain->N2N)
  {
    bfam_dictionary_allprefixed_ptr(
        domain->N2N, "", bfam_subdomain_dgx_clear_interpolation_dict, NULL);
    bfam_dictionary_clear(domain->N2N);
    bfam_free(domain->N2N);
  }
  domain->N2N = NULL;

  if (domain->dgx_ops)
  {
    bfam_dictionary_allprefixed_ptr(
        domain->dgx_ops, "", bfam_subdomain_dgx_clear_dgx_ops_dict, NULL);
    bfam_dictionary_clear(domain->dgx_ops);
    bfam_free(domain->dgx_ops);
  }
  domain->dgx_ops = NULL;

  bfam_domain_free(&domain->base);
}

static void bfam_domain_pxest_dgx_print_stats(bfam_domain_pxest_t *domain)
{
  bfam_domain_t *dbase = &domain->base;
  bfam_subdomain_t **subdomains =
      bfam_malloc(dbase->num_subdomains * sizeof(bfam_subdomain_t **));

  bfam_locidx_t num_subdomains = 0;

  const char *volume[] = {"_volume", NULL};

  bfam_domain_get_subdomains(dbase, BFAM_DOMAIN_AND, volume,
                             dbase->num_subdomains, subdomains,
                             &num_subdomains);

#define GRID_PTS 0
#define ELEMENTS 1
#define NUM_VALS 2

  bfam_gloidx_t vals_loc[NUM_VALS];
  bfam_gloidx_t vals_glo[NUM_VALS];

  for (size_t i = 0; i < NUM_VALS; ++i)
  {
    vals_loc[i] = 0;
    vals_glo[i] = 0;
  }

  for (bfam_locidx_t s = 0; s < num_subdomains; ++s)
  {
    bfam_subdomain_dgx_t *sub = (bfam_subdomain_dgx_t *)subdomains[s];

    vals_loc[GRID_PTS] += sub->K * sub->Np;
    vals_loc[ELEMENTS] += sub->K;
  }

  BFAM_MPI_CHECK(MPI_Reduce(vals_loc, vals_glo, (int)NUM_VALS, BFAM_GLOIDX_MPI,
                            MPI_SUM, 0, dbase->comm));

  BFAM_ROOT_INFO("Domain Stats --- Elements: %jd",
                 (intmax_t)vals_glo[ELEMENTS]);
  BFAM_ROOT_INFO("Domain Stats --- Grid Pts: %jd",
                 (intmax_t)vals_glo[GRID_PTS]);

#undef GRID_PTS
#undef ELEMENTS
#undef NUM_VALS

  bfam_free(subdomains);
}

#ifdef BFAM_DEBUG
static void
bfam_domain_pxest_tree_to_glueid_check(p4est_connectivity_t *conn,
                                       const p4est_locidx_t *tree_to_glueid)
{
  for (p4est_topidx_t t = 0; t < conn->num_trees; ++t)
  {
    for (int f = 0; f < P4EST_FACES; ++f)
    {
      p4est_topidx_t nt = conn->tree_to_tree[P4EST_FACES * t + f];
      int nf = (int)conn->tree_to_face[P4EST_FACES * t + f] % P4EST_FACES;

      bfam_locidx_t id = tree_to_glueid[P4EST_FACES * t + f];
      bfam_locidx_t nid = tree_to_glueid[P4EST_FACES * nt + nf];

      BFAM_ABORT_IF_NOT(
          id == nid, "tree_to_glueid invalid for (%jd,%d)->%jd!=%jd<-(%jd,%d)",
          (intmax_t)t, f, (intmax_t)id, (intmax_t)nid, (intmax_t)nt, nf);
    }
  }
}
#endif

typedef struct bfam_domain_glueid_iter_data
{
  const bfam_locidx_t *tree_to_glueid;
  bfam_locidx_t *quad_to_glueid;
} bfam_domain_glueid_iter_data_t;

static void set_quad_to_glueid(p4est_iter_face_info_t *info, void *arg)
{
  bfam_domain_glueid_iter_data_t *data = (bfam_domain_glueid_iter_data_t *)arg;

  int limit = (int)info->sides.elem_count;
  p4est_iter_face_side_t *fside;

  if (!info->tree_boundary)
    return;

  for (int i = 0; i < limit; ++i)
  {
    fside = p4est_iter_fside_array_index_int(&info->sides, i);

    int face = (int)fside->face;
    p4est_topidx_t treeid = fside->treeid;

    p4est_tree_t *tree = p4est_tree_array_index(info->p4est->trees, treeid);
    p4est_locidx_t offset = tree->quadrants_offset;

    BFAM_ASSERT(treeid >= 0 && treeid < info->p4est->connectivity->num_trees);

    if (!fside->is_hanging)
    {
      p4est_locidx_t qid = fside->is.full.quadid + offset;
      if (!fside->is.full.is_ghost)
      {
        BFAM_ASSERT(qid >= 0 && qid < info->p4est->local_num_quadrants);
        BFAM_LDEBUG("Glueid adding face (%jd %jd):(%jd %jd)", (intmax_t)treeid,
                    (intmax_t)face, (intmax_t)qid, (intmax_t)face);
        data->quad_to_glueid[P4EST_FACES * qid + face] =
            data->tree_to_glueid[P4EST_FACES * treeid + face];
      }
    }
    else
    {
      for (int h = 0; h < P4EST_HALF; ++h)
      {
        p4est_locidx_t qid = fside->is.hanging.quadid[h] + offset;
        if (!fside->is.hanging.is_ghost[h])
        {
          BFAM_ASSERT(qid >= 0 && qid < info->p4est->local_num_quadrants);
          BFAM_LDEBUG("Glueid adding face (%jd %jd):(%jd %jd)",
                      (intmax_t)treeid, (intmax_t)face, (intmax_t)qid,
                      (intmax_t)face);
          data->quad_to_glueid[P4EST_FACES * qid + face] =
              data->tree_to_glueid[P4EST_FACES * treeid + face];
        }
      }
    }
  }

  return;
}

void bfam_domain_pxest_quad_to_glueid(p4est_t *pxest,
                                      const bfam_locidx_t *tree_to_glueid,
                                      bfam_locidx_t *quad_to_glueid)
{
#ifdef BFAM_DEBUG
  bfam_domain_pxest_tree_to_glueid_check(pxest->connectivity, tree_to_glueid);
#endif

  for (p4est_locidx_t i = 0; i < pxest->local_num_quadrants * P4EST_FACES; ++i)
    quad_to_glueid[i] = -1;

  bfam_domain_glueid_iter_data_t data = {tree_to_glueid, quad_to_glueid};
  p4est_iterate(pxest, NULL, &data, NULL, set_quad_to_glueid,
#if DIM == 3
                NULL,
#endif
                NULL);
}

static bfam_locidx_t bfam_domain_pxest_num_parallel_faces(p4est_mesh_t *mesh)
{
  const p4est_locidx_t K = mesh->local_num_quadrants;

  bfam_locidx_t numParallelFaces = 0;

  for (p4est_locidx_t k = 0; k < K; ++k)
  {
    for (int f = 0; f < P4EST_FACES; ++f)
    {
      const p4est_locidx_t ck = mesh->quad_to_quad[P4EST_FACES * k + f];
      const int cf = mesh->quad_to_face[P4EST_FACES * k + f];

      if (cf >= 0)
      {
        /*
         * neighbor is same or double size
         */
        if (ck >= mesh->local_num_quadrants)
          ++numParallelFaces;
      }
      else
      {
        /*
         * two neighbors half size
         */
        p4est_locidx_t *cks;
        cks = sc_array_index(mesh->quad_to_half, ck);

        for (int h = 0; h < P4EST_HALF; ++h)
          if (cks[h] >= mesh->local_num_quadrants)
            ++numParallelFaces;
      }
    }
  }

  BFAM_LDEBUG("Counted %jd parallel faces.", (intmax_t)numParallelFaces);

  return numParallelFaces;
}

/** Build the parallel face mapping array.
 *
 * This array is used to determine the order in which data is sent
 * and received around the forest.  The order can be obtained by sorting
 * the mapping using \c bfam_subdomain_face_send_cmp
 * and \c bfam_domain_pxest_face_recv_cmp comparison
 * functions.
 *
 * \param [in]  mesh             p4est mesh to build the mapping for.
 * \param [in]  glueID           glue id of each face in the mesh
 * \param [in]  numParallelFaces the number of parallel faces in the
 *                               p4est mesh.
 * \param [out] mapping          the mapping array that will be filled.
 *
 */
static void bfam_domain_pxest_parallel_face_mapping(
    p4est_ghost_t *ghost, p4est_mesh_t *mesh, bfam_locidx_t *glueID,
    bfam_locidx_t numParallelFaces, bfam_subdomain_face_map_entry_t *mapping)
{
  const p4est_locidx_t K = mesh->local_num_quadrants;

  for (p4est_locidx_t k = 0, sk = 0; k < K; ++k)
  {
    for (int f = 0; f < P4EST_FACES; ++f)
    {
      const p4est_locidx_t ck = mesh->quad_to_quad[P4EST_FACES * k + f];
      const int cf = mesh->quad_to_face[P4EST_FACES * k + f];

      const bfam_locidx_t glueid = (glueID) ? glueID[P4EST_FACES * k + f] : -1;

      if (cf >= 0)
      {
        /*
         * neighbor is same or double size
         */
        if (ck >= mesh->local_num_quadrants)
        {
          p4est_locidx_t ghostid = ck - mesh->local_num_quadrants;
          p4est_quadrant_t *ghostquad =
              p4est_quadrant_array_index(&ghost->ghosts, (size_t)ghostid);
          p4est_locidx_t ghostp = mesh->ghost_to_proc[ghostid];
          p4est_locidx_t ghostk = ghostquad->p.piggy3.local_num;
          int ghostf = cf;
          int ghosth = 0;

          if (ghostf >= BFAM_PXEST_RELATIVE_ORIENTATIONS * P4EST_FACES)
          {
            ghostf -= BFAM_PXEST_RELATIVE_ORIENTATIONS * P4EST_FACES;

            ghosth =
                ghostf / (BFAM_PXEST_RELATIVE_ORIENTATIONS * P4EST_FACES) + 1;
            ghostf = ghostf % (BFAM_PXEST_RELATIVE_ORIENTATIONS * P4EST_FACES);

            BFAM_ASSERT(ghosth >= 0 && ghosth <= P4EST_HALF);
          }

          int o = ghostf / P4EST_FACES;

          ghostf = ghostf % P4EST_FACES;

          const int bo = BFAM_PXEST_ORIENTATION(f, ghostf, o);
          if (ghosth > 0)
            ghosth = BFAM_PXEST_BOHTONH(bo, ghosth - 1) + 1;

          mapping[sk].np = ghostp;
          mapping[sk].ns = 0;
          mapping[sk].nk = ghostk;
          mapping[sk].nf = (int8_t)ghostf;
          mapping[sk].nh = (int8_t)ghosth;
          mapping[sk].gi = ghostid;
          mapping[sk].i = sk;
          mapping[sk].s = 0;
          mapping[sk].k = k;
          mapping[sk].f = (int8_t)f;
          mapping[sk].h = (int8_t)0;
          mapping[sk].o = (int8_t)bo;
          mapping[sk].id = glueid;
          ++sk;
        }
      }
      else
      {
        /*
         * two neighbors half size
         */
        p4est_locidx_t *cks;
        cks = sc_array_index(mesh->quad_to_half, ck);

        for (int8_t h = 0; h < P4EST_HALF; ++h)
        {
          if (cks[h] >= mesh->local_num_quadrants)
          {
            p4est_locidx_t ghostid = cks[h] - mesh->local_num_quadrants;
            p4est_quadrant_t *ghostquad =
                p4est_quadrant_array_index(&ghost->ghosts, (size_t)ghostid);
            p4est_locidx_t ghostp = mesh->ghost_to_proc[ghostid];
            p4est_locidx_t ghostk = ghostquad->p.piggy3.local_num;
            int ghostf =
                ((BFAM_PXEST_RELATIVE_ORIENTATIONS * P4EST_FACES) + cf) %
                P4EST_FACES;
            int ghosth = 0;
            int o = ((BFAM_PXEST_RELATIVE_ORIENTATIONS * P4EST_FACES) + cf) /
                    P4EST_FACES;

            const int bo = BFAM_PXEST_ORIENTATION(f, ghostf, o);

            mapping[sk].np = ghostp;
            mapping[sk].ns = 0;
            mapping[sk].nk = ghostk;
            mapping[sk].nf = (int8_t)ghostf;
            mapping[sk].nh = (int8_t)ghosth;
            mapping[sk].gi = ghostid;
            mapping[sk].i = sk;
            mapping[sk].s = 0;
            mapping[sk].k = k;
            mapping[sk].f = (int8_t)f;
            mapping[sk].h = (int8_t)(BFAM_PXEST_BOHTONH_INV(bo, h) + 1);
            mapping[sk].o = (int8_t)bo;
            mapping[sk].id = glueid;
            ++sk;
          }
        }
      }
    }
  }

  qsort(mapping, numParallelFaces, sizeof(bfam_subdomain_face_map_entry_t),
        bfam_subdomain_face_send_cmp);

#ifdef BFAM_DEBUG
  {
    BFAM_LDEBUG("mapping: %5s %5s %5s %5s %5s %5s %5s %5s %5s %5s %5s", "np",
                "nk", "nf", "nh", "gi", "i", "k", "f", "h", "o", "id");
    for (bfam_locidx_t i = 0; i < numParallelFaces; ++i)
    {
      BFAM_LDEBUG(
          "mapping: %5jd %5jd %5jd %5jd %5jd %5jd %5jd %5jd %5jd %5jd %5jd",
          (intmax_t)mapping[i].np, (intmax_t)mapping[i].nk,
          (intmax_t)mapping[i].nf, (intmax_t)mapping[i].nh,
          (intmax_t)mapping[i].gi, (intmax_t)mapping[i].i,
          (intmax_t)mapping[i].k, (intmax_t)mapping[i].f,
          (intmax_t)mapping[i].h, (intmax_t)mapping[i].o,
          (intmax_t)mapping[i].id);
    }
  }
#endif
}

/** Count the number of processor neighbors.
 *
 * \param [in] numParallelFaces number of parallel faces in the mapping
 * \param [in] mapping          parallel face mapping array
 *
 * \returns the number of processor neighbors.
 */
static bfam_locidx_t bfam_domain_pxest_parallel_face_num_neighbors(
    bfam_locidx_t numParallelFaces, bfam_subdomain_face_map_entry_t *mapping)
{
  qsort(mapping, numParallelFaces, sizeof(bfam_subdomain_face_map_entry_t),
        bfam_subdomain_face_send_cmp);

  bfam_locidx_t numNeighbors = 0;

  if (numParallelFaces != 0)
  {
    numNeighbors = 1;
    for (bfam_locidx_t i = 1; i < numParallelFaces; ++i)
      if (mapping[i - 1].np != mapping[i].np)
        ++numNeighbors;
  }

  return numNeighbors;
}

/** Determine the neighbor ranks and number of faces per neighbor.
 *
 * \param [in]  numParallelFaces number of parallel faces in the mapping
 * \param [in]  mapping          parallel face mapping array
 * \param [in]  numNeighbors     the number of processor neighbors
 * \param [out] numNeighborFaces number of faces per neighbor
 * \param [out] neighborRank     rank of the neighbor
 */
static void bfam_domain_pxest_num_neighbor_faces(
    bfam_locidx_t numParallelFaces, bfam_subdomain_face_map_entry_t *mapping,
    bfam_locidx_t numNeighbors, bfam_locidx_t *numNeighborFaces,
    bfam_locidx_t *neighborRank)
{
  qsort(mapping, numParallelFaces, sizeof(bfam_subdomain_face_map_entry_t),
        bfam_subdomain_face_send_cmp);

  if (numParallelFaces != 0)
  {
    BFAM_ASSERT(numNeighbors > 0);

    numNeighborFaces[0] = 1;
    neighborRank[0] = mapping[0].np;

    for (bfam_locidx_t i = 1, sk = 0; i < numParallelFaces; ++i)
    {
      if (mapping[i - 1].np == mapping[i].np)
      {
        ++numNeighborFaces[sk];
      }
      else
      {
        ++sk;
        numNeighborFaces[sk] = 1;
        neighborRank[sk] = mapping[i].np;
      }
    }
  }

#ifdef BFAM_DEBUG
  {
    /*
     * Sanity check on the number of parallel faces per processor
     */
    bfam_locidx_t faces = 0;
    for (bfam_locidx_t n = 0; n < numNeighbors; ++n)
      faces += numNeighborFaces[n];

    BFAM_LDEBUG(" XXX  NumNeighbors %jd", (intmax_t)numNeighbors);
    for (bfam_locidx_t n = 0; n < numNeighbors; ++n)
      BFAM_LDEBUG(" XXX     neighborFaces[%jd] = %jd", (intmax_t)n,
                  (intmax_t)numNeighborFaces[n]);

    for (bfam_locidx_t n = 0; n < numNeighbors; ++n)
      BFAM_LDEBUG(" XXX     neighborRanks[%jd] = %jd", (intmax_t)n,
                  (intmax_t)neighborRank[n]);

    BFAM_LDEBUG(" XXX  faces: %jd =?= %jd", (intmax_t)faces,
                (intmax_t)numParallelFaces);

    BFAM_ASSERT(faces == numParallelFaces);
  }
#endif
}

/** Send and receive \a entries values per parallel face.
 *
 * \param [in]  numParallelFaces number of parallel faces in the mapping
 * \param [in]  mapping          parallel face mapping array
 * \param [in]  numNeighbors     the number of processor neighbors
 * \param [in]  numNeighborFaces number of faces per neighbor
 * \param [in]  neighborRank     rank of the neighbor
 * \param [in]  entries          number of entries to send per face
 * \param [in]  sendBuffer       values to send
 * \param [out] recvBuffer       values received
 *
 */
static void bfam_domain_pxest_parallel_face_send_recv(
    MPI_Comm comm, bfam_locidx_t numParallelFaces,
    bfam_subdomain_face_map_entry_t *mapping, bfam_locidx_t numNeighbors,
    bfam_locidx_t *numNeighborFaces, bfam_locidx_t *neighborRank,
    bfam_locidx_t entries, bfam_locidx_t *sendBuffer, bfam_locidx_t *recvBuffer)
{
  const int tag = 666;

  MPI_Request *sendRequest =
      bfam_malloc(2 * numNeighbors * sizeof(MPI_Request));
  MPI_Request *recvRequest = sendRequest + numNeighbors;

  MPI_Status *status = bfam_malloc(2 * numNeighbors * sizeof(MPI_Status));

  for (bfam_locidx_t n = 0; n < numNeighbors; ++n)
  {
    sendRequest[n] = MPI_REQUEST_NULL;
    recvRequest[n] = MPI_REQUEST_NULL;
  }

  /*
   * Post receives
   */
  for (bfam_locidx_t n = 0, sk = 0; n < numNeighbors; ++n)
  {
    const bfam_locidx_t count = entries * numNeighborFaces[n];
    BFAM_MPI_CHECK(MPI_Irecv(&recvBuffer[sk], count, BFAM_LOCIDX_MPI,
                             neighborRank[n], tag, comm, &recvRequest[n]));
    sk += count;
  }

  /*
   * Post sends
   */
  for (bfam_locidx_t n = 0, sk = 0; n < numNeighbors; ++n)
  {
    const bfam_locidx_t count = entries * numNeighborFaces[n];
    BFAM_MPI_CHECK(MPI_Isend(&sendBuffer[sk], count, BFAM_LOCIDX_MPI,
                             neighborRank[n], tag, comm, &sendRequest[n]));
    sk += count;
  }

  /*
   * Wait for the communication
   */
  BFAM_MPI_CHECK(MPI_Waitall(2 * numNeighbors, sendRequest, status));

  bfam_free(status);
  bfam_free(sendRequest);
}

#ifdef BFAM_DEBUG
/** Debug function to test the mapping.
 *
 * \param [in] numParallelFaces number of parallel faces in the mapping
 * \param [in] mapping          parallel face mapping array
 * \param [in] numNeighbors     the number of processor neighbors
 * \param [in] numNeighborFaces number of faces per neighbor
 * \param [in] neighborRank     rank of the neighbor
 *
 * \note This function aborts if it found an error in the mapping
 *
 */
static void bfam_domain_pxest_parallel_face_mapping_check(
    MPI_Comm comm, bfam_locidx_t numParallelFaces,
    bfam_subdomain_face_map_entry_t *mapping, bfam_locidx_t numNeighbors,
    bfam_locidx_t *numNeighborFaces, bfam_locidx_t *neighborRank)
{
  const int entries = 9;

  bfam_locidx_t *recvBuffer =
      bfam_malloc_aligned(numParallelFaces * entries * sizeof(bfam_locidx_t));
  bfam_locidx_t *sendBuffer =
      bfam_malloc_aligned(numParallelFaces * entries * sizeof(bfam_locidx_t));

  /*
   * Sort the mapping in send order
   */
  qsort(mapping, numParallelFaces, sizeof(bfam_subdomain_face_map_entry_t),
        bfam_subdomain_face_send_cmp);

  /*
   * Fill Send buffers
   */
  for (bfam_locidx_t i = 0; i < numParallelFaces; ++i)
  {
    sendBuffer[entries * i + 0] = mapping[i].ns;
    sendBuffer[entries * i + 1] = mapping[i].nk;
    sendBuffer[entries * i + 2] = mapping[i].nf;
    sendBuffer[entries * i + 3] = mapping[i].nh;

    sendBuffer[entries * i + 4] = mapping[i].s;
    sendBuffer[entries * i + 5] = mapping[i].k;
    sendBuffer[entries * i + 6] = mapping[i].f;
    sendBuffer[entries * i + 7] = mapping[i].h;

    sendBuffer[entries * i + 8] = mapping[i].id;
  }

  bfam_domain_pxest_parallel_face_send_recv(
      comm, numParallelFaces, mapping, numNeighbors, numNeighborFaces,
      neighborRank, entries, sendBuffer, recvBuffer);

  /*
   * Sort the mapping in recv order
   */
  qsort(mapping, numParallelFaces, sizeof(bfam_subdomain_face_map_entry_t),
        bfam_subdomain_face_recv_cmp);

  /*
   * Check the receive buffers
   */
  for (bfam_locidx_t i = 0; i < numParallelFaces; ++i)
  {
    BFAM_ASSERT(recvBuffer[entries * i + 0] == mapping[i].s);
    BFAM_ASSERT(recvBuffer[entries * i + 1] == mapping[i].k);
    BFAM_ASSERT(recvBuffer[entries * i + 2] == mapping[i].f);
    BFAM_ASSERT(recvBuffer[entries * i + 3] == mapping[i].h);

    BFAM_ASSERT(recvBuffer[entries * i + 4] == mapping[i].ns);
    BFAM_ASSERT(recvBuffer[entries * i + 5] == mapping[i].nk);
    BFAM_ASSERT(recvBuffer[entries * i + 6] == mapping[i].nf);
    BFAM_ASSERT(recvBuffer[entries * i + 7] == mapping[i].nh);

    BFAM_ASSERT(recvBuffer[entries * i + 8] == mapping[i].id);
  }

  bfam_free_aligned(sendBuffer);
  bfam_free_aligned(recvBuffer);
}
#endif

/** Get neighbor's subdomain information.
 *
 * \param [in]     comm             MPI Communicator
 * \param [in]     numParallelFaces number of parallel faces in the mapping
 * \param [in,out] mapping          parallel face mapping array (subdomain
 *                                  information is filled in the mapping).
 * \param [in]     numNeighbors     the number of processor neighbors
 * \param [in]     numNeighborFaces number of faces per neighbor
 * \param [in]     neighborRank     rank of the neighbor
 * \param [in]     subdomainID      array of subdomain ids for each local
 *                                  element
 * \param [in]     N                array of orders for each subdomain
 * \param [out]    ghostSubdomainID array of subdomain ids for each ghost
 *                                  element in the mesh
 * \param [out]     ghostN          array of orders for each ghost element in
 *                                  the mesh
 */
static void bfam_domain_pxest_fill_ghost_subdomain_ids(
    MPI_Comm comm, bfam_locidx_t numParallelFaces,
    bfam_subdomain_face_map_entry_t *mapping, bfam_locidx_t numNeighbors,
    bfam_locidx_t *numNeighborFaces, bfam_locidx_t *neighborRank,
    bfam_locidx_t *subdomainID, int *N, bfam_locidx_t *ghostSubdomainID,
    bfam_locidx_t *ghostN)
{
  const int entries = 2;

  bfam_locidx_t *recvBuffer =
      bfam_malloc_aligned(numParallelFaces * entries * sizeof(bfam_locidx_t));
  bfam_locidx_t *sendBuffer =
      bfam_malloc_aligned(numParallelFaces * entries * sizeof(bfam_locidx_t));

  /*
   * Sort the mapping in send order
   */
  qsort(mapping, numParallelFaces, sizeof(bfam_subdomain_face_map_entry_t),
        bfam_subdomain_face_send_cmp);

  /*
   * Fill Send buffers
   */
  for (bfam_locidx_t i = 0; i < numParallelFaces; ++i)
  {
    bfam_locidx_t s = subdomainID[mapping[i].k];
    sendBuffer[entries * i + 0] = s;
    sendBuffer[entries * i + 1] = N[s];
  }

  bfam_domain_pxest_parallel_face_send_recv(
      comm, numParallelFaces, mapping, numNeighbors, numNeighborFaces,
      neighborRank, entries, sendBuffer, recvBuffer);

  /*
   * Sort the mapping in recv order
   */
  qsort(mapping, numParallelFaces, sizeof(bfam_subdomain_face_map_entry_t),
        bfam_subdomain_face_recv_cmp);

  /*
   * Fill mapping with subdomain id
   */
  for (bfam_locidx_t i = 0; i < numParallelFaces; ++i)
    mapping[i].s = subdomainID[mapping[i].k];

  /*
   * Fill the ghost information
   */
  for (bfam_locidx_t i = 0; i < numParallelFaces; ++i)
  {
    mapping[i].ns = recvBuffer[entries * i + 0];
    ghostSubdomainID[mapping[i].gi] = recvBuffer[entries * i + 0];
    ghostN[mapping[i].gi] = recvBuffer[entries * i + 1];
  }

  bfam_free_aligned(sendBuffer);
  bfam_free_aligned(recvBuffer);
}

/** Count the number of local inter-subdomain faces.
 *
 * Note that faces on nonconforming interfaces are always counted even if
 * they connect to the same subdomain.  This is because we will use a
 * glue grid to handle the nonconforming interfaces.
 *
 * \param [in]  mesh          p4est mesh of the elements
 * \param [in]  subdomainID   subdomain id of each element in the mesh
 * \param [in]  glueID        glue id of each face in the mesh
 *
 * \return number of total local inter-subdomain faces.
 *
 */
static bfam_locidx_t bfam_domain_pxest_num_inter_subdomain_faces(
    p4est_mesh_t *mesh, bfam_locidx_t *subdomainID, bfam_locidx_t *glueID)
{
  const p4est_locidx_t K = mesh->local_num_quadrants;

  bfam_locidx_t numInterSubdomainFaces = 0;

  for (p4est_locidx_t k = 0; k < K; ++k)
  {
    const bfam_locidx_t idk = subdomainID[k];

    for (int f = 0; f < P4EST_FACES; ++f)
    {
      const p4est_locidx_t ck = mesh->quad_to_quad[P4EST_FACES * k + f];
      const int cf = mesh->quad_to_face[P4EST_FACES * k + f];

      const bfam_locidx_t glueid = (glueID) ? glueID[P4EST_FACES * k + f] : -1;

      if (cf >= 0)
      {
        /*
         * Neighbor is same or double size
         */
        if (ck < mesh->local_num_quadrants)
        {
          /*
           * Neighbor is on the same processor
           */
          const bfam_locidx_t idnk = subdomainID[ck];

          int hanging;

          if (cf >= BFAM_PXEST_RELATIVE_ORIENTATIONS * P4EST_FACES)
            hanging = 1;
          else
            hanging = 0;

          /*
           * Count intra and inter subdomain faces.
           *
           * Only count same subdomain to same subdomain if it is a hanging
           * face.
           */
          if ((idnk == idk && hanging) || idnk != idk ||
              (glueid >= 0 && (ck != k || cf != f)))
            ++numInterSubdomainFaces;
        }
      }
      else
      {
        p4est_locidx_t *cks;
        cks = sc_array_index(mesh->quad_to_half, ck);
        for (int8_t h = 0; h < P4EST_HALF; ++h)
          if (cks[h] < mesh->local_num_quadrants)
            ++numInterSubdomainFaces;
      }
    }
  }

  return numInterSubdomainFaces;
}

/** Build the inter subdomain face mapping array.
 *
 * This array is used to determine the order in which data is sent
 * and received around the subdomains.  The order can be obtained by sorting
 * the mapping using \c bfam_subdomain_face_send_cmp
 * and \c bfam_subdomain_face_recv_cmp comparison
 * functions.
 *
 * \param [in]  rank                   local MPI rank.
 * \param [in]  mesh                   pxest mesh to build the mapping for.
 * \param [in]  subdomainID            subdomain id of each element in the mesh.
 * \param [in]  glueID                 glue id of each face in the mesh
 * \param [in]  numInterSubdomainFaces the number of inter subdomain faces in
 *                                     the pxest mesh.
 * \param [out] mapping                the mapping array that will be filled.
 *
 */
static void bfam_domain_pxest_inter_subdomain_face_mapping(
    bfam_locidx_t rank, p4est_mesh_t *mesh, bfam_locidx_t *subdomainID,
    bfam_locidx_t *glueID, bfam_locidx_t numInterSubdomainFaces,
    bfam_subdomain_face_map_entry_t *mapping)
{
  const p4est_locidx_t K = mesh->local_num_quadrants;

  for (p4est_locidx_t k = 0, sk = 0; k < K; ++k)
  {
    const bfam_locidx_t idk = subdomainID[k];

    for (int f = 0; f < P4EST_FACES; ++f)
    {
      const p4est_locidx_t ck = mesh->quad_to_quad[P4EST_FACES * k + f];
      const int cf = mesh->quad_to_face[P4EST_FACES * k + f];

      const bfam_locidx_t glueid = (glueID) ? glueID[P4EST_FACES * k + f] : -1;

      if (cf >= 0)
      {
        /*
         * Neighbor is same or double size
         */
        if (ck < mesh->local_num_quadrants)
        {
          /*
           * Neighbor is on the same processor
           */
          const bfam_locidx_t idnk = subdomainID[ck];

          int hanging;

          if (cf >= BFAM_PXEST_RELATIVE_ORIENTATIONS * P4EST_FACES)
            hanging = 1;
          else
            hanging = 0;

          int nf = cf;
          int nh = 0;

          if (nf >= BFAM_PXEST_RELATIVE_ORIENTATIONS * P4EST_FACES)
          {
            nf -= BFAM_PXEST_RELATIVE_ORIENTATIONS * P4EST_FACES;

            nh = nf / (BFAM_PXEST_RELATIVE_ORIENTATIONS * P4EST_FACES) + 1;
            nf = nf % (BFAM_PXEST_RELATIVE_ORIENTATIONS * P4EST_FACES);
          }

          int o = nf / P4EST_FACES;

          nf = nf % P4EST_FACES;
          const int bo = BFAM_PXEST_ORIENTATION(f, nf, o);

          if (nh > 0)
          {
            nh = BFAM_PXEST_BOHTONH(bo, nh - 1) + 1;
          }

          /*
           * Count intra and inter subdomain faces.
           *
           * Only count same subdomain to same subdomain if it is a hanging
           * face.
           */
          if ((idnk == idk && hanging) || idnk != idk ||
              (glueid >= 0 && (ck != k || cf != f)))
          {
            BFAM_ASSERT(sk < numInterSubdomainFaces);

            mapping[sk].np = rank;

            mapping[sk].ns = idnk;
            mapping[sk].nk = ck;
            mapping[sk].nf = (int8_t)nf;
            mapping[sk].nh = (int8_t)nh;

            mapping[sk].s = idk;
            mapping[sk].k = k;
            mapping[sk].f = (int8_t)f;
            mapping[sk].h = (int8_t)0;
            mapping[sk].o = (int8_t)bo;

            mapping[sk].i = -1;
            mapping[sk].gi = -1;

            mapping[sk].id = glueid;
            ++sk;
          }
        }
      }
      else
      {
        p4est_locidx_t *cks;
        cks = sc_array_index(mesh->quad_to_half, ck);
        for (int8_t h = 0; h < P4EST_HALF; ++h)
        {
          if (cks[h] < mesh->local_num_quadrants)
          {
            BFAM_ASSERT(sk < numInterSubdomainFaces);

            const bfam_locidx_t idnk = subdomainID[cks[h]];

            int nf = ((BFAM_PXEST_RELATIVE_ORIENTATIONS * P4EST_FACES) + cf) %
                     P4EST_FACES;
            int nh = 0;
            int o = ((BFAM_PXEST_RELATIVE_ORIENTATIONS * P4EST_FACES) + cf) /
                    P4EST_FACES;

            const int bo = BFAM_PXEST_ORIENTATION(f, nf, o);

            mapping[sk].np = rank;

            mapping[sk].ns = idnk;
            mapping[sk].nk = cks[h];
            mapping[sk].nf = (int8_t)nf;
            mapping[sk].nh = (int8_t)nh;

            mapping[sk].s = idk;
            mapping[sk].k = k;
            mapping[sk].f = (int8_t)f;
            mapping[sk].h = (int8_t)(BFAM_PXEST_BOHTONH_INV(bo, h) + 1);
            mapping[sk].o = (int8_t)bo;

            mapping[sk].i = -1;
            mapping[sk].gi = -1;

            mapping[sk].id = glueid;

            ++sk;
          }
        }
      }
    }
  }
}

/** Count the number of boundary faces.
 *
 * \param [in]  mesh pxest mesh of the elements
 *
 * \return number of boundary faces.
 *
 */
static bfam_locidx_t bfam_domain_pxest_num_boundary_faces(p4est_mesh_t *mesh)
{
  const p4est_locidx_t K = mesh->local_num_quadrants;

  bfam_locidx_t numBoundaryFaces = 0;

  for (p4est_locidx_t k = 0; k < K; ++k)
  {
    for (int f = 0; f < P4EST_FACES; ++f)
    {
      const p4est_locidx_t ck = mesh->quad_to_quad[P4EST_FACES * k + f];
      const int cf = mesh->quad_to_face[P4EST_FACES * k + f];

      if (k == ck && f == cf)
        ++numBoundaryFaces;
    }
  }

  return numBoundaryFaces;
}

/** Build the boundary subdomain face mapping array.
 *
 * \param [in]  mesh             p4est mesh to build the mapping for.
 * \param [in]  subdomainID      subdomain id of each element in the mesh.
 * \param [in]  glueID           glue id of each face in the mesh
 * \param [in]  numBoundaryFaces the number of boundary subdomain faces in
 *                               the p4est mesh.
 * \param [out] mapping          the mapping array that will be filled.
 *
 */
static void bfam_domain_pxest_boundary_subdomain_face_mapping(
    p4est_mesh_t *mesh, bfam_locidx_t *subdomainID, bfam_locidx_t *glueID,
    bfam_locidx_t numBoundaryFaces, bfam_subdomain_face_map_entry_t *mapping)
{
  const p4est_locidx_t K = mesh->local_num_quadrants;

  for (p4est_locidx_t k = 0, sk = 0; k < K; ++k)
  {
    const bfam_locidx_t idk = subdomainID[k];

    for (int f = 0; f < P4EST_FACES; ++f)
    {
      const p4est_locidx_t ck = mesh->quad_to_quad[P4EST_FACES * k + f];
      const int cf = mesh->quad_to_face[P4EST_FACES * k + f];

      const bfam_locidx_t glueid = (glueID) ? glueID[P4EST_FACES * k + f] : -1;

      if (k == ck && f == cf)
      {
        BFAM_ASSERT(sk < numBoundaryFaces);

        mapping[sk].np = -1;

        mapping[sk].ns = idk;
        mapping[sk].nk = k;
        mapping[sk].nf = (int8_t)f;
        mapping[sk].nh = (int8_t)0;

        mapping[sk].s = idk;
        mapping[sk].k = k;
        mapping[sk].f = (int8_t)f;
        mapping[sk].h = (int8_t)0;
        mapping[sk].o = (int8_t)0;

        mapping[sk].i = -1;
        mapping[sk].gi = -1;

        mapping[sk].id = glueid;
        ++sk;
      }
    }
  }
}

/** Mark quadrants in pxest for adaptation.
 *
 * This marks all of the quadrants for refinement and coarsening.
 *
 * \param [in,out] domain coming in the domain's pxest needs to be in sync with
 *                        the subdomains; and returning the pxest user data will
 *                        be updated with the new refinement flags based on
 *                        info in the subdomain (any currently set refinement
 *                        flags will be zeroed.)
 */
void bfam_domain_pxest_mark_elements(bfam_domain_pxest_t *domain)
{
  p4est_t *pxest = domain->pxest;
  /*
   * Fill the quadrant user data refinement flags
   */
  for (p4est_topidx_t t = pxest->first_local_tree; t <= pxest->last_local_tree;
       ++t)
  {
    p4est_tree_t *tree = p4est_tree_array_index(pxest->trees, t);
    sc_array_t *quadrants = &tree->quadrants;
    size_t num_quads = quadrants->elem_count;

    for (size_t zz = 0; zz < num_quads; ++zz)
    {
      p4est_quadrant_t *quad = p4est_quadrant_array_index(quadrants, zz);
      bfam_pxest_user_data_t *ud = quad->p.user_data;

      /*
       * Here we make the assumption that all subdomains associated with a pxest
       * quadrant are of dgx type
       */
      bfam_subdomain_dgx_t *subdomain =
          (bfam_subdomain_dgx_t *)bfam_domain_get_subdomain_by_num(
              (bfam_domain_t *)domain, ud->subd_id);

      /* set the new order of the element */
      ud->N = subdomain->padapt[ud->elem_id];

      /* reset adaption flags */
      ud->flags = subdomain->hadapt[ud->elem_id];

      BFAM_VERBOSE(
          "Marking pxest t:%03jd q:%03jd <-- s:%03jd k:%03jd N:%02d h:%02x",
          (intmax_t)t, (intmax_t)zz, (intmax_t)ud->subd_id,
          (intmax_t)ud->elem_id, ud->N, ud->flags);
    }
  }
}

int bfam_domain_pxest_quadrant_coarsen(p4est_t *p4est,
                                       p4est_topidx_t which_tree,
                                       p4est_quadrant_t *quadrants[])
{
  bfam_pxest_user_data_t *ud = quadrants[0]->p.user_data;
  bfam_locidx_t root_id = ud->root_id;

  BFAM_VERBOSE("Coarsen Callback: root_id %jd coarsen flag %02x",
               (intmax_t)root_id, BFAM_FLAG_COARSEN);

  int N = ud->N;
  bfam_locidx_t subd_id = ud->subd_id;
  for (int k = 0; k < P4EST_CHILDREN; ++k)
  {
    ud = quadrants[k]->p.user_data;
    BFAM_VERBOSE(
        "Coarsen Callback: child[%d] s:%03jd k:%03jd h:%02x r:%jd N:%jd", k,
        (intmax_t)ud->subd_id, (intmax_t)ud->elem_id, ud->flags,
        (intmax_t)ud->root_id, (intmax_t)ud->N);

    /* only coarsen if we all want to coarsen */
    if (!(ud->flags & BFAM_FLAG_COARSEN))
      return 0;

    /* only coarsen if we have the same root ids */
    if (ud->root_id != root_id)
      return 0;

    if (N != ud->N)
      return 0;

    /* Make sure all elements are in the same subdomain to coarsen */
    if (subd_id != ud->subd_id)
      return 0;
  }

#if 0
  /* TODO: FIXME: Problem with this logic / or data */
  /* Only coarsen if the parent faces have the same glue id */
  for (unsigned int f = 0; f < P4EST_FACES; ++f)
  {
    bfam_locidx_t glue_id = ud->glue_id[(f % 2) << (f / 2)];
    for (unsigned int c = 0; c < P4EST_CHILDREN; ++c)
    {
      ud = quadrants[c]->p.user_data;
      if (((c >> (f / 2)) % 2 == f % 2) && (glue_id != ud->glue_id[f]))
        return 0;
    }
  }
#endif
  BFAM_VERBOSE("  Coarsen!");
  return 1;
}

int bfam_domain_pxest_quadrant_refine(p4est_t *p4est, p4est_topidx_t which_tree,
                                      p4est_quadrant_t *quadrant)
{
  bfam_pxest_user_data_t *ud = quadrant->p.user_data;

  BFAM_VERBOSE("Refinement Callback: s:%03jd k%03jd h:%02x r:%jd",
               (intmax_t)ud->subd_id, (intmax_t)ud->elem_id, ud->flags,
               (intmax_t)ud->root_id);

  if (!(ud->flags & BFAM_FLAG_REFINE))
    return 0;

  BFAM_VERBOSE("  Refine!");
  return 1;
}
void bfam_domain_pxest_quadrant_init(p4est_t *p4est, p4est_topidx_t which_tree,
                                     p4est_quadrant_t *quadrant)
{
  memset(quadrant->p.user_data, 0, sizeof(bfam_pxest_user_data_t));
}

void bfam_domain_pxest_quadrant_replace(p4est_t *p4est,
                                        p4est_topidx_t which_tree,
                                        int num_outgoing,
                                        p4est_quadrant_t *outgoing[],
                                        int num_incoming,
                                        p4est_quadrant_t *incoming[])
{
  BFAM_ASSERT(num_outgoing != 1 || num_incoming != 1);
  if (num_outgoing == 1)
  {
    /* Refining: copy data to all children */
    for (int c = 0; c < num_incoming; ++c)
    {
      bfam_pxest_user_data_t *in_ud = incoming[c]->p.user_data;
      memcpy(in_ud, outgoing[0]->p.user_data, sizeof(bfam_pxest_user_data_t));

      /* remove glue from internal faces */
      in_ud->glue_id[(c / 1 + 1) % 2 + 0] = -1;
      in_ud->glue_id[(c / 2 + 1) % 2 + 2] = -1;
#if DIM == 3
      in_ud->glue_id[(c / 4 + 1) % 2 + 4] = -1;
#endif
    }
  }
  else
  {
    /* Coarsening: copy data from the first child */
    BFAM_ASSERT(num_incoming == 1);
    bfam_pxest_user_data_t *in_ud = incoming[0]->p.user_data;
    bfam_pxest_user_data_t *out_ud;

    memcpy(in_ud, outgoing[0]->p.user_data, sizeof(bfam_pxest_user_data_t));
    for (int c = 0; c < num_outgoing; ++c)
      in_ud->N = BFAM_MAX(
          in_ud->N, ((bfam_pxest_user_data_t *)outgoing[c]->p.user_data)->N);

    for (unsigned int f = 0; f < P4EST_FACES; ++f)
    {
      /* grab a parent on the face and use their glue id */
      unsigned int parent = (f % 2) << (f / 2);
      out_ud = outgoing[parent]->p.user_data;
      in_ud->glue_id[f] = out_ud->glue_id[f];
    }
  }

  /* mark elements as new */
  for (int c = 0; c < num_incoming; ++c)
  {
    bfam_pxest_user_data_t *ud = incoming[c]->p.user_data;
    ud->flags |= BFAM_FLAG_ADAPTED;
  }
}

static int bfam_domain_pxest_select_N(uint8_t pflags, int N_old, int N_req)
{
  int N_new;
  if (N_req > 0)
  {
    if ((pflags & (BFAM_FLAG_COARSEN | BFAM_FLAG_REFINE)) ==
        (BFAM_FLAG_COARSEN | BFAM_FLAG_REFINE))
      N_new = N_req;
    else if (pflags & BFAM_FLAG_COARSEN)
      N_new = BFAM_MIN(N_old, N_req);
    else if (pflags & BFAM_FLAG_REFINE)
      N_new = BFAM_MAX(N_old, N_req);
    else
      N_new = N_old;
  }
  else
    N_new = N_old;

  return N_new;
}

void bfam_domain_pxest_compute_split(p4est_t *pxest, uint8_t pflags,
                                     bfam_locidx_t *num_subdomains,
                                     bfam_locidx_t **subdomain_id,
                                     bfam_locidx_t **roots, int **N,
                                     bfam_locidx_t **glue_id)
{
  char key[BFAM_BUFSIZ];
  bfam_dictionary_t rootN_to_sub;

  bfam_dictionary_init(&rootN_to_sub);

  bfam_locidx_t K = 0, k = 0;
  *num_subdomains = 0;

  /* find subdomains */
  for (p4est_topidx_t t = pxest->first_local_tree; t <= pxest->last_local_tree;
       ++t)
  {
    int retval;
    p4est_tree_t *tree = p4est_tree_array_index(pxest->trees, t);
    sc_array_t *quadrants = &tree->quadrants;
    size_t num_quads = quadrants->elem_count;

    for (size_t zz = 0; zz < num_quads; ++zz)
    {
      p4est_quadrant_t *quad = p4est_quadrant_array_index(quadrants, zz);
      bfam_pxest_user_data_t *ud = quad->p.user_data;

      /* Change order if we are refining or coarsening */
      const int N_new = bfam_domain_pxest_select_N(pflags, ud->Nold, ud->N);

      snprintf(key, BFAM_BUFSIZ, "%jd_%d", (intmax_t)ud->root_id, (int)N_new);

      retval =
          bfam_dictionary_insert_locidx(&rootN_to_sub, key, *num_subdomains);
      BFAM_ABORT_IF(retval == 0, "Can't insert '%s' into dictionary", key);

      /* If we have a new subdomain increment it */
      if (retval == 2)
        ++(*num_subdomains);

      ++K;
    }
  }

  BFAM_ASSERT((bfam_locidx_t)rootN_to_sub.num_entries == *num_subdomains);

  /* compute new split */
  *subdomain_id = bfam_malloc_aligned(K * sizeof(bfam_locidx_t));
  *roots = bfam_malloc_aligned(*num_subdomains * sizeof(bfam_locidx_t));
  *N = bfam_malloc_aligned(*num_subdomains * sizeof(bfam_locidx_t));
  *glue_id = bfam_malloc_aligned(P4EST_FACES * K * sizeof(bfam_locidx_t));

  k = 0;
  for (p4est_topidx_t t = pxest->first_local_tree; t <= pxest->last_local_tree;
       ++t)
  {
    int retval;
    p4est_tree_t *tree = p4est_tree_array_index(pxest->trees, t);
    sc_array_t *quadrants = &tree->quadrants;
    size_t num_quads = quadrants->elem_count;

    for (size_t zz = 0; zz < num_quads; ++zz, ++k)
    {
      p4est_quadrant_t *quad = p4est_quadrant_array_index(quadrants, zz);
      bfam_pxest_user_data_t *ud = quad->p.user_data;

      const int N_new = bfam_domain_pxest_select_N(pflags, ud->Nold, ud->N);

      snprintf(key, BFAM_BUFSIZ, "%jd_%d", (intmax_t)ud->root_id, (int)N_new);

      retval = bfam_dictionary_get_value_locidx(&rootN_to_sub, key,
                                                &(*subdomain_id)[k]);
      BFAM_ABORT_IF(retval == 0, "rootN key `%s` does not exist", key);
      (*N)[(*subdomain_id)[k]] = N_new;
      (*roots)[(*subdomain_id)[k]] = ud->root_id;
      for (int f = 0; f < P4EST_FACES; ++f)
        (*glue_id)[k * P4EST_FACES + f] = ud->glue_id[f];
    }
  }

  bfam_dictionary_clear(&rootN_to_sub);
}

typedef struct
{
  bfam_subdomain_dgx_t *subdomain_dst;
  bfam_domain_pxest_t *domain_src;
  bfam_domain_pxest_transfer_maps_t *maps;
  bfam_dictionary_t *N2N;
} bfam_subdomain_dgx_transfer_field_data_t;

static int quadrant_compare(const p4est_quadrant_t *quad_dst,
                            const p4est_quadrant_t *quad_src)
{
  if (p4est_quadrant_is_parent(quad_dst, quad_src))
  {
    /* h-coarsening */
    return -1;
  }
  else if (p4est_quadrant_is_parent(quad_src, quad_dst))
  {
    /* h-refining */
    return 1;
  }
  else if (p4est_quadrant_is_equal(quad_src, quad_dst))
  {
    /* h-same */
    return 0;
  }
  else
  {
    if (p4est_quadrant_overlaps(quad_dst, quad_src))
      BFAM_ABORT("Transfer aborted: more than one level of refinement "
                 "between quadrants");
    else
      BFAM_ABORT("Transfer aborted: Strange Things are Afoot at the Circle "
                 "K, we should never reach here");
    return 0;
  }
}

void bfam_domain_pxest_transfer_maps_init(
    bfam_domain_pxest_transfer_maps_t *maps, bfam_domain_pxest_t *domain_dst,
    bfam_domain_pxest_t *domain_src)
{
  p4est_t *pxest_dst = domain_dst->pxest;
  p4est_t *pxest_src = domain_src->pxest;

  p4est_topidx_t t;
  p4est_locidx_t k_src;
  p4est_locidx_t k_dst;
  p4est_locidx_t coarse_k_dst;

  BFAM_ABORT_IF(pxest_dst->first_local_tree != pxest_src->first_local_tree ||
                    pxest_dst->last_local_tree != pxest_src->last_local_tree,
                "Transfer aborted: Non-nested domains (Trees don't match)");

  /*
   * Count coarsened
   */
  size_t num_coarsened = 0;
  for (t = pxest_dst->first_local_tree; t <= pxest_dst->last_local_tree; ++t)
  {
    p4est_tree_t *tree_dst = p4est_tree_array_index(pxest_dst->trees, t);
    p4est_tree_t *tree_src = p4est_tree_array_index(pxest_src->trees, t);

    sc_array_t *quadrants_dst = &tree_dst->quadrants;
    sc_array_t *quadrants_src = &tree_src->quadrants;

    size_t z_src = 0;
    for (size_t z_dst = 0; z_dst < quadrants_dst->elem_count;)
    {
      BFAM_ASSERT(z_dst < quadrants_dst->elem_count);
      BFAM_ASSERT(z_src < quadrants_src->elem_count);

      p4est_quadrant_t *quad_dst =
          p4est_quadrant_array_index(quadrants_dst, z_dst);
      p4est_quadrant_t *quad_src =
          p4est_quadrant_array_index(quadrants_src, z_src);

      switch (quadrant_compare(quad_dst, quad_src))
      {
      case -1: /* h-coarsened */
        ++num_coarsened;
        z_src += P4EST_CHILDREN;
        z_dst += 1;
        break;
      case 1: /* h-refined */
        z_src += 1;
        z_dst += P4EST_CHILDREN;
        break;
      case 0: /* h-same */
        z_src += 1;
        z_dst += 1;
        break;
      default:
        BFAM_ABORT("Never should be reached");
      }
    }
  }

  BFAM_VERBOSE("Building transfer maps with %zd coarsened elements",
               num_coarsened);

  maps->num_dst = pxest_dst->local_num_quadrants;
  maps->dst_to_adapt_flags =
      bfam_malloc_aligned(maps->num_dst * sizeof(uint8_t));
  maps->dst_to_dst_chld_id =
      bfam_malloc_aligned(maps->num_dst * sizeof(int8_t));
  maps->dst_to_src_subd_id =
      bfam_malloc_aligned(maps->num_dst * sizeof(bfam_locidx_t));
  maps->dst_to_src_elem_id =
      bfam_malloc_aligned(maps->num_dst * sizeof(bfam_locidx_t));

  maps->num_coarse_dst = (bfam_locidx_t)num_coarsened;
  maps->coarse_dst_to_src_chld_id =
      bfam_malloc_aligned(P4EST_CHILDREN * num_coarsened * sizeof(int8_t));
  maps->coarse_dst_to_src_subd_id = bfam_malloc_aligned(
      P4EST_CHILDREN * num_coarsened * sizeof(bfam_locidx_t));
  maps->coarse_dst_to_src_elem_id = bfam_malloc_aligned(
      P4EST_CHILDREN * num_coarsened * sizeof(bfam_locidx_t));

  /*
   * Fill Maps
   */
  k_src = 0;
  k_dst = 0;
  coarse_k_dst = 0;
  for (t = pxest_dst->first_local_tree; t <= pxest_dst->last_local_tree; ++t)
  {
    p4est_tree_t *tree_dst = p4est_tree_array_index(pxest_dst->trees, t);
    p4est_tree_t *tree_src = p4est_tree_array_index(pxest_src->trees, t);

    sc_array_t *quadrants_dst = &tree_dst->quadrants;
    sc_array_t *quadrants_src = &tree_src->quadrants;

    size_t z_src = 0;
    for (size_t z_dst = 0; z_dst < quadrants_dst->elem_count;)
    {
      BFAM_ASSERT(z_dst < quadrants_dst->elem_count);
      BFAM_ASSERT(z_src < quadrants_src->elem_count);

      p4est_quadrant_t *quad_dst =
          p4est_quadrant_array_index(quadrants_dst, z_dst);
      p4est_quadrant_t *quad_src =
          p4est_quadrant_array_index(quadrants_src, z_src);

      BFAM_ASSERT(p4est_quadrant_overlaps(quad_dst, quad_src));

      bfam_pxest_user_data_t *ud_src = quad_src->p.user_data;

      switch (quadrant_compare(quad_dst, quad_src))
      {
      case -1: /* h-coarsened */
        maps->dst_to_adapt_flags[k_dst] = BFAM_FLAG_COARSEN;
        maps->dst_to_dst_chld_id[k_dst] =
            (int8_t)p4est_quadrant_child_id(quad_dst);

        /* Store index into coarse maps */
        maps->dst_to_src_subd_id[k_dst] = -1 - coarse_k_dst;
        maps->dst_to_src_elem_id[k_dst] = -1 - coarse_k_dst;

        for (int c = 0; c < P4EST_CHILDREN; ++c)
        {
          size_t c_id = P4EST_CHILDREN * coarse_k_dst + c;
          BFAM_ASSERT(c_id < P4EST_CHILDREN * num_coarsened);
          quad_src = p4est_quadrant_array_index(quadrants_src, z_src + c);
          ud_src = quad_src->p.user_data;

          BFAM_ASSERT(ud_src->subd_id >= 0);
          BFAM_ASSERT(ud_src->elem_id >= 0);

          maps->coarse_dst_to_src_chld_id[c_id] =
              (int8_t)p4est_quadrant_child_id(quad_src);
          maps->coarse_dst_to_src_subd_id[c_id] = ud_src->subd_id;
          maps->coarse_dst_to_src_elem_id[c_id] = ud_src->elem_id;
        }

        coarse_k_dst += 1;
        k_src += P4EST_CHILDREN;
        k_dst += 1;
        z_src += P4EST_CHILDREN;
        z_dst += 1;
        break;
      case 1: /* h-refined */
        BFAM_ASSERT(ud_src->subd_id >= 0);
        BFAM_ASSERT(ud_src->elem_id >= 0);

        for (int c = 0; c < P4EST_CHILDREN; ++c)
        {
          quad_dst = p4est_quadrant_array_index(quadrants_dst, z_dst + c);

          maps->dst_to_adapt_flags[k_dst + c] = BFAM_FLAG_REFINE;
          maps->dst_to_dst_chld_id[k_dst + c] =
              (int8_t)p4est_quadrant_child_id(quad_dst);
          maps->dst_to_src_subd_id[k_dst + c] = ud_src->subd_id;
          maps->dst_to_src_elem_id[k_dst + c] = ud_src->elem_id;
        }
        k_src += 1;
        k_dst += P4EST_CHILDREN;
        z_src += 1;
        z_dst += P4EST_CHILDREN;
        break;
      case 0: /* h-same */
        BFAM_ASSERT(ud_src->subd_id >= 0);
        BFAM_ASSERT(ud_src->elem_id >= 0);

        maps->dst_to_adapt_flags[k_dst] = 0;
        maps->dst_to_dst_chld_id[k_dst] =
            (int8_t)p4est_quadrant_child_id(quad_dst);
        maps->dst_to_src_subd_id[k_dst] = ud_src->subd_id;
        maps->dst_to_src_elem_id[k_dst] = ud_src->elem_id;

        k_src += 1;
        k_dst += 1;
        z_src += 1;
        z_dst += 1;
        break;
      default:
        BFAM_ABORT("Never should be reached");
      }
    }
  }
}

void bfam_domain_pxest_transfer_maps_free(
    bfam_domain_pxest_transfer_maps_t *maps)
{
  bfam_free_aligned(maps->dst_to_adapt_flags);
  bfam_free_aligned(maps->dst_to_dst_chld_id);
  bfam_free_aligned(maps->dst_to_src_subd_id);
  bfam_free_aligned(maps->dst_to_src_elem_id);
  bfam_free_aligned(maps->coarse_dst_to_src_chld_id);
  bfam_free_aligned(maps->coarse_dst_to_src_subd_id);
  bfam_free_aligned(maps->coarse_dst_to_src_elem_id);
}

// }}}

// {{{ subdomain dgx

#define BFAM_LOAD_FIELD_RESTRICT_ALIGNED(field, prefix, base, dictionary)      \
  bfam_real_t *restrict field;                                                 \
  {                                                                            \
    char bfam_load_field_name[BFAM_BUFSIZ];                                    \
    snprintf(bfam_load_field_name, BFAM_BUFSIZ, "%s%s", (prefix), (base));     \
    field = bfam_dictionary_get_value_ptr(dictionary, bfam_load_field_name);   \
    BFAM_ASSERT(field != NULL);                                                \
  }                                                                            \
  BFAM_ASSUME_ALIGNED(field, 32);
#define BFAM_LOAD_FIELD_ALIGNED(field, prefix, base, dictionary)               \
  bfam_real_t *field;                                                          \
  {                                                                            \
    char bfam_load_field_name[BFAM_BUFSIZ];                                    \
    snprintf(bfam_load_field_name, BFAM_BUFSIZ, "%s%s", (prefix), (base));     \
    field = bfam_dictionary_get_value_ptr(dictionary, bfam_load_field_name);   \
    BFAM_ASSERT(field != NULL);                                                \
  }                                                                            \
  BFAM_ASSUME_ALIGNED(field, 32);

static void init_interpolator(bfam_subdomain_dgx_interpolator_t *interp_a2b,
                              const int N_a, const int N_b)
{
  const bfam_locidx_t num_prj = 5;
  const bfam_locidx_t Np_a = N_a + 1;
  const bfam_locidx_t Np_b = N_b + 1;

  interp_a2b->N_src = N_a;
  interp_a2b->N_dst = N_b;
  interp_a2b->num_prj = num_prj;
  interp_a2b->prj = bfam_malloc(num_prj * sizeof(bfam_real_t *));

  /* Storage for the projection operators */
  for (bfam_locidx_t k = 0; k < num_prj; k++)
  {
    if (k == 0 && N_a == N_b)
      interp_a2b->prj[k] = NULL;
    else
      interp_a2b->prj[k] =
          bfam_malloc_aligned(Np_a * Np_b * sizeof(bfam_real_t));
  }

  interp_a2b->mass_prj = bfam_malloc(num_prj * sizeof(bfam_real_t *));

  interp_a2b->wi_mass_prj = bfam_malloc(num_prj * sizeof(bfam_real_t *));

  /* Storage for the projection operators */
  for (bfam_locidx_t k = 0; k < num_prj; k++)
    interp_a2b->mass_prj[k] =
        bfam_malloc_aligned(Np_a * Np_b * sizeof(bfam_real_t));
  for (bfam_locidx_t k = 0; k < num_prj; k++)
    interp_a2b->wi_mass_prj[k] =
        bfam_malloc_aligned(Np_a * Np_b * sizeof(bfam_real_t));
}

static void fill_grid_data(bfam_locidx_t N, bfam_long_real_t *lr,
                           bfam_long_real_t *lw, bfam_long_real_t *lV,
                           bfam_long_real_t *M)
{
  bfam_jacobi_gauss_lobatto_quadrature(0, 0, N, lr, lw);
  bfam_jacobi_p_vandermonde(0, 0, N, N + 1, lr, lV);
  bfam_jacobi_p_mass(0, 0, N, lV, M);
}

/* Build interpolation and projection between space a and g. Space a is the
 * lower order space and space g is the higher order space.
 */
static void fill_interp_proj_data(bfam_locidx_t N_a, bfam_locidx_t N_g,
                                  bfam_long_real_t *V_a, bfam_long_real_t *M_a,
                                  bfam_long_real_t *M_g, bfam_long_real_t *lr_g,
                                  bfam_long_real_t *I_a2g,
                                  bfam_long_real_t *P_g2a)
{
  BFAM_ASSERT(N_g >= N_a);
  bfam_locidx_t Np_a = N_a + 1;
  bfam_locidx_t Np_g = N_g + 1;

  bfam_jacobi_p_interpolation(0, 0, N_a, Np_g, lr_g, V_a, I_a2g);

  bfam_long_real_t MP_g2a[Np_g * Np_a];
  bfam_long_real_t Vt_MP_g2a[Np_g * Np_a];
  for (bfam_locidx_t n = 0; n < Np_g * Np_a; n++)
  {
    MP_g2a[n] = 0;
    Vt_MP_g2a[n] = 0;
    P_g2a[n] = 0;
  }

  /* MP_g2a = M_a * P_g2a = I_a2g^{T} * M_g */
  /* [Np_a X Np_g] = [Np_a X Np_g] [Np_g X Np_g] */
  bfam_util_mTmmult(Np_a, Np_g, Np_g, I_a2g, Np_g, M_g, Np_g, MP_g2a, Np_a);

  /* Vt_MP_g2a = V_a^{-1} * P_g2a = V_a^{T} MP_g2a */
  /* [Np_a X Np_g] = [Np_a X Np_a] [Np_a X Np_g] */
  bfam_util_mTmmult(Np_a, Np_g, Np_a, V_a, Np_a, MP_g2a, Np_a, Vt_MP_g2a, Np_a);

  /* P_g2a = V_a * V_a^{T} I_a2g^{T} * M_g */
  /* [Np_a X Np_g] = [Np_a X Np_a] [Np_a X Np_g] */
  bfam_util_mmmult(Np_a, Np_g, Np_a, V_a, Np_a, Vt_MP_g2a, Np_a, P_g2a, Np_a);

  /*
  if (N_g != N_a)
  {
    BFAM_INFO("%d -> %d", (int)N_g, (int)N_a);
    for (int i = 0; i < Np_a; i++)
      for (int j = 0; j < Np_g; j++)
        BFAM_INFO("P_g2a[%d][%d] = %+" BFAM_REAL_PRIe, i, j,
                  (bfam_real_t)P_g2a[i * Np_g + j]);
  }
  */
}

static void fill_hanging_data(bfam_locidx_t N, bfam_locidx_t Np,
                              bfam_long_real_t *lr_f, bfam_long_real_t *M,
                              bfam_long_real_t *lV, bfam_long_real_t *P_b2f,
                              bfam_long_real_t *P_t2f, bfam_long_real_t *I_f2b,
                              bfam_long_real_t *I_f2t)
{
  const bfam_long_real_t HALF = BFAM_LONG_REAL(0.5);

  /* Grid for the top and bottom */
  bfam_long_real_t lr_t[Np];
  bfam_long_real_t lr_b[Np];

  for (bfam_locidx_t k = 0; k < Np; k++)
  {
    lr_t[k] = HALF * (lr_f[k] + 1);
    lr_b[k] = HALF * (lr_f[k] - 1);
  }

  /* Mass matrix for the top and bottom */
  bfam_long_real_t Mh[Np * Np];
  for (bfam_locidx_t k = 0; k < Np * Np; k++)
    Mh[k] = HALF * M[k];

  fill_interp_proj_data(N, N, lV, M, Mh, lr_t, I_f2t, P_t2f);
  fill_interp_proj_data(N, N, lV, M, Mh, lr_b, I_f2b, P_b2f);
}

static void multiply_projections(const int N_b, const int N_a, const int N_g,
                                 bfam_long_real_t *P_g2b,
                                 bfam_long_real_t *P_g2g,
                                 bfam_long_real_t *P_a2g, bfam_real_t *P_a2b,
                                 bfam_long_real_t *M_b, bfam_real_t *MP_a2b,
                                 bfam_long_real_t *w_b, bfam_real_t *wiMP_a2b)
{
  const int Np_b = N_b + 1;
  const int Np_g = N_g + 1;
  const int Np_a = N_a + 1;

  BFAM_ASSERT(MP_a2b);

  /* In the case the target is NULL return */
  if (!P_a2b)
  {
    BFAM_ASSERT(Np_a == Np_b);
    for (bfam_locidx_t n = 0; n < Np_a * Np_b; n++)
      MP_a2b[n] = (bfam_real_t)(M_b[n]);
    for (bfam_locidx_t j = 0; j < Np_a; j++)
      for (bfam_locidx_t i = 0; i < Np_b; i++)
        wiMP_a2b[i + j * Np_b] = (bfam_real_t)(M_b[i + j * Np_b] / w_b[i]);
    return;
  }

  /* First set up the long storage for multiplication */

  bfam_long_real_t tmp[Np_b * Np_a];
  bfam_long_real_t *l_P = tmp;
  bfam_long_real_t l_MP[Np_b * Np_a];

  for (bfam_locidx_t n = 0; n < Np_a * Np_b; n++)
  {
    tmp[n] = 0;
    l_MP[n] = 0;
  }
  /* No glue to b (thus b is glue) */
  if (!P_g2b && P_g2g && P_a2g)
    bfam_util_mmmult(Np_b, Np_a, Np_g, P_g2g, Np_b, P_a2g, Np_g, l_P, Np_b);
  else if (P_g2b && P_g2g && !P_a2g)
    bfam_util_mmmult(Np_b, Np_a, Np_g, P_g2b, Np_b, P_g2g, Np_g, l_P, Np_b);
  else if (!P_g2b && !P_g2g && P_a2g)
    l_P = P_a2g;
  else if (!P_g2b && P_g2g && !P_a2g)
    l_P = P_g2g;
  else if (P_g2b && !P_g2g && !P_a2g)
    l_P = P_g2b;
  else
    BFAM_ABORT("Case of all NULL or all not NULL is not handled");

  /* MP_a2b = M_b * P_a2b */
  /* [Np_b X Np_a] = [Np_b X Np_b] [Np_b X Np_a]  */
  bfam_util_mmmult(Np_b, Np_a, Np_b, M_b, Np_b, l_P, Np_b, l_MP, Np_b);
  for (bfam_locidx_t n = 0; n < Np_a * Np_b; n++)
  {
    P_a2b[n] = (bfam_real_t)l_P[n];
    MP_a2b[n] = (bfam_real_t)(l_MP[n]);
  }
  for (bfam_locidx_t j = 0; j < Np_a; j++)
    for (bfam_locidx_t i = 0; i < Np_b; i++)
      wiMP_a2b[i + j * Np_b] = (bfam_real_t)(l_MP[i + j * Np_b] / w_b[i]);
}

static void create_interpolators(bfam_subdomain_dgx_interpolator_t *interp_a2b,
                                 bfam_subdomain_dgx_interpolator_t *interp_b2a,
                                 const int N_a, const int N_b)
{
  /* Figure out the glue space and the number of points */
  const int N_g = BFAM_MAX(N_a, N_b);
  const int Np_g = N_g + 1;
  const bfam_locidx_t Np_a = N_a + 1;
  const bfam_locidx_t Np_b = N_b + 1;

  /* initialize the interpolators */
  init_interpolator(interp_a2b, N_a, N_b);
  if (interp_b2a)
    init_interpolator(interp_b2a, N_b, N_a);

  /* Set up the reference grids for the side a */
  bfam_long_real_t lr_a[Np_a];
  bfam_long_real_t lw_a[Np_a];
  bfam_long_real_t V_a[Np_a * Np_a];
  bfam_long_real_t M_a[Np_a * Np_a];
  fill_grid_data(N_a, lr_a, lw_a, V_a, M_a);

  /* Set up the reference grids for the side b */
  bfam_long_real_t lr_b[Np_b];
  bfam_long_real_t lw_b[Np_b];
  bfam_long_real_t V_b[Np_b * Np_b];
  bfam_long_real_t M_b[Np_b * Np_b];
  fill_grid_data(N_b, lr_b, lw_b, V_b, M_b);

  /* Set up the reference grids for the glue space */
  bfam_long_real_t lr_g[Np_g];
  bfam_long_real_t lw_g[Np_g];
  bfam_long_real_t V_g[Np_g * Np_g];
  bfam_long_real_t M_g[Np_g * Np_g];
  fill_grid_data(N_g, lr_g, lw_g, V_g, M_g);

  /* Interpolate to the hanging faces */
  bfam_long_real_t *prj_g[5];
  prj_g[0] = NULL;
  for (bfam_locidx_t k = 1; k < 5; k++)
    prj_g[k] = bfam_malloc_aligned(sizeof(bfam_long_real_t) * Np_g * Np_g);

  fill_hanging_data(N_g, Np_g, lr_g, M_g, V_g, prj_g[1], prj_g[2], prj_g[3],
                    prj_g[4]);

  /* Interpolate to the intermediate space and projection back */
  bfam_long_real_t *I_a2g = NULL;
  bfam_long_real_t *P_g2a = NULL;
  bfam_long_real_t *I_b2g = NULL;
  bfam_long_real_t *P_g2b = NULL;
  if (N_a < N_g && N_b == N_g)
  {
    I_a2g = bfam_malloc_aligned(sizeof(bfam_long_real_t) * Np_g * Np_a);
    P_g2a = bfam_malloc_aligned(sizeof(bfam_long_real_t) * Np_g * Np_a);
    fill_interp_proj_data(N_a, N_g, V_a, M_a, M_g, lr_g, I_a2g, P_g2a);
  }
  if (N_a == N_g && N_b < N_g)
  {
    I_b2g = bfam_malloc_aligned(sizeof(bfam_long_real_t) * Np_g * Np_b);
    P_g2b = bfam_malloc_aligned(sizeof(bfam_long_real_t) * Np_g * Np_b);
    fill_interp_proj_data(N_b, N_g, V_b, M_b, M_g, lr_g, I_b2g, P_g2b);
  }

  for (bfam_locidx_t k = 0; k < 5; k++)
    multiply_projections(N_b, N_a, N_g, P_g2b, prj_g[k], I_a2g,
                         interp_a2b->prj[k], M_b, interp_a2b->mass_prj[k], lw_b,
                         interp_a2b->wi_mass_prj[k]);
  if (interp_b2a)
    for (bfam_locidx_t k = 0; k < 5; k++)
      multiply_projections(N_a, N_b, N_g, P_g2a, prj_g[k], I_b2g,
                           interp_b2a->prj[k], M_a, interp_b2a->mass_prj[k],
                           lw_a, interp_b2a->wi_mass_prj[k]);

  if (I_a2g)
    bfam_free_aligned(I_a2g);
  if (P_g2a)
    bfam_free_aligned(P_g2a);
  if (I_b2g)
    bfam_free_aligned(I_b2g);
  if (P_g2b)
    bfam_free_aligned(P_g2b);
  bfam_free_aligned(prj_g[1]);
  bfam_free_aligned(prj_g[2]);
  bfam_free_aligned(prj_g[3]);
  bfam_free_aligned(prj_g[4]);
}

/** return a pointer to the interpolator for the given orders (created and
 *  stored if it does not exsit)
 *
 * \param [in,out] N2N    dictionary for interpolator pointers
 * \param [in] N_src      source order
 * \param [in] N_dst      destination order
 * \param [in] inDim      dimensions of the subdomain
 *
 * \return pointer to bfam_subdomain_dgx_interpolator_t
 */
static bfam_subdomain_dgx_interpolator_t *
bfam_subdomain_dgx_get_interpolator(bfam_dictionary_t *N2N, const int N_src,
                                    const int N_dst, const int inDIM)
{
  char str[BFAM_BUFSIZ];
  snprintf(str, BFAM_BUFSIZ, "%d_to_%d", N_src, N_dst);

  bfam_subdomain_dgx_interpolator_t *interp =
      (bfam_subdomain_dgx_interpolator_t *)bfam_dictionary_get_value_ptr(N2N,
                                                                         str);
  if (!interp)
  {
    interp = bfam_malloc(sizeof(bfam_subdomain_dgx_interpolator_t));
    bfam_subdomain_dgx_interpolator_t *interp2 = NULL;
    if (N_src != N_dst)
      interp2 = bfam_malloc(sizeof(bfam_subdomain_dgx_interpolator_t));
    create_interpolators(interp, interp2, N_src, N_dst);
    BFAM_VERBOSE(">>>>>> Interpolator `%s' created", str);
    int rval = bfam_dictionary_insert_ptr(N2N, str, interp);
    BFAM_ABORT_IF_NOT(rval != 1, "Error inserting `%s` in N2N dict", str);

    if (interp2)
    {
      char str2[BFAM_BUFSIZ];
      snprintf(str2, BFAM_BUFSIZ, "%d_to_%d", N_dst, N_src);
      BFAM_VERBOSE(">>>>>> Interpolator `%s' created", str2);
      rval = bfam_dictionary_insert_ptr(N2N, str2, interp2);
      BFAM_ABORT_IF_NOT(rval != 1, "Error inserting `%s` in N2N dict", str);
    }
  }

  return interp;
}

typedef struct bfam_subdomain_dgx_get_put_data
{
  bfam_subdomain_dgx_t *sub;
  bfam_real_t *buffer;
  size_t size;
  size_t field;
} bfam_subdomain_dgx_get_put_data_t;

typedef struct bfam_subdomain_comm_args
{
  const char **scalars_m; /* \c NULL terminated array of scalars to
                           * send */

  const char **vectors_m; /* \c NULL terminated array of vectors to
                           * send
                           * \note will suffix n for normal component
                           *       and p[1-3] for perpendicular
                           *       components */

  const char **vector_components_m; /* \c NULL terminated array of vectors
                                     * components
                                     * \note must be three components per
                                     *       vector to send and a \c NULL entry
                                     *       will lead to 0 being used for that
                                     *       component */

  const char **tensors_m; /* \c NULL terminated array of tensor to
                           * send
                           * \note will suffix n for normal component
                           *       and p[1-3] for perpendicular
                           *       components */

  const char **tensor_components_m; /* \c NULL terminated array of symetric
                                     * tensor components in order
                                     * {11,22,33,12,13,23}
                                     * \note must be six components per tensor
                                     *       to send and \c NULL entry will
                                     *       lead to 0 being used for that
                                     *       component */

  const char **face_scalars_m; /* \c NULL terminated array of face scalars
                                *  to send */
  const char **scalars_p;
  const char **vectors_p;
  const char **vector_components_p;
  const char **tensors_p;
  const char **tensor_components_p;
  const char **face_scalars_p;

  /**< user specified glue grid communication info:
   *   recv_sz and send_sz should be added to and not reset
   */
  void (*user_comm_info)(struct bfam_subdomain *thisSubdomain, size_t *send_sz,
                         size_t *recv_sz, void *args);

  /**< user specified put data into the send buffer */
  void (*user_put_send_buffer)(struct bfam_subdomain *thisSubdomain,
                               void *buffer, size_t send_sz, void *args);

  /**< use specified get data from the recv buffer */
  void (*user_get_recv_buffer)(struct bfam_subdomain *thisSubdomain,
                               void *buffer, size_t recv_sz, void *args);

  /**< all the user to pass data */
  void *user_data;

  /**< callback function for custom user prefix */
  void (*user_prefix_function)(struct bfam_subdomain *thisSubdomain,
                               char *prefix, size_t buf_siz, void *user_data);

} bfam_subdomain_comm_args_t;

static void bfam_subdomain_dgx_comm_info(bfam_subdomain_t *thisSubdomain,
                                         int *rank, bfam_locidx_t *sort,
                                         int num_sort, size_t *send_sz,
                                         size_t *recv_sz, void *comm_args)
{
  BFAM_ASSERT(num_sort > 2);
  bfam_subdomain_dgx_t *sub = (bfam_subdomain_dgx_t *)thisSubdomain;

  BFAM_ASSERT(sub->base.glue_m && sub->base.glue_p);

  *rank = sub->base.glue_p->rank;
  sort[0] = sub->base.glue_p->id; /* neighbor ID */
  sort[1] = sub->base.glue_m->id; /* my ID */
  sort[2] = sub->base.uid;        /* my user ID */

  size_t send_num = sub->base.glue_m->fields.num_entries * sub->K * sub->Np;
  size_t recv_num = sub->base.glue_p->fields.num_entries * sub->K * sub->Np;
  BFAM_ASSERT(send_num == recv_num);

  if (comm_args != NULL)
  {
    bfam_subdomain_comm_args_t *args = (bfam_subdomain_comm_args_t *)comm_args;

    int count = 0;
    for (int i = 0; args->scalars_m[i] != NULL; i++)
      count++;
    for (int i = 0; args->vectors_m[i] != NULL; i++)
      count += 4;
    for (int i = 0; args->tensors_m[i] != NULL; i++)
      count += 4;
    for (int i = 0; args->face_scalars_m[i] != NULL; i++)
      count++;
    send_num = count * sub->K * sub->Np;

    count = 0;
    for (int i = 0; args->scalars_p[i] != NULL; i++)
      count++;
    for (int i = 0; args->vectors_p[i] != NULL; i++)
      count += 4;
    for (int i = 0; args->tensors_p[i] != NULL; i++)
      count += 4;
    for (int i = 0; args->face_scalars_p[i] != NULL; i++)
      count++;
    recv_num = count * sub->K * sub->Np;
  }

  *send_sz = send_num * sizeof(bfam_real_t);
  *recv_sz = recv_num * sizeof(bfam_real_t);

  if (comm_args != NULL)
  {
    bfam_subdomain_comm_args_t *args = (bfam_subdomain_comm_args_t *)comm_args;

    if (args->user_comm_info)
      args->user_comm_info(thisSubdomain, send_sz, recv_sz, comm_args);
  }

  BFAM_LDEBUG(
      " rank %3d   ns %3jd   ms %3jd   uid %3jd   send_sz %3zd   recv_sz %3zd",
      *rank, (intmax_t)sort[0], (intmax_t)sort[1], (intmax_t)sort[2], *send_sz,
      *recv_sz);
}

static void bfam_subdomain_dgx_vtk_interp(bfam_locidx_t K, int N_d,
                                          bfam_real_t *restrict d, int N_s,
                                          const bfam_real_t *restrict s,
                                          const bfam_real_t *restrict interp,
                                          int inDIM)
{
  BFAM_ASSUME_ALIGNED(d, 32);
  BFAM_ASSUME_ALIGNED(s, 32);
  BFAM_ASSUME_ALIGNED(interp, 32);
  for (bfam_locidx_t elem = 0; elem < K; ++elem)
  {
    const int Np_d = bfam_ipow(N_d + 1, inDIM);
    const int Np_s = bfam_ipow(N_s + 1, inDIM);

    const bfam_locidx_t o_d = elem * Np_d;
    const bfam_locidx_t o_s = elem * Np_s;

    for (bfam_locidx_t n = 0; n < Np_d; n++)
      d[o_d + n] = 0;
    if (inDIM == 1)
      for (int l = 0; l < N_s + 1; l++)
        for (int i = 0; i < N_d + 1; i++)
          d[o_d + i] += interp[(N_d + 1) * l + i] * s[o_s + l];
    else if (inDIM == 2)
      for (int m = 0; m < N_s + 1; m++)
        for (int l = 0; l < N_s + 1; l++)
          for (int j = 0; j < N_d + 1; j++)
            for (int i = 0; i < N_d + 1; i++)
              d[o_d + j * (N_d + 1) + i] += interp[(N_d + 1) * m + j] *
                                            interp[(N_d + 1) * l + i] *
                                            s[o_s + m * (N_s + 1) + l];
    else if (inDIM == 3)
      for (int n = 0; n < N_s + 1; n++)
        for (int m = 0; m < N_s + 1; m++)
          for (int l = 0; l < N_s + 1; l++)
            for (int k = 0; k < N_d + 1; k++)
              for (int j = 0; j < N_d + 1; j++)
                for (int i = 0; i < N_d + 1; i++)
                  d[o_d + k * (N_d + 1) * (N_d + 1) + j * (N_d + 1) + i] +=
                      interp[(N_d + 1) * l + i] * interp[(N_d + 1) * m + j] *
                      interp[(N_d + 1) * n + k] *
                      s[o_s + n * (N_s + 1) * (N_s + 1) + m * (N_s + 1) + l];
    else
      BFAM_ABORT("Cannot handle dim = %d", inDIM);
  }
}

/** Utility function to write binary data in VTK format.
 *
 * Currently this is just a wrapper to call a similar function in libsc.
 *
 * \param [in]  compressed boolean specifying if the binary data should be
 *                         compressed.
 * \param [out] file       stream to write the data to.
 * \param [in]  data       data to write out.
 * \param [in]  size       size of the data in bytes.
 *
 * \returns 0 on success and -1 on file error.
 */
static int bfam_vtk_write_binary_data(int compressed, FILE *file, char *data,
                                      size_t size)
{
  if (compressed)
    return sc_vtk_write_compressed(file, data, size);
  else
    return sc_vtk_write_binary(file, data, size);
}

/** Utility function for writing a scalar data array.
 *
 * \param [out] file            stream to write the scalar to.
 * \param [in]  name            name of the scalar.
 * \param [in]  writeBinary     boolean indicating if the data should be written
 *                              in binary.
 * \param [in]  writeCompressed boolean indicating if the data should be
 *                              compressed.
 * \param [in]  Ntotal          length of the scalar.
 * \param [in]  s               scalar data.
 */
static void bfam_vtk_write_real_scalar_data_array(FILE *file, const char *name,
                                                  int writeBinary,
                                                  int writeCompressed,
                                                  bfam_locidx_t Ntotal,
                                                  const bfam_real_t *s)
{
  const char *format;

  if (writeBinary)
    format = "binary";
  else
    format = "ascii";

  fprintf(file, "        <DataArray type=\"%s\" Name=\"%s\" format=\"%s\">\n",
          BFAM_REAL_VTK, name, format);
  if (writeBinary)
  {
    fprintf(file, "          ");
    int rval = bfam_vtk_write_binary_data(writeCompressed, file, (char *)s,
                                          Ntotal * sizeof(bfam_real_t));
    fprintf(file, "\n");
    if (rval)
      BFAM_WARNING("Error encoding %s", name);
  }
  else
  {
    for (bfam_locidx_t n = 0; n < Ntotal; ++n)
      fprintf(file, "         %" BFAM_REAL_FMTe "\n", s[n]);
  }

  fprintf(file, "        </DataArray>\n");
}

/** Utility function for writing a vector data array.
 *
 * \param [out] file            stream to write the vector to.
 * \param [in]  name            name of the vector.
 * \param [in]  writeBinary     boolean indicating if the data should be written
 *                              in binary.
 * \param [in]  writeCompressed boolean indicating if the data should be
 *                              compressed.
 * \param [in]  Ntotal          length of the vector components.
 * \param [in]  v1              1st component of the vector.
 * \param [in]  v2              2nd component of the vector.
 * \param [in]  v3              3rd component of the vector.
 */
static void bfam_vtk_write_real_vector_data_array(
    FILE *file, const char *name, int writeBinary, int writeCompressed,
    bfam_locidx_t Ntotal, const bfam_real_t *v1, const bfam_real_t *v2,
    const bfam_real_t *v3)
{
  const char *format;

  if (writeBinary)
    format = "binary";
  else
    format = "ascii";

  fprintf(file, "        <DataArray type=\"%s\" Name=\"%s\""
                " NumberOfComponents=\"3\" format=\"%s\">\n",
          BFAM_REAL_VTK, name, format);
  if (writeBinary)
  {
    size_t vSize = 3 * Ntotal * sizeof(bfam_real_t);
    bfam_real_t *v = bfam_malloc_aligned(vSize);

    for (bfam_locidx_t n = 0; n < Ntotal; ++n)
    {
      v[3 * n + 0] = (v1 != NULL) ? v1[n] : 0;
      v[3 * n + 1] = (v2 != NULL) ? v2[n] : 0;
      v[3 * n + 2] = (v3 != NULL) ? v3[n] : 0;
    }

    fprintf(file, "          ");
    int rval =
        bfam_vtk_write_binary_data(writeCompressed, file, (char *)v, vSize);
    fprintf(file, "\n");
    if (rval)
      BFAM_WARNING("Error encoding %s", name);

    bfam_free_aligned(v);
  }
  else
  {
    for (bfam_locidx_t n = 0; n < Ntotal; ++n)
    {
      bfam_real_t a = (v1 != NULL) ? v1[n] : 0;
      bfam_real_t b = (v2 != NULL) ? v2[n] : 0;
      bfam_real_t c = (v3 != NULL) ? v3[n] : 0;

      fprintf(file, "         %" BFAM_REAL_FMTe " %" BFAM_REAL_FMTe
                    " %" BFAM_REAL_FMTe "\n",
              a, b, c);
    }
  }

  fprintf(file, "        </DataArray>\n");
}

static int bfam_subdomain_dgx_vtk_write_vtu_piece(
    bfam_subdomain_t *subdomain, FILE *file, bfam_real_t time,
    const char **scalars, const char **vectors, const char **components,
    int writeBinary, int writeCompressed, int rank, bfam_locidx_t id,
    int Np_write)
{
  bfam_subdomain_dgx_t *sub = (bfam_subdomain_dgx_t *)subdomain;

  BFAM_LDEBUG("Handling vtk for subdomain %s", subdomain->name);

  const char *format;

  if (writeBinary)
    format = "binary";
  else
    format = "ascii";

  const bfam_locidx_t K = sub->K;
  int N_vtk = sub->N;
  int Np_vtk = sub->Np;
  bfam_real_t *interp = NULL;

  bfam_real_t *restrict stor1 = NULL;
  bfam_real_t *restrict stor2 = NULL;
  bfam_real_t *restrict stor3 = NULL;

  if (Np_write > 0)
  {
    BFAM_ABORT_IF_NOT(Np_write > 1, "Np_write = %d is not valid", Np_write);

    N_vtk = Np_write - 1;
    Np_vtk = bfam_ipow(Np_write, sub->dim);

    interp =
        bfam_malloc_aligned(sizeof(bfam_real_t) * (sub->N + 1) * (N_vtk + 1));

    bfam_long_real_t *cal_interp = bfam_malloc_aligned(
        sizeof(bfam_long_real_t) * (sub->N + 1) * (N_vtk + 1));
    bfam_long_real_t *lr =
        bfam_malloc_aligned(sizeof(bfam_long_real_t) * Np_write);

    for (int r = 0; r < Np_write; r++)
      lr[r] = -1 + 2 * (bfam_long_real_t)r / (Np_write - 1);

    bfam_jacobi_p_interpolation(0, 0, sub->N, Np_write, lr, sub->lV,
                                cal_interp);

    for (int n = 0; n < (sub->N + 1) * (N_vtk + 1); n++)
      interp[n] = (bfam_real_t)cal_interp[n];

    stor1 = bfam_malloc_aligned(sizeof(bfam_real_t) * Np_vtk * K);
    stor2 = bfam_malloc_aligned(sizeof(bfam_real_t) * Np_vtk * K);
    stor3 = bfam_malloc_aligned(sizeof(bfam_real_t) * Np_vtk * K);

    bfam_free_aligned(lr);
    bfam_free_aligned(cal_interp);
  }

  const int Ncorners = sub->Ng[sub->numg - 1];

  const bfam_locidx_t Ncells = K * bfam_ipow(N_vtk, sub->dim);
  const bfam_locidx_t Ntotal = K * Np_vtk;

  bfam_real_t *restrict x =
      bfam_dictionary_get_value_ptr(&subdomain->fields, "_grid_x0");
  bfam_real_t *restrict y =
      bfam_dictionary_get_value_ptr(&subdomain->fields, "_grid_x1");
  bfam_real_t *restrict z =
      bfam_dictionary_get_value_ptr(&subdomain->fields, "_grid_x2");

  if (interp == NULL)
  {
    stor1 = x;
    stor2 = y;
    stor3 = z;
  }
  else
  {
    bfam_subdomain_dgx_vtk_interp(K, N_vtk, stor1, sub->N, x, interp, sub->dim);
    bfam_subdomain_dgx_vtk_interp(K, N_vtk, stor2, sub->N, y, interp, sub->dim);
    if (z != NULL)
      bfam_subdomain_dgx_vtk_interp(K, N_vtk, stor3, sub->N, z, interp,
                                    sub->dim);
    else
      for (bfam_locidx_t k = 0; k < Np_vtk * K; k++)
        stor3[k] = 0;
  }

  fprintf(file, "    <Piece NumberOfPoints=\"%jd\" NumberOfCells=\"%jd\">\n",
          (intmax_t)Ntotal, (intmax_t)Ncells);

  /*
   * Points
   */
  fprintf(file, "      <Points>\n");

  bfam_vtk_write_real_vector_data_array(file, "Position", writeBinary,
                                        writeCompressed, Ntotal, stor1, stor2,
                                        stor3);

  fprintf(file, "      </Points>\n");

  /*
   * Cells
   */
  fprintf(file, "      <Cells>\n");

  /*
   * Connectivity
   */
  fprintf(file, "        <DataArray type=\"%s\" Name=\"connectivity\""
                " format=\"%s\">\n",
          BFAM_LOCIDX_VTK, format);
  if (writeBinary)
  {
    size_t cellsSize = Ncells * Ncorners * sizeof(bfam_locidx_t);
    bfam_locidx_t *cells = bfam_malloc_aligned(cellsSize);

    if (sub->dim == 1)
      for (bfam_locidx_t k = 0, i = 0; k < K; ++k)
      {
        for (int n = 0; n < N_vtk; ++n)
        {
          cells[i++] = Np_vtk * k + (n + 0);
          cells[i++] = Np_vtk * k + (n + 1);
        }
      }
    else if (sub->dim == 2)
      for (bfam_locidx_t k = 0, i = 0; k < K; ++k)
      {
        for (int m = 0; m < N_vtk; ++m)
        {
          for (int n = 0; n < N_vtk; ++n)
          {
            cells[i++] = Np_vtk * k + (N_vtk + 1) * (m + 0) + (n + 0);
            cells[i++] = Np_vtk * k + (N_vtk + 1) * (m + 0) + (n + 1);
            cells[i++] = Np_vtk * k + (N_vtk + 1) * (m + 1) + (n + 0);
            cells[i++] = Np_vtk * k + (N_vtk + 1) * (m + 1) + (n + 1);
          }
        }
      }
    else if (sub->dim == 3)
      for (bfam_locidx_t k = 0, i = 0; k < K; ++k)
      {
        for (int l = 0; l < N_vtk; ++l)
        {
          for (int m = 0; m < N_vtk; ++m)
          {
            for (int n = 0; n < N_vtk; ++n)
            {
              cells[i++] = Np_vtk * k + (N_vtk + 1) * (N_vtk + 1) * (l + 0) +
                           (N_vtk + 1) * (m + 0) + (n + 0);
              cells[i++] = Np_vtk * k + (N_vtk + 1) * (N_vtk + 1) * (l + 0) +
                           (N_vtk + 1) * (m + 0) + (n + 1);
              cells[i++] = Np_vtk * k + (N_vtk + 1) * (N_vtk + 1) * (l + 0) +
                           (N_vtk + 1) * (m + 1) + (n + 0);
              cells[i++] = Np_vtk * k + (N_vtk + 1) * (N_vtk + 1) * (l + 0) +
                           (N_vtk + 1) * (m + 1) + (n + 1);
              cells[i++] = Np_vtk * k + (N_vtk + 1) * (N_vtk + 1) * (l + 1) +
                           (N_vtk + 1) * (m + 0) + (n + 0);
              cells[i++] = Np_vtk * k + (N_vtk + 1) * (N_vtk + 1) * (l + 1) +
                           (N_vtk + 1) * (m + 0) + (n + 1);
              cells[i++] = Np_vtk * k + (N_vtk + 1) * (N_vtk + 1) * (l + 1) +
                           (N_vtk + 1) * (m + 1) + (n + 0);
              cells[i++] = Np_vtk * k + (N_vtk + 1) * (N_vtk + 1) * (l + 1) +
                           (N_vtk + 1) * (m + 1) + (n + 1);
            }
          }
        }
      }
    else
      BFAM_ABORT("not implemented for dim = %d", sub->dim);

    fprintf(file, "          ");
    int rval = bfam_vtk_write_binary_data(writeCompressed, file, (char *)cells,
                                          cellsSize);
    fprintf(file, "\n");
    if (rval)
      BFAM_WARNING("Error encoding cells");

    bfam_free_aligned(cells);
  }
  else
  {
    if (sub->dim == 1)
      for (bfam_locidx_t k = 0; k < K; ++k)
        for (int n = 0; n < N_vtk; ++n)
          fprintf(file, "          %8jd %8jd\n", (intmax_t)Np_vtk * k + (n + 0),
                  (intmax_t)Np_vtk * k + (n + 1));
    else if (sub->dim == 2)
      for (bfam_locidx_t k = 0; k < K; ++k)
        for (int m = 0; m < N_vtk; ++m)
          for (int n = 0; n < N_vtk; ++n)
            fprintf(file, "          %8jd %8jd %8jd %8jd\n",
                    (intmax_t)Np_vtk * k + (N_vtk + 1) * (m + 0) + (n + 0),
                    (intmax_t)Np_vtk * k + (N_vtk + 1) * (m + 0) + (n + 1),
                    (intmax_t)Np_vtk * k + (N_vtk + 1) * (m + 1) + (n + 0),
                    (intmax_t)Np_vtk * k + (N_vtk + 1) * (m + 1) + (n + 1));
    else if (sub->dim == 3)
      for (bfam_locidx_t k = 0; k < K; ++k)
        for (int l = 0; l < N_vtk; ++l)
          for (int m = 0; m < N_vtk; ++m)
            for (int n = 0; n < N_vtk; ++n)
              fprintf(
                  file, "          %8jd %8jd %8jd %8jd %8jd %8jd %8jd %8jd\n",
                  (intmax_t)Np_vtk * k + (N_vtk + 1) * (N_vtk + 1) * (l + 0) +
                      (N_vtk + 1) * (m + 0) + (n + 0),
                  (intmax_t)Np_vtk * k + (N_vtk + 1) * (N_vtk + 1) * (l + 0) +
                      (N_vtk + 1) * (m + 0) + (n + 1),
                  (intmax_t)Np_vtk * k + (N_vtk + 1) * (N_vtk + 1) * (l + 0) +
                      (N_vtk + 1) * (m + 1) + (n + 0),
                  (intmax_t)Np_vtk * k + (N_vtk + 1) * (N_vtk + 1) * (l + 0) +
                      (N_vtk + 1) * (m + 1) + (n + 1),
                  (intmax_t)Np_vtk * k + (N_vtk + 1) * (N_vtk + 1) * (l + 1) +
                      (N_vtk + 1) * (m + 0) + (n + 0),
                  (intmax_t)Np_vtk * k + (N_vtk + 1) * (N_vtk + 1) * (l + 1) +
                      (N_vtk + 1) * (m + 0) + (n + 1),
                  (intmax_t)Np_vtk * k + (N_vtk + 1) * (N_vtk + 1) * (l + 1) +
                      (N_vtk + 1) * (m + 1) + (n + 0),
                  (intmax_t)Np_vtk * k + (N_vtk + 1) * (N_vtk + 1) * (l + 1) +
                      (N_vtk + 1) * (m + 1) + (n + 1));
    else
      BFAM_ABORT("not implemented for dim = %d %d", sub->dim);
  }
  fprintf(file, "        </DataArray>\n");

  /*
   * Offsets
   */
  fprintf(file, "        <DataArray type=\"%s\" Name=\"offsets\""
                " format=\"%s\">\n",
          BFAM_LOCIDX_VTK, format);
  fprintf(file, "          ");
  if (writeBinary)
  {
    size_t offsetsSize = Ncells * sizeof(bfam_locidx_t);
    bfam_locidx_t *offsets = bfam_malloc_aligned(offsetsSize);

    for (bfam_locidx_t i = 1; i <= Ncells; ++i)
      offsets[i - 1] = Ncorners * i;

    int rval = bfam_vtk_write_binary_data(writeCompressed, file,
                                          (char *)offsets, offsetsSize);
    if (rval)
      BFAM_WARNING("Error encoding offsets");

    bfam_free_aligned(offsets);
  }
  else
  {
    for (bfam_locidx_t i = 1, sk = 1; i <= Ncells; ++i, ++sk)
    {
      fprintf(file, " %8jd", (intmax_t)(Ncorners * i));
      if (!(sk % 20) && i != Ncells)
        fprintf(file, "\n          ");
    }
  }
  fprintf(file, "\n");
  fprintf(file, "        </DataArray>\n");

  /*
   * Types
   */
  fprintf(file, "        <DataArray type=\"UInt8\" Name=\"types\""
                " format=\"%s\">\n",
          format);
  fprintf(file, "          ");
  if (writeBinary)
  {
    size_t typesSize = Ncells * sizeof(uint8_t);
    uint8_t *types = bfam_malloc_aligned(typesSize);

    if (sub->dim == 1)
      for (bfam_locidx_t i = 0; i < Ncells; ++i)
        types[i] = 3; /* VTK_LINE */
    else if (sub->dim == 2)
      for (bfam_locidx_t i = 0; i < Ncells; ++i)
        types[i] = 8; /* VTK_PIXEL */
    else if (sub->dim == 3)
      for (bfam_locidx_t i = 0; i < Ncells; ++i)
        types[i] = 11; /* VTK_VOXEL */
    else
      BFAM_ABORT("cannot handle dim = %d", sub->dim);

    int rval = bfam_vtk_write_binary_data(writeCompressed, file, (char *)types,
                                          typesSize);
    if (rval)
      BFAM_WARNING("Error encoding types");

    bfam_free_aligned(types);
  }
  else
  {
    for (bfam_locidx_t i = 0, sk = 1; i < Ncells; ++i, ++sk)
    {
      if (sub->dim == 1)
        fprintf(file, " 3"); /* VTK_LINE */
      else if (sub->dim == 2)
        fprintf(file, " 8"); /* VTK_PIXEL */
      else if (sub->dim == 3)
        fprintf(file, " 11"); /* VTK_VOXEL */
      else
        BFAM_ABORT("cannot handle dim = %d", sub->dim);
      if (!(sk % 20) && i != (Ncells - 1))
        fprintf(file, "\n         ");
    }
  }
  fprintf(file, "\n");
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "      </Cells>\n");

  /*
   * Cell Data
   */
  fprintf(file,
          "      <CellData Scalars=\"time,mpirank,subdomain_id,root_id\">\n");
  fprintf(file, "        <DataArray type=\"%s\" Name=\"time\""
                " format=\"%s\">\n",
          BFAM_REAL_VTK, format);
  fprintf(file, "          ");
  if (writeBinary)
  {
    size_t timesize = Ncells * sizeof(bfam_real_t);
    bfam_real_t *times = bfam_malloc_aligned(timesize);

    for (bfam_locidx_t i = 0; i < Ncells; ++i)
      times[i] = time;

    int rval = bfam_vtk_write_binary_data(writeCompressed, file, (char *)times,
                                          timesize);
    if (rval)
      BFAM_WARNING("Error encoding times");

    bfam_free_aligned(times);
  }
  else
  {
    for (bfam_locidx_t i = 0, sk = 1; i < Ncells; ++i, ++sk)
    {
      fprintf(file, " %" BFAM_REAL_FMTe, time);
      if (!(sk % 8) && i != (Ncells - 1))
        fprintf(file, "\n         ");
    }
  }
  fprintf(file, "\n");
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "        <DataArray type=\"%s\" Name=\"mpirank\""
                " format=\"%s\">\n",
          BFAM_LOCIDX_VTK, format);
  fprintf(file, "          ");
  if (writeBinary)
  {
    size_t ranksSize = Ncells * sizeof(bfam_locidx_t);
    bfam_locidx_t *ranks = bfam_malloc_aligned(ranksSize);

    for (bfam_locidx_t i = 0; i < Ncells; ++i)
      ranks[i] = rank;

    int rval = bfam_vtk_write_binary_data(writeCompressed, file, (char *)ranks,
                                          ranksSize);
    if (rval)
      BFAM_WARNING("Error encoding ranks");

    bfam_free_aligned(ranks);
  }
  else
  {
    for (bfam_locidx_t i = 0, sk = 1; i < Ncells; ++i, ++sk)
    {
      fprintf(file, " %6jd", (intmax_t)rank);
      if (!(sk % 8) && i != (Ncells - 1))
        fprintf(file, "\n         ");
    }
  }
  fprintf(file, "\n");
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "        <DataArray type=\"%s\" Name=\"subdomain_id\""
                " format=\"%s\">\n",
          BFAM_LOCIDX_VTK, format);
  fprintf(file, "          ");
  if (writeBinary)
  {
    size_t idsSize = Ncells * sizeof(bfam_locidx_t);
    bfam_locidx_t *ids = bfam_malloc_aligned(idsSize);

    for (bfam_locidx_t i = 0; i < Ncells; ++i)
      ids[i] = subdomain->id;

    int rval =
        bfam_vtk_write_binary_data(writeCompressed, file, (char *)ids, idsSize);
    if (rval)
      BFAM_WARNING("Error encoding ids");

    bfam_free_aligned(ids);
  }
  else
  {
    for (bfam_locidx_t i = 0, sk = 1; i < Ncells; ++i, ++sk)
    {
      fprintf(file, " %6jd", (intmax_t)subdomain->id);
      if (!(sk % 8) && i != (Ncells - 1))
        fprintf(file, "\n         ");
    }
  }
  fprintf(file, "\n");
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "        <DataArray type=\"%s\" Name=\"root_id\""
                " format=\"%s\">\n",
          BFAM_LOCIDX_VTK, format);
  fprintf(file, "          ");
  if (writeBinary)
  {
    size_t idsSize = Ncells * sizeof(bfam_locidx_t);
    bfam_locidx_t *ids = bfam_malloc_aligned(idsSize);

    for (bfam_locidx_t i = 0; i < Ncells; ++i)
      ids[i] = subdomain->uid;

    int rval =
        bfam_vtk_write_binary_data(writeCompressed, file, (char *)ids, idsSize);
    if (rval)
      BFAM_WARNING("Error encoding ids");

    bfam_free_aligned(ids);
  }
  else
  {
    for (bfam_locidx_t i = 0, sk = 1; i < Ncells; ++i, ++sk)
    {
      fprintf(file, " %6jd", (intmax_t)subdomain->uid);
      if (!(sk % 8) && i != (Ncells - 1))
        fprintf(file, "\n         ");
    }
  }
  fprintf(file, "\n");
  fprintf(file, "        </DataArray>\n");

  fprintf(file, "      </CellData>\n");

  char pointscalars[BFAM_BUFSIZ];
  bfam_util_strcsl(pointscalars, scalars);

  char pointvectors[BFAM_BUFSIZ];
  bfam_util_strcsl(pointvectors, vectors);

  fprintf(file, "      <PointData Scalars=\"%s\" Vectors=\"%s\">\n",
          pointscalars, pointvectors);

  if (scalars)
  {
    for (size_t s = 0; scalars[s]; ++s)
    {
      bfam_real_t *sdata =
          bfam_dictionary_get_value_ptr(&subdomain->fields, scalars[s]);
      BFAM_ABORT_IF(sdata == NULL, "VTK: Field %s not in subdomain %s",
                    scalars[s], subdomain->name);
      if (interp == NULL)
      {
        stor1 = sdata;
      }
      else
      {
        bfam_subdomain_dgx_vtk_interp(K, N_vtk, stor1, sub->N, sdata, interp,
                                      sub->dim);
      }

      bfam_vtk_write_real_scalar_data_array(file, scalars[s], writeBinary,
                                            writeCompressed, Ntotal, stor1);
    }
  }

  if (vectors)
  {
    for (size_t v = 0; vectors[v]; ++v)
    {

      bfam_real_t *v1 = bfam_dictionary_get_value_ptr(&subdomain->fields,
                                                      components[3 * v + 0]);
      bfam_real_t *v2 = bfam_dictionary_get_value_ptr(&subdomain->fields,
                                                      components[3 * v + 1]);
      bfam_real_t *v3 = bfam_dictionary_get_value_ptr(&subdomain->fields,
                                                      components[3 * v + 2]);

      BFAM_ABORT_IF(v1 == NULL, "VTK: Field %s not in subdomain %s",
                    components[3 * v + 0], subdomain->name);
      BFAM_ABORT_IF(v2 == NULL, "VTK: Field %s not in subdomain %s",
                    components[3 * v + 1], subdomain->name);
      BFAM_ABORT_IF(v3 == NULL, "VTK: Field %s not in subdomain %s",
                    components[3 * v + 2], subdomain->name);
      if (interp == NULL)
      {
        stor1 = v1;
        stor2 = v2;
        stor3 = v3;
      }
      else
      {
        bfam_subdomain_dgx_vtk_interp(K, N_vtk, stor1, sub->N, v1, interp,
                                      sub->dim);
        bfam_subdomain_dgx_vtk_interp(K, N_vtk, stor2, sub->N, v2, interp,
                                      sub->dim);
        bfam_subdomain_dgx_vtk_interp(K, N_vtk, stor3, sub->N, v3, interp,
                                      sub->dim);
      }

      bfam_vtk_write_real_vector_data_array(file, vectors[v], writeBinary,
                                            writeCompressed, Ntotal, stor1,
                                            stor2, stor3);
    }
  }

  fprintf(file, "      </PointData>\n");
  fprintf(file, "    </Piece>\n");

  if (interp != NULL)
  {
    bfam_free_aligned(interp);
    bfam_free_aligned(stor1);
    bfam_free_aligned(stor2);
    bfam_free_aligned(stor3);
  }
  return 1;
}

static int bfam_subdomain_dgx_field_add(bfam_subdomain_t *subdomain,
                                        const char *name)
{
  bfam_subdomain_dgx_t *s = (bfam_subdomain_dgx_t *)subdomain;

  if (bfam_dictionary_get_value_ptr(&s->base.fields, name))
    return 1;

  size_t fieldSize = s->Np * s->K * sizeof(bfam_real_t);
  bfam_real_t *field = bfam_malloc_aligned(fieldSize);
#ifdef BFAM_DEBUG
  for (int i = 0; i < s->Np * s->K; i++)
    field[i] = bfam_real_nan("");
#endif

  int rval = bfam_dictionary_insert_ptr(&s->base.fields, name, field);

  BFAM_ASSERT(rval != 1);

  if (rval == 0)
    bfam_free_aligned(field);

  return rval;
}

static inline int ***bfam_subdomain_dgx_gmask_set(const int numg, const int N,
                                                  int *Np, int *Ng, int *Ngp,
                                                  int inDIM)
{
  BFAM_ABORT_IF(inDIM > 3 || inDIM < 0,
                "bfam_subdomain_dgx_gmask_set cannot handle dim = %d", inDIM);

  if (inDIM == 0)
  {
    *Np = 1;
    return NULL;
  }

  /* this could probably be made generic for arbitrary dimensions, but until
   * that's needed... */
  switch (inDIM)
  {
  case 1:
    *Np = N + 1;

    /* just corners */
    Ng[0] = 2;
    Ngp[0] = 1;
    break;

  case 2:
    *Np = (N + 1) * (N + 1);

    /* edges */
    Ng[0] = 4;
    Ngp[0] = N + 1;

    /* corners */
    Ng[1] = 4;
    Ngp[1] = 1;
    break;

  case 3:
    *Np = (N + 1) * (N + 1) * (N + 1);

    /* faces */
    Ng[0] = 6;
    Ngp[0] = (N + 1) * (N + 1);

    /* edges */
    Ng[1] = 12;
    Ngp[1] = N + 1;

    /* corners */
    Ng[2] = 8;
    Ngp[2] = 1;
    break;

  default:
    BFAM_ABORT("cannot handle dim = %d", inDIM);
  }

  int ***gmask = bfam_malloc_aligned(numg * sizeof(int **));
  for (int g = 0; g < numg; g++)
  {
    gmask[g] = bfam_malloc_aligned(Ng[g] * sizeof(int *));
    for (int i = 0; i < Ng[g]; i++)
      gmask[g][i] = bfam_malloc_aligned(Ngp[g] * sizeof(int));
  }

  switch (inDIM)
  {
  case 1:
    gmask[0][0][0] = 0;
    gmask[0][1][0] = N;
    break;

  case 2:
    /* edges */
    for (int i = 0; i < N + 1; ++i)
      gmask[0][0][i] = i * (N + 1);
    for (int i = 0; i < N + 1; ++i)
      gmask[0][1][i] = (i + 1) * (N + 1) - 1;
    for (int i = 0; i < N + 1; ++i)
      gmask[0][2][i] = i;
    for (int i = 0; i < N + 1; ++i)
      gmask[0][3][i] = (N + 1) * N + i;

    /* corners */
    for (int j = 0; j < 2; ++j)
      for (int i = 0; i < 2; ++i)
        gmask[1][i + j * 2][0] = i * N + j * (N + 1) * N;
    break;

  case 3:
    /* This could all probably be cleaned up... */

    /* faces */
    {
      int n, i, j, k, f = -1;

      n = 0;
      i = 0;
      f++;
      for (k = 0; k < N + 1; k++)
        for (j = 0; j < N + 1; j++)
          gmask[0][f][n++] = i + j * (N + 1) + k * (N + 1) * (N + 1);

      n = 0;
      i = N;
      f++;
      for (k = 0; k < N + 1; k++)
        for (j = 0; j < N + 1; j++)
          gmask[0][f][n++] = i + j * (N + 1) + k * (N + 1) * (N + 1);

      n = 0;
      j = 0;
      f++;
      for (k = 0; k < N + 1; k++)
        for (i = 0; i < N + 1; i++)
          gmask[0][f][n++] = i + j * (N + 1) + k * (N + 1) * (N + 1);

      n = 0;
      j = N;
      f++;
      for (k = 0; k < N + 1; k++)
        for (i = 0; i < N + 1; i++)
          gmask[0][f][n++] = i + j * (N + 1) + k * (N + 1) * (N + 1);

      n = 0;
      k = 0;
      f++;
      for (j = 0; j < N + 1; j++)
        for (i = 0; i < N + 1; i++)
          gmask[0][f][n++] = i + j * (N + 1) + k * (N + 1) * (N + 1);

      n = 0;
      k = N;
      f++;
      for (j = 0; j < N + 1; j++)
        for (i = 0; i < N + 1; i++)
          gmask[0][f][n++] = i + j * (N + 1) + k * (N + 1) * (N + 1);
    }

    /* edges */
    {
      int n, i, j, k, e = 0;

      for (k = 0; k < N + 1; k += N)
        for (j = 0; j < N + 1; j += N)
        {
          n = 0;
          for (i = 0; i < N + 1; i++)
            gmask[1][e][n++] = i + j * (N + 1) + k * (N + 1) * (N + 1);
          e++;
        }
      for (k = 0; k < N + 1; k += N)
        for (i = 0; i < N + 1; i += N)
        {
          n = 0;
          for (j = 0; j < N + 1; j++)
            gmask[1][e][n++] = i + j * (N + 1) + k * (N + 1) * (N + 1);
          e++;
        }
      for (j = 0; j < N + 1; j += N)
        for (i = 0; i < N + 1; i += N)
        {
          n = 0;
          for (k = 0; k < N + 1; k++)
            gmask[1][e][n++] = i + j * (N + 1) + k * (N + 1) * (N + 1);
          e++;
        }
    }

    /* corners */
    for (int k = 0, c = 0; k < N + 1; k += N)
      for (int j = 0; j < N + 1; j += N)
        for (int i = 0; i < N + 1; i += N)
          gmask[2][c++][0] = i + j * (N + 1) + k * (N + 1) * (N + 1);

    break;

  default:
    BFAM_ABORT("cannot handle dim = %d", inDIM);
  }

  return gmask;
}

static void bfam_subdomain_dgx_buildmaps(
    int N, bfam_locidx_t K, int Np, int Nfp, int Nfaces,
    const bfam_locidx_t *EToE, const int8_t *EToF, int ***gmask,
    bfam_locidx_t *restrict vmapP, bfam_locidx_t *restrict vmapM, int inDIM)
{
  for (bfam_locidx_t k1 = 0, sk = 0; k1 < K; ++k1)
  {
    for (int8_t f1 = 0; f1 < Nfaces; ++f1)
    {
      bfam_locidx_t k2 = EToE[Nfaces * k1 + f1];
      int8_t f2 = (int8_t)(EToF[Nfaces * k1 + f1] % Nfaces);
      int8_t o = (int8_t)(EToF[Nfaces * k1 + f1] / Nfaces);

      for (int n = 0; n < Nfp; ++n)
      {
        vmapM[sk + n] = Np * k1 + gmask[0][f1][n];
      }

      switch (inDIM)
      {
      case 1:
        /* Orientation does not matter in 1D */
        for (int n = 0; n < Nfp; ++n)
        {
          vmapP[sk + n] = Np * k2 + gmask[0][f2][n];
        }
        break;
      case 2:
        for (int n = 0; n < Nfp; ++n)
        {
          if (o)
            vmapP[sk + n] = Np * k2 + gmask[0][f2][Nfp - 1 - n];
          else
            vmapP[sk + n] = Np * k2 + gmask[0][f2][n];
        }
        break;
      case 3:
      {
        const int Nrp = N + 1;
        BFAM_ASSERT(Nfp == Nrp * Nrp);

        int oidx = -1;

        const int8_t nr = bfam_p8est_FToF_code[f1][f2];
        const int8_t ns =
            (k1 == k2 && f1 == f2) ? 0 : bfam_p8est_code_to_perm[nr][o];

        for (int j = 0, n = 0; j < Nrp; ++j)
        {
          for (int i = 0; i < Nrp; ++i, ++n)
          {
            int ir = Nrp - (i + 1);
            int jr = Nrp - (j + 1);
            switch (ns)
            {
            case 0:
              oidx = i + j * Nrp;
              break;
            case 1:
              oidx = j + i * Nrp;
              break;
            case 2:
              oidx = ir + j * Nrp;
              break;
            case 3:
              oidx = jr + i * Nrp;
              break;
            case 4:
              oidx = j + ir * Nrp;
              break;
            case 5:
              oidx = i + jr * Nrp;
              break;
            case 6:
              oidx = jr + ir * Nrp;
              break;
            case 7:
              oidx = ir + jr * Nrp;
              break;
            default:
              BFAM_ABORT("invalid orientation %d", o);
            }
            vmapP[sk + n] = Np * k2 + gmask[0][f2][oidx];
          }
        }
      }
      break;
      default:
        BFAM_ABORT("cannot handle dim = %d", inDIM);
      }

      sk += Nfp;
    }
  }
}

/* set all subdomain values to something logical */
static void bfam_subdomain_dgx_null_all_values(bfam_subdomain_dgx_t *sub)
{
  sub->K = 0;
  sub->N = -1;
  sub->Np = 0;
  sub->Ngp = NULL;
  sub->numg = 0;
  sub->Ng = NULL;
  sub->r = NULL;
  sub->w = NULL;
  sub->wi = NULL;
  sub->Dr = NULL;
  sub->lDr = NULL;
  sub->lr = NULL;
  sub->lw = NULL;
  sub->lV = NULL;
  sub->K = 0;
  sub->vmapM = NULL;
  sub->vmapP = NULL;
  sub->gmask = NULL;
  sub->EToQ = NULL;
}

static int bfam_subdomain_dgx_free_fields(const char *key, void *val, void *arg)
{
  bfam_free_aligned(val);

  return 1;
}

static void bfam_subdomain_dgx_free_glue(bfam_subdomain_dgx_glue_data_t *glue)
{
  if (glue)
  {
    if (glue->interpolation)
      bfam_free_aligned(glue->interpolation);
    if (glue->massprojection)
      bfam_free_aligned(glue->massprojection);
    if (glue->projection)
      bfam_free_aligned(glue->projection);
    if (glue->exact_mass)
      bfam_free_aligned(glue->exact_mass);

    if (glue->EToEp)
      bfam_free_aligned(glue->EToEp);
    if (glue->EToHp)
      bfam_free_aligned(glue->EToHp);
    if (glue->EToEm)
      bfam_free_aligned(glue->EToEm);
    if (glue->EToFm)
      bfam_free_aligned(glue->EToFm);
    if (glue->EToHm)
      bfam_free_aligned(glue->EToHm);
    if (glue->EToOp)
      bfam_free_aligned(glue->EToOp);
    if (glue->mapOp)
    {
      for (int n = 0; n < glue->num_orient; n++)
        bfam_free_aligned(glue->mapOp[n]);
      bfam_free_aligned(glue->mapOp);
    }
    bfam_critbit0_clear(&glue->base.tags);
    bfam_free(glue);
  }
}

/** free up the memory allocated by the subdomain
 *
 * \param [in,out] subdomain subdomain to clean up
 */
static void bfam_subdomain_dgx_free(bfam_subdomain_t *thisSubdomain)
{
  bfam_subdomain_dgx_t *sub = (bfam_subdomain_dgx_t *)thisSubdomain;

  bfam_dictionary_allprefixed_ptr(&sub->base.fields, "",
                                  &bfam_subdomain_dgx_free_fields, NULL);
  if (sub->base.glue_p)
    bfam_dictionary_allprefixed_ptr(&sub->base.glue_p->fields, "",
                                    &bfam_subdomain_dgx_free_fields, NULL);
  if (sub->base.glue_m)
    bfam_dictionary_allprefixed_ptr(&sub->base.glue_m->fields, "",
                                    &bfam_subdomain_dgx_free_fields, NULL);
  bfam_dictionary_allprefixed_ptr(&sub->base.fields_face, "",
                                  &bfam_subdomain_dgx_free_fields, NULL);

  bfam_subdomain_free(thisSubdomain);

  if (sub->gmask)
  {
    for (int g = 0; g < sub->numg; g++)
    {
      for (int i = 0; i < sub->Ng[g]; i++)
        bfam_free_aligned(sub->gmask[g][i]);
      bfam_free_aligned(sub->gmask[g]);
    }
    bfam_free_aligned(sub->gmask);
    sub->gmask = NULL;
  }
  if (sub->Ng)
    bfam_free_aligned(sub->Ng);
  sub->Ng = NULL;
  if (sub->Ngp)
    bfam_free_aligned(sub->Ngp);
  sub->Ngp = NULL;

  if (sub->EToQ)
    bfam_free_aligned(sub->EToQ);
  sub->EToQ = NULL;
  if (sub->vmapP)
    bfam_free_aligned(sub->vmapP);
  sub->vmapP = NULL;
  if (sub->vmapM)
    bfam_free_aligned(sub->vmapM);
  sub->vmapM = NULL;

  if (sub->hadapt)
    bfam_free_aligned(sub->hadapt);
  sub->hadapt = NULL;
  if (sub->padapt)
    bfam_free_aligned(sub->padapt);
  sub->padapt = NULL;

  if (sub->q_id)
    bfam_free_aligned(sub->q_id);
  sub->q_id = NULL;

  if (sub->lvl)
    bfam_free_aligned(sub->lvl);
  sub->lvl = NULL;

  bfam_subdomain_dgx_null_all_values(sub);

  bfam_subdomain_dgx_free_glue(
      (bfam_subdomain_dgx_glue_data_t *)sub->base.glue_m);
  bfam_subdomain_dgx_free_glue(
      (bfam_subdomain_dgx_glue_data_t *)sub->base.glue_p);
}

static void
bfam_subdomain_dgx_generic_init(bfam_subdomain_dgx_t *subdomain,
                                const bfam_locidx_t id, const bfam_locidx_t uid,
                                const char *name, const int N,
                                const bfam_locidx_t K, bfam_dictionary_t *N2N,
                                bfam_dictionary_t *dgx_ops, const int inDIM)
{
  BFAM_ABORT_IF(inDIM < 0, "dimension %d is not possible in bfam", inDIM);
  BFAM_ABORT_IF(inDIM == 0 && N != 0,
                "if inDIM < 1 then N must be zero (i.e., constant");

  bfam_subdomain_init(&subdomain->base, id, uid, name);
  bfam_subdomain_add_tag(&subdomain->base, "_subdomain_dgx");
  char dim_str[BFAM_BUFSIZ];
  snprintf(dim_str, BFAM_BUFSIZ, "_dimension_%d", inDIM);
  bfam_subdomain_add_tag(&subdomain->base, dim_str);

  bfam_subdomain_dgx_null_all_values(subdomain);

  subdomain->dim = inDIM;

  subdomain->base.free = bfam_subdomain_dgx_free;
  subdomain->base.vtk_write_vtu_piece = bfam_subdomain_dgx_vtk_write_vtu_piece;
  subdomain->base.field_add = bfam_subdomain_dgx_field_add;
  subdomain->base.glue_comm_info = bfam_subdomain_dgx_comm_info;

  subdomain->numg = inDIM;
  const int numg = subdomain->numg;

  int *Ng = NULL;
  if (numg > 0)
    Ng = bfam_malloc_aligned(sizeof(int) * numg);
  subdomain->Ng = Ng;

  int *Ngp = NULL;
  if (numg > 0)
    Ngp = bfam_malloc_aligned(sizeof(int) * numg);
  subdomain->Ngp = Ngp;

  subdomain->gmask =
      bfam_subdomain_dgx_gmask_set(numg, N, &subdomain->Np, Ng, Ngp, inDIM);

  subdomain->hadapt = NULL;
  subdomain->padapt = NULL;
  subdomain->q_id = NULL;
  subdomain->lvl = NULL;

  if (inDIM > 0)
  {
    subdomain->K = K;
    subdomain->N = N;
    subdomain->hadapt = bfam_malloc_aligned(K * sizeof(uint8_t));
    subdomain->padapt = bfam_malloc_aligned(K * sizeof(int8_t));
    subdomain->q_id = bfam_malloc_aligned(K * sizeof(bfam_locidx_t));
    subdomain->lvl = bfam_malloc_aligned(K * sizeof(int8_t));
    for (bfam_locidx_t k = 0; k < K; ++k)
    {
      subdomain->hadapt[k] = BFAM_FLAG_SAME;
      subdomain->padapt[k] = (int8_t)N;
      subdomain->q_id[k] = -1;
      subdomain->lvl[k] = -1;
    }

    /*
     * TODO: query dictionary to see if we need to create this, otherwise used
     * stored values
     */
    BFAM_ASSERT(dgx_ops);

    char name[BFAM_BUFSIZ];
    snprintf(name, BFAM_BUFSIZ, "lr_%d", N);
    if (!bfam_dictionary_contains(dgx_ops, name))
    {

      const int Nrp = N + 1;
      bfam_long_real_t *lr =
          bfam_malloc_aligned(Nrp * sizeof(bfam_long_real_t));
      bfam_long_real_t *lw =
          bfam_malloc_aligned(Nrp * sizeof(bfam_long_real_t));
      bfam_long_real_t *lV =
          bfam_malloc_aligned(Nrp * Nrp * sizeof(bfam_long_real_t));
      bfam_long_real_t *lDr =
          bfam_malloc_aligned(Nrp * Nrp * sizeof(bfam_long_real_t));
      bfam_real_t *Dr = bfam_malloc_aligned(Nrp * Nrp * sizeof(bfam_real_t));
      bfam_real_t *r = bfam_malloc_aligned(Nrp * sizeof(bfam_real_t));
      bfam_real_t *w = bfam_malloc_aligned(Nrp * sizeof(bfam_real_t));
      bfam_real_t *wi = bfam_malloc_aligned(Nrp * sizeof(bfam_real_t));

      bfam_jacobi_gauss_lobatto_quadrature(0, 0, N, lr, lw);
      bfam_jacobi_p_vandermonde(0, 0, N, Nrp, lr, lV);
      bfam_jacobi_p_differentiation(0, 0, N, Nrp, lr, lV, lDr);

      /* store the volume stuff */
      for (int n = 0; n < Nrp; ++n)
      {
        r[n] = (bfam_real_t)lr[n];
        w[n] = (bfam_real_t)lw[n];
        wi[n] = (bfam_real_t)(1.0l / lw[n]);
      }
      for (int n = 0; n < Nrp * Nrp; ++n)
      {
        Dr[n] = (bfam_real_t)lDr[n];
      }

      int BFAM_UNUSED_VAR rval = 1;

      snprintf(name, BFAM_BUFSIZ, "lr_%d", N);
      rval = bfam_dictionary_insert_ptr(dgx_ops, name, lr);
      BFAM_ASSERT(rval != 1);

      snprintf(name, BFAM_BUFSIZ, "lw_%d", N);
      rval = bfam_dictionary_insert_ptr(dgx_ops, name, lw);
      BFAM_ASSERT(rval != 1);

      snprintf(name, BFAM_BUFSIZ, "lV_%d", N);
      rval = bfam_dictionary_insert_ptr(dgx_ops, name, lV);
      BFAM_ASSERT(rval != 1);

      snprintf(name, BFAM_BUFSIZ, "lDr_%d", N);
      rval = bfam_dictionary_insert_ptr(dgx_ops, name, lDr);
      BFAM_ASSERT(rval != 1);

      snprintf(name, BFAM_BUFSIZ, "Dr_%d", N);
      rval = bfam_dictionary_insert_ptr(dgx_ops, name, Dr);
      BFAM_ASSERT(rval != 1);

      snprintf(name, BFAM_BUFSIZ, "r_%d", N);
      rval = bfam_dictionary_insert_ptr(dgx_ops, name, r);
      BFAM_ASSERT(rval != 1);

      snprintf(name, BFAM_BUFSIZ, "w_%d", N);
      rval = bfam_dictionary_insert_ptr(dgx_ops, name, w);
      BFAM_ASSERT(rval != 1);

      snprintf(name, BFAM_BUFSIZ, "wi_%d", N);
      rval = bfam_dictionary_insert_ptr(dgx_ops, name, wi);
      BFAM_ASSERT(rval != 1);
    }

    snprintf(name, BFAM_BUFSIZ, "lr_%d", N);
    subdomain->lr = bfam_dictionary_get_value_ptr(dgx_ops, name);
    BFAM_ASSERT(subdomain->lr != NULL);

    snprintf(name, BFAM_BUFSIZ, "lw_%d", N);
    subdomain->lw = bfam_dictionary_get_value_ptr(dgx_ops, name);
    BFAM_ASSERT(subdomain->lw != NULL);

    snprintf(name, BFAM_BUFSIZ, "lV_%d", N);
    subdomain->lV = bfam_dictionary_get_value_ptr(dgx_ops, name);
    BFAM_ASSERT(subdomain->lV != NULL);

    snprintf(name, BFAM_BUFSIZ, "lDr_%d", N);
    subdomain->lDr = bfam_dictionary_get_value_ptr(dgx_ops, name);
    BFAM_ASSERT(subdomain->lDr != NULL);

    snprintf(name, BFAM_BUFSIZ, "Dr_%d", N);
    subdomain->Dr = bfam_dictionary_get_value_ptr(dgx_ops, name);
    BFAM_ASSERT(subdomain->Dr != NULL);

    snprintf(name, BFAM_BUFSIZ, "r_%d", N);
    subdomain->r = bfam_dictionary_get_value_ptr(dgx_ops, name);
    BFAM_ASSERT(subdomain->r != NULL);

    snprintf(name, BFAM_BUFSIZ, "w_%d", N);
    subdomain->w = bfam_dictionary_get_value_ptr(dgx_ops, name);
    BFAM_ASSERT(subdomain->w != NULL);

    snprintf(name, BFAM_BUFSIZ, "wi_%d", N);
    subdomain->wi = bfam_dictionary_get_value_ptr(dgx_ops, name);
    BFAM_ASSERT(subdomain->wi != NULL);
  }
}

/** initializes a dgx subdomain
 *
 * \param [in,out] subdomain pointer to the subdomain to initialize
 * \param [in]     id        unique id number for this subdomain
 * \param [in]     uid       user id number for this subdomain
 * \param [in]     name      name of this subdomain
 * \param [in]     N         polynomial order of elements in each dimension
 * \param [in]     K         number of elements in the subdomain
 * \param [in]     EToQ      Mapping such that \c EToQ[k] gives the pxest
 *                           quadrant number of element \c k.
 * \param [in]     EToE      Mapping such that \c EToE[k*4+f] gives the element
 *                           number connected to face \c f of element \c k.
 * \param [in]     EToF      Mapping such that \c EToF[k*4+f] gives the
 *                           face+orientation number connected to face \c f of
 *                           element \c k.
 * \param [in]     N2N       dictionary of projections and interpolations
 * \param [in]     dgx_ops   dictionary of dgx operators (derivative, etc.)
 * \param [in]     dim       number of (computational) dimensions
 */
static void bfam_subdomain_dgx_init(bfam_subdomain_dgx_t *subdomain,
                                    const bfam_locidx_t id,
                                    const bfam_locidx_t uid, const char *name,
                                    const int N, const bfam_locidx_t K,
                                    const bfam_locidx_t *EToQ,
                                    const bfam_locidx_t *EToE,
                                    const int8_t *EToF, bfam_dictionary_t *N2N,
                                    bfam_dictionary_t *dgx_ops, const int inDIM)
{
  bfam_subdomain_dgx_generic_init(subdomain, id, uid, name, N, K, N2N, dgx_ops,
                                  inDIM);

  const int *Ng = subdomain->Ng;
  const int *Ngp = subdomain->Ngp;
  const int Np = subdomain->Np;

  if (EToQ)
  {
    subdomain->EToQ = bfam_malloc_aligned(K * sizeof(bfam_locidx_t));
    for (bfam_locidx_t k = 0; k < K; k++)
      subdomain->EToQ[k] = EToQ[k];
  }

  if (inDIM > 0)
  {
    /* store the face stuff */
    subdomain->vmapP =
        bfam_malloc_aligned(K * Ngp[0] * Ng[0] * sizeof(bfam_locidx_t));
    subdomain->vmapM =
        bfam_malloc_aligned(K * Ngp[0] * Ng[0] * sizeof(bfam_locidx_t));

    bfam_subdomain_dgx_buildmaps(N, K, Np, Ngp[0], Ng[0], EToE, EToF,
                                 subdomain->gmask, subdomain->vmapP,
                                 subdomain->vmapM, inDIM);
  }
}

/** create a dgx subdomain.
 *
 * \param [in] id      unique id number for this subdomain
 * \param [in] uid     user id number for this subdomain
 * \param [in] name    name of this subdomain
 * \param [in] N       polynomial order of elements in each dimension
 * \param [in] K       number of elements in the subdomain
 * \param [in] EToE    Mapping such that \c EToE[k*4+f] gives the element number
 *                     connected to face \c f of element \c k.
 * \param [in] EToF    Mapping such that \c EToF[k*4+f] gives the
 *                     face+orientation number connected to face \c f of element
 *                     \c k.
 * \param [in] N2N     dictionary of projections and interpolations
 * \param [in] dgx_ops dictionary of dgx operators (derivative, etc.)
 * \param [in] dim     number of (computational) dimensions
 *
 * \return Initialized dgx subdomain
 *
 */
static bfam_subdomain_dgx_t *
bfam_subdomain_dgx_new(const bfam_locidx_t id, const bfam_locidx_t uid,
                       const char *name, const int N, const bfam_locidx_t K,
                       const bfam_locidx_t *EToQ, const bfam_locidx_t *EToE,
                       const int8_t *EToF, bfam_dictionary_t *N2N,
                       bfam_dictionary_t *dgx_ops, const int inDIM)
{
  bfam_subdomain_dgx_t *newSubdomain =
      bfam_malloc(sizeof(bfam_subdomain_dgx_t));

  bfam_subdomain_dgx_init(newSubdomain, id, uid, name, N, K, EToQ, EToE, EToF,
                          N2N, dgx_ops, inDIM);
  return newSubdomain;
}

static void bfam_subdomain_dgx_glue_generic_init(
    bfam_subdomain_dgx_glue_data_t *glue, const bfam_locidx_t rank,
    const bfam_locidx_t id, const bfam_locidx_t id_s,
    bfam_subdomain_dgx_t *sub_m, int inDIM)
{
  bfam_subdomain_glue_init(&glue->base, rank, id, id_s,
                           (bfam_subdomain_t *)sub_m);
  glue->EToEp = NULL;
  glue->EToHp = NULL;
  glue->EToEm = NULL;
  glue->EToFm = NULL;
  glue->EToHm = NULL;
  glue->EToOp = NULL;
  glue->mapOp = NULL;
  glue->num_orient = 0;
  glue->num_interp = 0;
  glue->interpolation = NULL;
  glue->projection = NULL;
  glue->massprojection = NULL;
  glue->exact_mass = NULL;
  glue->base.tags.root = NULL;
}

/** initializes a dg glue subdomain.
 *
 * \param [out]    subdomain    pointer to the subdomain to initialize
 * \param [in]     id           unique id number for this subdomain
 * \param [in]     uid          user id number for this subdomain
 * \param [in]     name         name of this subdomain
 * \param [in]     N_m          minus side polynomial order of elements in each
 *                              dimension
 * \param [in]     N_p          plus  side polynomial order of elements in each
 *                              dimension
 * \param [in]     N_g          polynomial order of each element in this glue
 * \param [in]     rank_m       minus side processor rank
 * \param [in]     rank_p       plus  side processor rank
 * \param [in]     id_m         minus side subdomain id
 * \param [in]     id_p         plus  side subdomain id
 * \param [in]     subdomain_m  minus side subdomain pointer
 * \param [in]     ktok_m       map: element number -> minus side element number
 * \param [in]     K            number of elements in the glue grid
 * \param [in,out] mapping      face mapping (might get sorted)
 * \param [in]     N2N          dictionary of projections and interpolations
 * \param [in]     dgx_ops      dictionary of dgx operators (derivative, etc.)
 * \param [in]     inDIM        dimension of the subdomain
 *
 */
static void bfam_subdomain_dgx_glue_init(
    bfam_subdomain_dgx_t *subdomain, const bfam_locidx_t id,
    const bfam_locidx_t uid, const char *name, const int N_m, const int N_p,
    const int N_g, const bfam_locidx_t rank_m, const bfam_locidx_t rank_p,
    const bfam_locidx_t id_m, const bfam_locidx_t id_p,
    bfam_subdomain_dgx_t *sub_m, bfam_locidx_t *ktok_m, const bfam_locidx_t K,
    bfam_subdomain_face_map_entry_t *mapping, bfam_dictionary_t *N2N,
    bfam_dictionary_t *dgx_ops, const int inDIM)
{
  BFAM_ASSERT(inDIM > 0);

  bfam_subdomain_dgx_generic_init(subdomain, id, uid, name, N_g, K, N2N,
                                  dgx_ops, inDIM);

#ifdef BFAM_DEBUG
  {
    BFAM_LDEBUG("glue mapping: %5s %5s %5s %5s %5s %5s %5s %5s %5s %5s %5s",
                "np", "nk", "nf", "nh", "gi", "i", "k", "f", "h", "o", "id");
    for (bfam_locidx_t k = 0; k < K; ++k)
    {
      BFAM_LDEBUG("glue mapping: %5jd %5jd %5jd %5jd %5jd %5jd %5jd %5jd %5jd "
                  "%5jd %5jd",
                  (intmax_t)mapping[k].np, (intmax_t)mapping[k].nk,
                  (intmax_t)mapping[k].nf, (intmax_t)mapping[k].nh,
                  (intmax_t)mapping[k].gi, (intmax_t)mapping[k].i,
                  (intmax_t)mapping[k].k, (intmax_t)mapping[k].f,
                  (intmax_t)mapping[k].h, (intmax_t)mapping[k].o,
                  (intmax_t)mapping[k].id);
    }
  }
#endif

  BFAM_ABORT_IF(N_g < N_m || N_g < N_p,
                "glue space must currently be a higher order space"
                " N_g = %d, N_m = %d, N_p = %d",
                N_g, N_m, N_p);

  bfam_subdomain_add_tag(&subdomain->base, "_subdomain_dgx");

  BFAM_ASSERT(subdomain->base.glue_m == NULL);

  subdomain->base.glue_m = bfam_malloc(sizeof(bfam_subdomain_dgx_glue_data_t));
  bfam_subdomain_dgx_glue_data_t *glue_m =
      (bfam_subdomain_dgx_glue_data_t *)subdomain->base.glue_m;
  bfam_subdomain_dgx_glue_generic_init(
      glue_m, rank_m, id_m, (bfam_locidx_t)(imaxabs(id_m) - 1), sub_m, inDIM);

  subdomain->base.glue_p = bfam_malloc(sizeof(bfam_subdomain_dgx_glue_data_t));
  bfam_subdomain_dgx_glue_data_t *glue_p =
      (bfam_subdomain_dgx_glue_data_t *)subdomain->base.glue_p;
  bfam_subdomain_dgx_glue_generic_init(
      glue_p, rank_p, id_p, (bfam_locidx_t)(imaxabs(id_p) - 1), NULL, inDIM);

  const int num_interp = 3;
  glue_m->num_interp = num_interp;

  const int N = subdomain->N;
  const int Nrp = N + 1;

  bfam_subdomain_dgx_interpolator_t *proj_m2g =
      bfam_subdomain_dgx_get_interpolator(N2N, N_m, N, inDIM);
  bfam_subdomain_dgx_interpolator_t *proj_g2m =
      bfam_subdomain_dgx_get_interpolator(N2N, N, N_m, inDIM);

  BFAM_ASSERT(proj_m2g);
  BFAM_ASSERT(proj_g2m);

  glue_m->interpolation =
      bfam_malloc_aligned(glue_m->num_interp * sizeof(bfam_real_t *));
  glue_m->interpolation[0] = proj_m2g->prj[0];
  glue_m->interpolation[1] = proj_m2g->prj[3];
  glue_m->interpolation[2] = proj_m2g->prj[4];

  glue_m->projection =
      bfam_malloc_aligned(glue_m->num_interp * sizeof(bfam_real_t *));
  glue_m->projection[0] = proj_g2m->prj[0];
  glue_m->projection[1] = proj_g2m->prj[1];
  glue_m->projection[2] = proj_g2m->prj[2];

  glue_m->massprojection =
      bfam_malloc_aligned(glue_m->num_interp * sizeof(bfam_real_t *));
  glue_m->massprojection[0] = proj_g2m->mass_prj[0];
  glue_m->massprojection[1] = proj_g2m->mass_prj[1];
  glue_m->massprojection[2] = proj_g2m->mass_prj[2];

  glue_p->same_order = (N_p == N) && (N_m == N);
  glue_p->EToEp = bfam_malloc_aligned(K * sizeof(bfam_locidx_t));
  glue_p->EToHp = bfam_malloc_aligned(K * sizeof(int8_t));
  glue_m->EToEm = bfam_malloc_aligned(K * sizeof(bfam_locidx_t));
  glue_m->EToFm = bfam_malloc_aligned(K * sizeof(int8_t));
  glue_m->EToHm = bfam_malloc_aligned(K * sizeof(int8_t));
  glue_p->EToOp = bfam_malloc_aligned(K * sizeof(int8_t));
  if (inDIM == 1)
    glue_p->num_orient = 2;
  else if (inDIM == 2)
    glue_p->num_orient = 8;
  else
    BFAM_ABORT("Cannot handle dim = %d", inDIM);

  glue_p->mapOp =
      bfam_malloc_aligned(glue_p->num_orient * sizeof(bfam_locidx_t *));
  if (inDIM == 1)
  {
    for (int n = 0; n < glue_p->num_orient; n++)
      glue_p->mapOp[n] = bfam_malloc_aligned(Nrp * sizeof(bfam_locidx_t));
    for (int n = 0; n < Nrp; n++)
    {
      glue_p->mapOp[0][n] = n;
      glue_p->mapOp[1][n] = Nrp - (n + 1);
    }
  }
  else if (inDIM == 2)
  {
    for (int n = 0; n < glue_p->num_orient; n++)
      glue_p->mapOp[n] = bfam_malloc_aligned(Nrp * Nrp * sizeof(bfam_locidx_t));
    for (int j = 0; j < Nrp; j++)
      for (int i = 0; i < Nrp; i++)
      {
        int ir = Nrp - (i + 1);
        int jr = Nrp - (j + 1);

        glue_p->mapOp[0][i + j * Nrp] = i + j * Nrp;
        glue_p->mapOp[1][i + j * Nrp] = j + i * Nrp;
        glue_p->mapOp[2][i + j * Nrp] = ir + j * Nrp;
        glue_p->mapOp[3][i + j * Nrp] = j + ir * Nrp;
        glue_p->mapOp[4][i + j * Nrp] = jr + i * Nrp;
        glue_p->mapOp[5][i + j * Nrp] = i + jr * Nrp;
        glue_p->mapOp[6][i + j * Nrp] = jr + ir * Nrp;
        glue_p->mapOp[7][i + j * Nrp] = ir + jr * Nrp;
      }
  }
  else
    BFAM_ABORT("Cannot handle dim = %d", inDIM);

  qsort(mapping, K, sizeof(bfam_subdomain_face_map_entry_t),
        bfam_subdomain_face_send_cmp);

  for (bfam_locidx_t k = 0; k < K; ++k)
    mapping[k].i = k;

  qsort(mapping, K, sizeof(bfam_subdomain_face_map_entry_t),
        bfam_subdomain_face_recv_cmp);

  for (bfam_locidx_t k = 0; k < K; ++k)
  {
    glue_m->EToEm[k] = ktok_m[mapping[k].k];
    glue_m->EToFm[k] = mapping[k].f;
    glue_m->EToHm[k] = mapping[k].h;
    glue_p->EToOp[k] = mapping[k].o;

    glue_p->EToEp[k] = mapping[k].i;
    glue_p->EToHp[k] = mapping[k].nh;
  }

#ifdef BFAM_DEBUG
  for (bfam_locidx_t k = 0; k < K; ++k)
    BFAM_ASSERT(mapping[k].s == glue_m->base.id_s &&
                mapping[k].ns == glue_p->base.id_s);
  BFAM_LDEBUG("    k EToEp EToHp EToEm EToFm EToHm EToOp");
  for (bfam_locidx_t k = 0; k < K; ++k)
  {
    BFAM_LDEBUG("%5d %5d %5d %5d %5d %5d %5d", k, glue_p->EToEp[k],
                glue_p->EToHp[k], glue_m->EToEm[k], glue_m->EToFm[k],
                glue_m->EToHm[k], glue_p->EToOp[k]);
  }
#endif

  for (int i = 0; i < num_interp; ++i)
  {
    // bfam_free_aligned(lr[i]);
    // bfam_free_aligned(interpolation[i]);
    // bfam_free_aligned(massprojection[i]);
    // bfam_free_aligned(projection[i]);
  }

  // bfam_free_aligned(lV);
  // bfam_free_aligned(mass);

  // bfam_free_aligned(lr);
  // bfam_free_aligned(lw);

  // bfam_free_aligned(interpolation);
  // bfam_free_aligned(massprojection);
  // bfam_free_aligned(projection);

  // bfam_free_aligned(sub_m_lr);
  // bfam_free_aligned(sub_m_lw);
  // bfam_free_aligned(sub_m_V);
  // bfam_free_aligned(Vt_MP);
}

/** create a dgx glue subdomain.
 *
 * \param [in]     id          unique id number for this subdomain
 * \param [in]     uid         user id number for this subdomain
 * \param [in]     name        name of this subdomain
 * \param [in]     N_m         minus side polynomial order of elements in each
 *                             dimension
 * \param [in]     N_p         plus  side polynomial order of elements in each
 *                             dimension
 * \param [in]     N_g         polynomial order of each element in this glue
 * \param [in]     rank_m      minus side processor rank
 * \param [in]     rank_p      plus  side processor rank
 * \param [in]     id_m        minus side subdomain id
 * \param [in]     id_p        plus  side subdomain id
 * \param [in]     subdomain_m minus side subdomain pointer
 * \param [in]     ktok_m      map: element number -> minus side element number
 * \param [in]     K           number of elements in the glue grid
 * \param [in,out] mapping     face mapping (might get sorted)
 * \param [in]     N2N         dictionary of projections and interpolations
 * \param [in]     dgx_ops     dictionary of dgx operators (derivative, etc.)
 * \param [in]     inDIM       dimension of the subdomain
 *
 * \return Initialized dg quad glue subdomain
 *
 */
static bfam_subdomain_dgx_t *bfam_subdomain_dgx_glue_new(
    const bfam_locidx_t id, const bfam_locidx_t uid, const char *name,
    const int N_m, const int N_p, const int N_g, const bfam_locidx_t rank_m,
    const bfam_locidx_t rank_p, const bfam_locidx_t id_m,
    const bfam_locidx_t id_p, bfam_subdomain_dgx_t *sub_m,
    bfam_locidx_t *ktok_m, const bfam_locidx_t K,
    bfam_subdomain_face_map_entry_t *mapping, bfam_dictionary_t *N2N,
    bfam_dictionary_t *dgx_ops, const int inDIM)
{
  bfam_subdomain_dgx_t *newSubdomain =
      bfam_malloc(sizeof(bfam_subdomain_dgx_t));

  bfam_subdomain_dgx_glue_init(newSubdomain, id, uid, name, N_m, N_p, N_g,
                               rank_m, rank_p, id_m, id_p, sub_m, ktok_m, K,
                               mapping, N2N, dgx_ops, inDIM);
  return newSubdomain;
}

void bfam_domain_pxest_split_dgx_subdomains(
    bfam_domain_pxest_t *domain, bfam_locidx_t num_subdomains,
    bfam_locidx_t *subdomainID, bfam_locidx_t *roots, int *N,
    bfam_locidx_t *glueID, bfam_glue_order_t glue_order, void *go_user_args)
{
  BFAM_ROOT_LDEBUG("Begin splitting p4est domain into subdomains.");
  const int HF = P4EST_HALF * P4EST_FACES;

  p4est_t *pxest = domain->pxest;
  p4est_ghost_t *ghost = p4est_ghost_new(pxest, BFAM_PXEST_CONNECT);
  p4est_mesh_t *mesh = p4est_mesh_new(pxest, ghost, BFAM_PXEST_CONNECT);
  p4est_nodes_t *nodes = p4est_nodes_new(pxest, NULL);

  p4est_locidx_t *subK = bfam_calloc(num_subdomains, sizeof(p4est_locidx_t));
  p4est_locidx_t *subk = bfam_calloc(num_subdomains, sizeof(p4est_locidx_t));
  char **name = bfam_malloc(num_subdomains * sizeof(char *));
  bfam_locidx_t **EToE = bfam_malloc(num_subdomains * sizeof(bfam_locidx_t *));
  bfam_locidx_t **EToQ = bfam_malloc(num_subdomains * sizeof(bfam_locidx_t *));
  int8_t **EToF = bfam_malloc(num_subdomains * sizeof(int8_t *));

  bfam_locidx_t *sub_to_actual_sub_id =
      bfam_malloc(num_subdomains * sizeof(bfam_locidx_t));

  bfam_locidx_t *ktosubk =
      bfam_malloc(mesh->local_num_quadrants * sizeof(bfam_locidx_t));

  bfam_locidx_t *ghostSubdomainID =
      bfam_malloc(mesh->ghost_num_quadrants * sizeof(bfam_locidx_t));

  bfam_locidx_t *ghostN = bfam_malloc(mesh->ghost_num_quadrants * sizeof(int));

  bfam_locidx_t numParallelFaces = bfam_domain_pxest_num_parallel_faces(mesh);

  bfam_subdomain_face_map_entry_t *pfmapping = bfam_malloc_aligned(
      numParallelFaces * sizeof(bfam_subdomain_face_map_entry_t));

  bfam_domain_pxest_parallel_face_mapping(ghost, mesh, glueID, numParallelFaces,
                                          pfmapping);

  bfam_locidx_t numNeighbors = bfam_domain_pxest_parallel_face_num_neighbors(
      numParallelFaces, pfmapping);

  bfam_locidx_t *numNeighborFaces =
      bfam_malloc_aligned(numNeighbors * sizeof(bfam_locidx_t));

  bfam_locidx_t *neighborRank =
      bfam_malloc_aligned(numNeighbors * sizeof(bfam_locidx_t));

  bfam_domain_pxest_num_neighbor_faces(numParallelFaces, pfmapping,
                                       numNeighbors, numNeighborFaces,
                                       neighborRank);

#ifdef BFAM_DEBUG
  bfam_domain_pxest_parallel_face_mapping_check(
      pxest->mpicomm, numParallelFaces, pfmapping, numNeighbors,
      numNeighborFaces, neighborRank);
#endif

  bfam_domain_pxest_fill_ghost_subdomain_ids(
      pxest->mpicomm, numParallelFaces, pfmapping, numNeighbors,
      numNeighborFaces, neighborRank, subdomainID, N, ghostSubdomainID, ghostN);

  bfam_locidx_t numInterSubdomainFaces =
      bfam_domain_pxest_num_inter_subdomain_faces(mesh, subdomainID, glueID);

  BFAM_LDEBUG("numInterSubdomainFaces = %jd", (intmax_t)numInterSubdomainFaces);

  bfam_subdomain_face_map_entry_t *ifmapping = bfam_malloc_aligned(
      numInterSubdomainFaces * sizeof(bfam_subdomain_face_map_entry_t));

  int rank;
  BFAM_MPI_CHECK(MPI_Comm_rank(pxest->mpicomm, &rank));

  bfam_domain_pxest_inter_subdomain_face_mapping(
      rank, mesh, subdomainID, glueID, numInterSubdomainFaces, ifmapping);

  bfam_locidx_t numBoundaryFaces = bfam_domain_pxest_num_boundary_faces(mesh);

  BFAM_LDEBUG("numBoundaryFaces = %jd", (intmax_t)numBoundaryFaces);

  bfam_subdomain_face_map_entry_t *bfmapping = bfam_malloc_aligned(
      numBoundaryFaces * sizeof(bfam_subdomain_face_map_entry_t));

  bfam_domain_pxest_boundary_subdomain_face_mapping(
      mesh, subdomainID, glueID, numBoundaryFaces, bfmapping);

  /*
   * Count the number of elements in each new subdomain
   */
  for (p4est_locidx_t k = 0; k < pxest->local_num_quadrants; ++k)
  {
    bfam_locidx_t id = (bfam_locidx_t)subdomainID[k];

    BFAM_ABORT_IF(id < 0 || id >= num_subdomains, "Bad Subdomain id: %jd",
                  (intmax_t)id);

    ++subK[id];
  }

  for (bfam_locidx_t id = 0; id < num_subdomains; ++id)
  {
    name[id] = bfam_malloc(BFAM_BUFSIZ * sizeof(char));
    snprintf(name[id], BFAM_BUFSIZ, "dgx_dim_%1d_%05jd", (int)DIM,
             (intmax_t)id);

    EToQ[id] = bfam_malloc(subK[id] * sizeof(bfam_locidx_t));
    EToE[id] = bfam_malloc(subK[id] * P4EST_FACES * sizeof(bfam_locidx_t));
    EToF[id] = bfam_malloc(subK[id] * P4EST_FACES * sizeof(int8_t));
  }

  const p4est_locidx_t K = mesh->local_num_quadrants;

  BFAM_ASSERT(K == pxest->local_num_quadrants);

  for (bfam_locidx_t id = 0; id < num_subdomains; ++id)
  {
    subk[id] = 0;
  }
  for (p4est_locidx_t k = 0; k < K; ++k)
  {
    const bfam_locidx_t idk = subdomainID[k];
    ktosubk[k] = subk[idk];
    ++subk[idk];
  }

  /*
   * Here we are decoding the p4est_mesh_t structure.  See p4est_mesh.h
   * for more details on how the data is stored.
   */

  /*
   * First build up the volume grids
   */
  for (bfam_locidx_t id = 0; id < num_subdomains; ++id)
  {
    subk[id] = 0;
  }
  for (p4est_locidx_t k = 0; k < K; ++k)
  {
    const bfam_locidx_t idk = subdomainID[k];

    EToQ[idk][subk[idk]] = k;

    for (int f = 0; f < P4EST_FACES; ++f)
    {
      const p4est_locidx_t ck = mesh->quad_to_quad[P4EST_FACES * k + f];
      const int cf = mesh->quad_to_face[P4EST_FACES * k + f];
      const bfam_locidx_t glueid = (glueID) ? glueID[P4EST_FACES * k + f] : -1;

      int nk = k;
      int nf = f;

      if (glueid < 0 && cf >= 0 && cf < HF && ck < K && idk == subdomainID[ck])
      {
        /*
         * Neighbor is local to the subdomain and is the same size and does
         * not have a glue grid inbetween.
         */
        nk = ck;
        nf = cf;
      }

      EToE[idk][P4EST_FACES * subk[idk] + f] = ktosubk[nk];
      EToF[idk][P4EST_FACES * subk[idk] + f] = (int8_t)nf;
    }

    ++subk[idk];
  }

  bfam_subdomain_dgx_t **subdomains =
      bfam_malloc(num_subdomains * sizeof(bfam_subdomain_dgx_t **));
  for (bfam_locidx_t id = 0; id < num_subdomains; ++id)
  {
    if (roots)
      subdomains[id] = bfam_subdomain_dgx_new(
          id, roots[id], name[id], N[id], subK[id], EToQ[id], EToE[id],
          EToF[id], domain->N2N, domain->dgx_ops, DIM);
    else
      subdomains[id] = bfam_subdomain_dgx_new(
          id, -1, name[id], N[id], subK[id], EToQ[id], EToE[id], EToF[id],
          domain->N2N, domain->dgx_ops, DIM);

    bfam_subdomain_add_tag((bfam_subdomain_t *)subdomains[id], "_volume");
    char root_id_tag[BFAM_BUFSIZ];
    if (roots)
      snprintf(root_id_tag, BFAM_BUFSIZ, "_volume_id_%jd", (intmax_t)roots[id]);
    else
      snprintf(root_id_tag, BFAM_BUFSIZ, "_volume_id_0");
    bfam_subdomain_add_tag((bfam_subdomain_t *)subdomains[id], root_id_tag);

    sub_to_actual_sub_id[id] = bfam_domain_add_subdomain(
        (bfam_domain_t *)domain, (bfam_subdomain_t *)subdomains[id]);
  }

  /*
   * Sort the local mapping
   */
  qsort(ifmapping, numInterSubdomainFaces,
        sizeof(bfam_subdomain_face_map_entry_t), bfam_subdomain_face_recv_cmp);

  /*
   * Setup the local glue grids
   */
  bfam_locidx_t numGlue = 0;

  for (bfam_locidx_t ifk = 0; ifk < numInterSubdomainFaces;)
  {
    /*
     * Count the number of element in the glue grid
     */
    bfam_locidx_t Kglue = 1;
    while (ifk + Kglue < numInterSubdomainFaces &&
           ifmapping[ifk + Kglue].np == ifmapping[ifk + Kglue - 1].np &&
           ifmapping[ifk + Kglue].ns == ifmapping[ifk + Kglue - 1].ns &&
           ifmapping[ifk + Kglue].s == ifmapping[ifk + Kglue - 1].s &&
           ifmapping[ifk + Kglue].id == ifmapping[ifk + Kglue - 1].id)
      ++Kglue;

    const bfam_locidx_t id_m = ifmapping[ifk].s;
    const bfam_locidx_t id_p = ifmapping[ifk].ns;
    const bfam_locidx_t glueid = ifmapping[ifk].id;

    int repeat = (id_m == id_p);
    repeat = 0;

    for (int r = 0; r <= repeat; ++r)
    {
      const bfam_locidx_t id = num_subdomains + numGlue;

      const bfam_locidx_t rank_m = rank;
      const bfam_locidx_t rank_p = rank;

      char glueName[BFAM_BUFSIZ];
      snprintf(glueName, BFAM_BUFSIZ,
               "dgx_dim_%1d_glue_%05jd_%05jd_%05jd_%05jd", (int)BDIM,
               (intmax_t)id, (intmax_t)id_m, (intmax_t)id_p, (intmax_t)glueid);

      /*
       * For subdomains that connect to themselves we need to distinguish
       * between them based on id.  So we have decided to use a minus sign
       * to distinguish between the two different glue grids.
       */
      // bfam_locidx_t sign_m = 1;
      bfam_locidx_t sign_p = 1;
      if (id_m == id_p)
      {
        // sign_m = (r) ? -1 :  1;
        sign_p = (r) ? 1 : -1;
      }

      int N_g = BFAM_MAX(N[id_m], N[id_p]);
      if (glue_order && glueid >= 0)
        N_g = glue_order(N[id_m], N[id_p], glueid, go_user_args);
      bfam_subdomain_dgx_t *glue = bfam_subdomain_dgx_glue_new(
          id, glueid, glueName, N[id_m], N[id_p], N_g, rank_m, rank_p,
          sign_p * (id_m + 1), sign_p * (id_p + 1), subdomains[id_m], ktosubk,
          Kglue, ifmapping + ifk, domain->N2N, domain->dgx_ops, DIM - 1);

      bfam_subdomain_add_tag((bfam_subdomain_t *)glue, "_glue");
      bfam_subdomain_add_tag((bfam_subdomain_t *)glue, "_glue_local");
      char glue_num_tag[BFAM_BUFSIZ];
      snprintf(glue_num_tag, BFAM_BUFSIZ, "_glue_%jd_%jd",
               (intmax_t)BFAM_MIN(id_m, id_p), (intmax_t)BFAM_MAX(id_m, id_p));
      bfam_subdomain_add_tag((bfam_subdomain_t *)glue, glue_num_tag);

      char glue_id_tag[BFAM_BUFSIZ];
      snprintf(glue_id_tag, BFAM_BUFSIZ, "_glue_id_%jd", (intmax_t)glueid);
      bfam_subdomain_add_tag((bfam_subdomain_t *)glue, glue_id_tag);

      bfam_domain_add_subdomain((bfam_domain_t *)domain,
                                (bfam_subdomain_t *)glue);

      numGlue += 1;
    }
    ifk += Kglue;
  }

  /*
   * Sort the boundary mapping
   */
  qsort(bfmapping, numBoundaryFaces, sizeof(bfam_subdomain_face_map_entry_t),
        bfam_subdomain_face_recv_cmp);

  /*
   * Setup the boundary glue grids
   */
  for (bfam_locidx_t bfk = 0; bfk < numBoundaryFaces;)
  {
    /*
     * Count the number of elements in the glue grid
     */
    bfam_locidx_t Kglue = 1;
    while (bfk + Kglue < numBoundaryFaces &&
           bfmapping[bfk + Kglue].np == bfmapping[bfk + Kglue - 1].np &&
           bfmapping[bfk + Kglue].ns == bfmapping[bfk + Kglue - 1].ns &&
           bfmapping[bfk + Kglue].s == bfmapping[bfk + Kglue - 1].s &&
           bfmapping[bfk + Kglue].id == bfmapping[bfk + Kglue - 1].id)
      ++Kglue;

    const bfam_locidx_t id = num_subdomains + numGlue;

    const bfam_locidx_t id_m = bfmapping[bfk].s;
    const bfam_locidx_t id_p = bfmapping[bfk].ns;

    const bfam_locidx_t rank_m = rank;
    const bfam_locidx_t rank_p = bfmapping[bfk].np;

    const bfam_locidx_t glueid = bfmapping[bfk].id;

    char glueName[BFAM_BUFSIZ];
    snprintf(glueName, BFAM_BUFSIZ,
             "dgx_dim_%1d_glue_b_%05jd_%05jd_%05jd_%05jd", (int)BDIM,
             (intmax_t)id, (intmax_t)id_m, (intmax_t)id_p, (intmax_t)glueid);

    int N_g = BFAM_MAX(N[id_m], N[id_p]);
    if (glue_order && glueid >= 0)
      N_g = glue_order(N[id_m], N[id_p], glueid, go_user_args);
    bfam_subdomain_dgx_t *glue = bfam_subdomain_dgx_glue_new(
        id, glueid, glueName, N[id_m], N[id_p], N_g, rank_m, rank_p, id_m + 1,
        id_p + 1, subdomains[id_m], ktosubk, Kglue, bfmapping + bfk,
        domain->N2N, domain->dgx_ops, DIM - 1);

    bfam_subdomain_add_tag((bfam_subdomain_t *)glue, "_glue");
    bfam_subdomain_add_tag((bfam_subdomain_t *)glue, "_glue_boundary");

    char glue_id_tag[BFAM_BUFSIZ];
    snprintf(glue_id_tag, BFAM_BUFSIZ, "_glue_id_%jd", (intmax_t)glueid);
    bfam_subdomain_add_tag((bfam_subdomain_t *)glue, glue_id_tag);

    bfam_domain_add_subdomain((bfam_domain_t *)domain,
                              (bfam_subdomain_t *)glue);

    numGlue += 1;
    bfk += Kglue;
  }

  /*
   * Sort the parallel mapping
   */
  qsort(pfmapping, numParallelFaces, sizeof(bfam_subdomain_face_map_entry_t),
        bfam_subdomain_face_recv_cmp);

  /*
   * Setup the parallel glue grids
   */
  for (bfam_locidx_t pfk = 0; pfk < numParallelFaces;)
  {
    /*
     * Count the number of elements in the glue grid
     */
    bfam_locidx_t Kglue = 1;
    while (pfk + Kglue < numParallelFaces &&
           pfmapping[pfk + Kglue].np == pfmapping[pfk + Kglue - 1].np &&
           pfmapping[pfk + Kglue].ns == pfmapping[pfk + Kglue - 1].ns &&
           pfmapping[pfk + Kglue].s == pfmapping[pfk + Kglue - 1].s &&
           pfmapping[pfk + Kglue].id == pfmapping[pfk + Kglue - 1].id)
      ++Kglue;

    const bfam_locidx_t id = num_subdomains + numGlue;

    const bfam_locidx_t id_m = pfmapping[pfk].s;
    const bfam_locidx_t id_p = pfmapping[pfk].ns;
    const bfam_locidx_t gid_p = pfmapping[pfk].gi;
    const bfam_locidx_t glueid = pfmapping[pfk].id;

    const bfam_locidx_t rank_m = rank;
    const bfam_locidx_t rank_p = pfmapping[pfk].np;

    char glueName[BFAM_BUFSIZ];
    snprintf(glueName, BFAM_BUFSIZ,
             "dgx_dim_%1d_glue_p_%05jd_%05jd_%05jd_%05jd", (int)BDIM,
             (intmax_t)id, (intmax_t)id_m, (intmax_t)id_p, (intmax_t)glueid);

    int N_g = BFAM_MAX(N[id_m], ghostN[gid_p]);
    if (glue_order && glueid >= 0)
      N_g = glue_order(N[id_m], ghostN[gid_p], glueid, go_user_args);
    bfam_subdomain_dgx_t *glue = bfam_subdomain_dgx_glue_new(
        id, glueid, glueName, N[id_m], ghostN[gid_p], N_g, rank_m, rank_p,
        id_m + 1, id_p + 1, subdomains[id_m], ktosubk, Kglue, pfmapping + pfk,
        domain->N2N, domain->dgx_ops, DIM - 1);

    bfam_subdomain_add_tag((bfam_subdomain_t *)glue, "_glue");
    bfam_subdomain_add_tag((bfam_subdomain_t *)glue, "_glue_parallel");
    char glue_num_tag[BFAM_BUFSIZ];
    snprintf(glue_num_tag, BFAM_BUFSIZ, "_glue_%jd_%jd",
             (intmax_t)BFAM_MIN(id_m, id_p), (intmax_t)BFAM_MAX(id_m, id_p));
    bfam_subdomain_add_tag((bfam_subdomain_t *)glue, glue_num_tag);

    char glue_id_tag[BFAM_BUFSIZ];
    snprintf(glue_id_tag, BFAM_BUFSIZ, "_glue_id_%jd", (intmax_t)glueid);
    bfam_subdomain_add_tag((bfam_subdomain_t *)glue, glue_id_tag);

    bfam_domain_add_subdomain((bfam_domain_t *)domain,
                              (bfam_subdomain_t *)glue);

    numGlue += 1;
    pfk += Kglue;
  }

  /*
   * Fill the quadrant user data
   */
  {
    p4est_topidx_t t;
    p4est_locidx_t k;
    for (t = pxest->first_local_tree, k = 0; t <= pxest->last_local_tree; ++t)
    {
      p4est_tree_t *tree = p4est_tree_array_index(pxest->trees, t);
      sc_array_t *quadrants = &tree->quadrants;
      size_t num_quads = quadrants->elem_count;

      for (size_t zz = 0; zz < num_quads; ++zz, ++k)
      {
        BFAM_ASSERT(k < K);
        bfam_pxest_user_data_t *ud;
        p4est_quadrant_t *quad;

        const bfam_locidx_t new_subd_id = sub_to_actual_sub_id[subdomainID[k]];
        const bfam_locidx_t new_elem_id = ktosubk[k];

        quad = p4est_quadrant_array_index(quadrants, zz);
        ud = quad->p.user_data;

        bfam_subdomain_dgx_t *new_subdomain =
            (bfam_subdomain_dgx_t *)bfam_domain_get_subdomain_by_num(
                (bfam_domain_t *)domain, new_subd_id);

        new_subdomain->hadapt[new_elem_id] = ud->flags;
        if (ud->N > 0)
          new_subdomain->padapt[new_elem_id] = ud->N;

        new_subdomain->q_id[new_elem_id] = k;
        new_subdomain->lvl[new_elem_id] = quad->level;

        ud->N = (int8_t)N[subdomainID[k]];
        ud->Nold = (int8_t)N[subdomainID[k]];
        ud->subd_id = new_subd_id;
        ud->elem_id = new_elem_id;

        ud->root_id =
            (roots) ? roots[subdomainID[k]] : BFAM_DEFAULT_SUBDOMAIN_ROOT;
        for (int f = 0; f < P4EST_FACES; ++f)
        {
          ud->glue_id[f] = (glueID) ? glueID[P4EST_FACES * k + f] : -1;
        }
      }
    }
  }

  /*
   * Start Cleanup
   */

  bfam_free(subdomains);

  for (bfam_locidx_t id = 0; id < num_subdomains; ++id)
  {
    bfam_free(name[id]);
    bfam_free(EToE[id]);
    bfam_free(EToQ[id]);
    bfam_free(EToF[id]);
  }

  bfam_free(name);
  bfam_free(EToE);
  bfam_free(EToQ);
  bfam_free(EToF);

  bfam_free(subK);
  bfam_free(subk);

  bfam_free(ktosubk);
  bfam_free(sub_to_actual_sub_id);

  bfam_free_aligned(bfmapping);
  bfam_free_aligned(ifmapping);

  bfam_free_aligned(neighborRank);
  bfam_free_aligned(numNeighborFaces);
  bfam_free_aligned(pfmapping);
  bfam_free(ghostN);
  bfam_free(ghostSubdomainID);

  p4est_nodes_destroy(nodes);
  p4est_mesh_destroy(mesh);
  p4est_ghost_destroy(ghost);

  BFAM_ROOT_LDEBUG("End splitting pxest domain into subdomains.");
  bfam_domain_pxest_dgx_print_stats(domain);
}

// }}}

// {{{ vtk

#include <sc.h>
#include <sc_io.h>

#define BFAM_VTK_VTU_FORMAT "%s_%05d.vtu"
#define BFAM_VTK_VTS_FORMAT "%s_%s_%s.vts"
#define BFAM_VTK_PVTS_FORMAT "%s_%s.pvts"

static void bfam_vtk_write_file_pvtu(int size, const char *directory,
                                     const char *prefix, const char **scalars,
                                     const char **vectors,
                                     const char **components, int binary,
                                     int compress)
{
  const int endian = bfam_endian();

  char filename[BFAM_BUFSIZ];
  snprintf(filename, BFAM_BUFSIZ, "%s.pvtu", prefix);

  const char *format;
  if (binary)
    format = "binary";
  else
    format = "ascii";

  BFAM_VERBOSE("Writing file: '%s'", filename);
  FILE *file = fopen(filename, "w");

  if (file == NULL)
  {
    BFAM_LERROR("Could not open %s for output!\n", filename);
    return;
  }

  fprintf(file, "<?xml version=\"1.0\"?>\n");
  fprintf(file, "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\"");

  if (binary && compress)
    fprintf(file, " compressor=\"vtkZLibDataCompressor\"");

  if (endian == BFAM_BIG_ENDIAN)
    fprintf(file, " byte_order=\"BigEndian\">\n");
  else
    fprintf(file, " byte_order=\"LittleEndian\">\n");

  fprintf(file, "  <PUnstructuredGrid GhostLevel=\"0\">\n");
  fprintf(file, "    <PPoints>\n");
  fprintf(file, "      <PDataArray type=\"%s\" Name=\"Position\""
                " NumberOfComponents=\"3\" format=\"%s\"/>\n",
          BFAM_REAL_VTK, format);
  fprintf(file, "    </PPoints>\n");

  fprintf(file,
          "    <PCellData Scalars=\"time,mpirank,subdomain_id,root_id\">\n");
  fprintf(file, "      <PDataArray type=\"%s\" Name=\"time\""
                " format=\"%s\"/>\n",
          BFAM_REAL_VTK, format);
  fprintf(file, "      <PDataArray type=\"%s\" Name=\"mpirank\""
                " format=\"%s\"/>\n",
          BFAM_LOCIDX_VTK, format);
  fprintf(file, "      <PDataArray type=\"%s\" Name=\"subdomain_id\""
                " format=\"%s\"/>\n",
          BFAM_LOCIDX_VTK, format);
  fprintf(file, "      <PDataArray type=\"%s\" Name=\"root_id\""
                " format=\"%s\"/>\n",
          BFAM_LOCIDX_VTK, format);
  fprintf(file, "    </PCellData>\n");

  char pointscalars[BFAM_BUFSIZ];
  bfam_util_strcsl(pointscalars, scalars);

  char pointvectors[BFAM_BUFSIZ];
  bfam_util_strcsl(pointvectors, vectors);

  fprintf(file, "    <PPointData Scalars=\"%s\" Vectors=\"%s\">\n",
          pointscalars, pointvectors);

  if (scalars)
    for (size_t s = 0; scalars[s]; ++s)
      fprintf(file, "      <PDataArray type=\"%s\" Name=\"%s\""
                    " format=\"%s\"/>\n",
              BFAM_REAL_VTK, scalars[s], format);

  if (vectors)
    for (size_t v = 0; vectors[v]; ++v)
      fprintf(file, "      <PDataArray type=\"%s\" Name=\"%s\""
                    " NumberOfComponents=\"3\" format=\"%s\"/>\n",
              BFAM_REAL_VTK, vectors[v], format);

  fprintf(file, "    </PPointData>\n");

  for (int s = 0; s < size; ++s)
  {
    if (directory)
      fprintf(file, "    <Piece Source=\"%s/" BFAM_VTK_VTU_FORMAT "\"/>\n",
              directory, prefix, s);
    else
      fprintf(file, "    <Piece Source=\"" BFAM_VTK_VTU_FORMAT "\"/>\n", prefix,
              s);
  }
  fprintf(file, "  </PUnstructuredGrid>\n");
  fprintf(file, "</VTKFile>\n");

  if (ferror(file))
  {
    BFAM_LERROR("Error writing to %s\n", filename);
  }

  if (fclose(file))
  {
    BFAM_LERROR("Error closing %s\n", filename);
  }
}

/** Utility function for writing an empty vtu file
 *
 * \param [out] file            stream to write the vector to.
 * \param [in]  name            name of the vector.
 * \param [in]  writeBinary     boolean indicating if the data should be written
 *                              in binary.
 */
static void bfam_vtk_write_vtu_empty(FILE *file, int writeBinary)
{
  const char *format;

  if (writeBinary)
    format = "binary";
  else
    format = "ascii";

  fprintf(file, "    <Piece NumberOfPoints=\"%jd\" NumberOfCells=\"%jd\">\n",
          (intmax_t)0, (intmax_t)0);

  /*
   * Points
   */
  fprintf(file, "      <Points>\n");
  fprintf(file, "        <DataArray type=\"%s\" Name=\"%s\""
                " NumberOfComponents=\"3\" format=\"%s\">\n",
          BFAM_REAL_VTK, "Position", format);
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "      </Points>\n");

  /*
   * Cells
   */
  fprintf(file, "      <Cells>\n");

  /*
   * Connectivity
   */
  fprintf(file, "        <DataArray type=\"%s\" Name=\"connectivity\""
                " format=\"%s\">\n",
          BFAM_LOCIDX_VTK, format);
  fprintf(file, "        </DataArray>\n");

  /*
   * Offsets
   */
  fprintf(file, "        <DataArray type=\"%s\" Name=\"offsets\""
                " format=\"%s\">\n",
          BFAM_LOCIDX_VTK, format);
  fprintf(file, "        </DataArray>\n");

  /*
   * Types
   */
  fprintf(file, "        <DataArray type=\"UInt8\" Name=\"types\""
                " format=\"%s\">\n",
          format);
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "      </Cells>\n");
  fprintf(file, "    </Piece>\n");
}

void bfam_vtk_write_file(bfam_domain_t *domain, bfam_domain_match_t match,
                         const char **tags, const char *directory,
                         const char *prefix, bfam_real_t time,
                         const char **scalars, const char **vectors,
                         const char **components, int binary, int compress,
                         int Np_write)
{
  const int endian = bfam_endian();

  bfam_locidx_t numElements = domain->num_subdomains;

  int rank, size;
  BFAM_MPI_CHECK(MPI_Comm_rank(domain->comm, &rank));
  BFAM_MPI_CHECK(MPI_Comm_size(domain->comm, &size));

  if (rank == 0)
    bfam_vtk_write_file_pvtu(size, directory, prefix, scalars, vectors,
                             components, binary, compress);

  bfam_subdomain_t **subdomains =
      bfam_malloc(numElements * sizeof(bfam_subdomain_t *));
  bfam_locidx_t num_subdomains;

  bfam_domain_get_subdomains(domain, match, tags, numElements, subdomains,
                             &num_subdomains);

  char filename[BFAM_BUFSIZ];
  if (directory)
    snprintf(filename, BFAM_BUFSIZ, "%s/" BFAM_VTK_VTU_FORMAT, directory,
             prefix, rank);
  else
    snprintf(filename, BFAM_BUFSIZ, BFAM_VTK_VTU_FORMAT, prefix, rank);

  BFAM_VERBOSE("Writing file: '%s'", filename);
  FILE *file = fopen(filename, "w");

  if (file == NULL)
  {
    BFAM_LERROR("Could not open %s for output!\n", filename);
    return;
  }

  fprintf(file, "<?xml version=\"1.0\"?>\n");
  fprintf(file, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"");

  if (binary && compress)
    fprintf(file, " compressor=\"vtkZLibDataCompressor\"");

  if (endian == BFAM_BIG_ENDIAN)
    fprintf(file, " byte_order=\"BigEndian\">\n");
  else
    fprintf(file, " byte_order=\"LittleEndian\">\n");

  fprintf(file, "  <UnstructuredGrid>\n");

  int files_written = 0;
  for (bfam_locidx_t s = 0; s < num_subdomains; ++s)
  {
    bfam_subdomain_t *subdomain = subdomains[s];

    if (subdomain->vtk_write_vtu_piece)
    {
      files_written += subdomain->vtk_write_vtu_piece(
          subdomains[s], file, time, scalars, vectors, components, binary,
          compress, rank, s, Np_write);
    }
    else
    {
      BFAM_WARNING("Subdomain: %s does not implement vtk_write_vtu_piece",
                   subdomain->name);
    }
  }

  if (files_written == 0)
    bfam_vtk_write_vtu_empty(file, binary);

  fprintf(file, "  </UnstructuredGrid>\n");
  fprintf(file, "</VTKFile>\n");

  if (ferror(file))
  {
    BFAM_LERROR("Error writing to %s\n", filename);
  }

  if (fclose(file))
  {
    BFAM_LERROR("Error closing %s\n", filename);
  }

  bfam_free(subdomains);
}

// }}}

// {{{ pcg32
// *Really* minimal PCG32 code / (c) 2014 M.E. O'Neill / pcg-random.org
// Licensed under Apache License 2.0 (NO WARRANTY, etc. see website)
// Copied from: http://www.pcg-random.org/
// Copied from: https://github.com/imneme/pcg-c-basic

static uint32_t bfam_pcg32_random_r(bfam_pcg32_random_t *rng)
{
  uint64_t oldstate = rng->state;
  // Advance internal state
  rng->state = oldstate * 6364136223846793005ULL + (rng->inc | 1);
  // Calculate output function (XSH RR), uses old state for max ILP
  uint32_t xorshifted = (uint32_t)(((oldstate >> 18u) ^ oldstate) >> 27u);
  uint32_t rot = (uint32_t)(oldstate >> 59u);
  return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

void bfam_pcg32_srandom_r(bfam_pcg32_random_t *rng, uint64_t initstate,
                          uint64_t initseq)
{
  rng->state = 0U;
  rng->inc = (initseq << 1u) | 1u;
  bfam_pcg32_random_r(rng);
  rng->state += initstate;
  bfam_pcg32_random_r(rng);
}

uint32_t bfam_pcg32_boundedrand_r(bfam_pcg32_random_t *rng, uint32_t bound)
{
  // To avoid bias, we need to make the range of the RNG a multiple of
  // bound, which we do by dropping output less than a threshold.
  // A naive scheme to calculate the threshold would be to do
  //
  //     uint32_t threshold = 0x100000000ull % bound;
  //
  // but 64-bit div/mod is slower than 32-bit div/mod (especially on
  // 32-bit platforms).  In essence, we do
  //
  //     uint32_t threshold = (0x100000000ull-bound) % bound;
  //
  // because this version will calculate the same modulus, but the LHS
  // value is less than 2^32.

  uint32_t threshold = -bound % bound;

  // Uniformity guarantees that this loop will terminate.  In practice, it
  // should usually terminate quickly; on average (assuming all bounds are
  // equally likely), 82.25% of the time, we can expect it to require just
  // one iteration.  In the worst case, someone passes a bound of 2^31 + 1
  // (i.e., 2147483649), which invalidates almost 50% of the range.  In
  // practice, bounds are typically small and only a tiny amount of the range
  // is eliminated.
  for (;;)
  {
    uint32_t r = bfam_pcg32_random_r(rng);
    if (r >= threshold)
      return r % bound;
  }
}
// }}}

#ifdef BFAM_USE_BFAMO

#ifdef STATS
stats_t _t;
#endif

const char *occa_modes[5] = {"SERIAL", "OpenMP", "OpenCL", "CUDA", NULL};

/* These names have to do with the type of h refinement */
const char *prj_base_names[6] = {"Pr",     /* Same in h */
                                 "Pr_cfb", /* coarsen from bottom */
                                 "Pr_cft", /* coarsen from top */
                                 "Pr_r2b", /* refine to bottom */
                                 "Pr_r2t", /* refine to top */
                                 NULL};

/* These names have to do with p refinement */
const char *bfamo_adapt_postfix[4] = {"R", /* refine */
                                      "S", /* same */
                                      "C", /*coarsen */
                                      NULL};

bfamo_sub_kernel_t bfamo_sub_kernel_null()
{
  bfamo_sub_kernel_t bsk;
  bsk.bk = NULL;
  bsk.dim = -1;
  bsk.items.x = -1;
  bsk.items.y = -1;
  bsk.items.z = -1;
  bsk.groups.x = -1;
  bsk.groups.y = -1;
  bsk.groups.z = -1;
  bsk.work1 = -1;
  bsk.work2 = -1;
  return bsk;
}

static inline void bfamo_kernel_set_stats(bfamo_kernel_t *bk, int flops1,
                                          int flops2, int memread)
{
#ifdef KERNEL_TIME
  BFAM_ASSERT(bk->flops1 == -1 || bk->flops1 == flops1);
  BFAM_ASSERT(bk->flops2 == -1 || bk->flops2 == flops2);
  BFAM_ASSERT(bk->memread == -1 || bk->memread == memread);

  bk->flops1 = flops1;
  bk->flops2 = flops2;
  bk->memread = memread;

  if (!sc_keyvalue_exists(_t.stat->kv, bk->stat_key))
    sc_statistics_add_empty(_t.stat, bk->stat_key);

  char memread_key[BFAM_BUFSIZ];
  snprintf(memread_key, BFAM_BUFSIZ, "Bandwidth (%d) :: %s", memread,
           bk->stat_key);
  size_t memread_key_len = strlen(memread_key) + 1;
  bk->memread_key = bfam_malloc(memread_key_len * sizeof(char));
  strncpy(bk->memread_key, memread_key, memread_key_len);
  if (!sc_keyvalue_exists(_t.stat->kv, bk->memread_key))
    sc_statistics_add_empty(_t.stat, bk->memread_key);

  char flops_key[BFAM_BUFSIZ];
  snprintf(flops_key, BFAM_BUFSIZ, "FLOPRate (%d, %d):: %s", flops1, flops2,
           bk->stat_key);
  size_t flops_key_len = strlen(flops_key) + 1;
  bk->flops_key = bfam_malloc(flops_key_len * sizeof(char));
  strncpy(bk->flops_key, flops_key, flops_key_len);
  if (!sc_keyvalue_exists(_t.stat->kv, bk->flops_key))
    sc_statistics_add_empty(_t.stat, bk->flops_key);
#endif
}

void bfamo_sub_kernel_set_stats(bfamo_sub_kernel_t *bsk, int work1, int work2,
                                int flops1, int flops2, int memread)
{
#ifdef KERNEL_TIME
  BFAM_ASSERT(bsk->work1 == -1);
  BFAM_ASSERT(bsk->work2 == -1);

  bsk->work1 = work1;
  bsk->work2 = work2;

  bfamo_kernel_set_stats(bsk->bk, flops1, flops2, memread);
#endif
}

bfamo_sub_kernel_t
bfamo_build_sub_kernel(bfam_dictionary_t *c_kernel_cache, const char *postfix,
                       bfamo_device_t *bdevice, const char *str,
                       const char *function_name, occaKernelInfo info)
{
  occaDevice device = bdevice->device;
  BFAM_ABORT_IF_NOT(postfix, "postfix not set for '%s' '%s'", str,
                    function_name);
  char key[BFAM_BUFSIZ];
  BFAM_ASSERT(sizeof(device) == sizeof(void *));
  snprintf(key, BFAM_BUFSIZ, "%p :: %s :: %s :: %s", (void *)device, str,
           function_name, postfix);

  bfamo_kernel_t *bk = bfam_dictionary_get_value_ptr(c_kernel_cache, key);

  if (!bk)
  {
    bk = bfam_malloc(sizeof(bfamo_kernel_t));
    bk->kernel = NULL;
    bk->dim = -1;
    bk->items.x = -1;
    bk->items.y = -1;
    bk->items.z = -1;
    bk->groups.x = -1;
    bk->groups.y = -1;
    bk->groups.z = -1;

    size_t key_len = strlen(key) + 1;
    bk->key = bfam_malloc(key_len * sizeof(char));
    strncpy(bk->key, key, key_len);
    BFAM_VERBOSE("Building kernel key '%s'", key);

    bk->flops_key = NULL;
    bk->memread_key = NULL;
    bk->stat_key = NULL;

    bk->kernel = occaDeviceBuildKernel(device, str, function_name, info);
    BFAMO_DICTIONARY_INSERT_PTR(c_kernel_cache, key, bk);

    bk->flops1 = -1;
    bk->flops2 = -1;
    bk->memread = -1;

#ifdef KERNEL_TIME
    /* We strip pointer since different MPI ranks will have different values */
    char stat_key[BFAM_BUFSIZ];
    snprintf(stat_key, BFAM_BUFSIZ, "%s :: %s :: %s", str, function_name,
             postfix);

    size_t stat_key_len = strlen(stat_key) + 1;
    bk->stat_key = bfam_malloc(stat_key_len * sizeof(char));
    strncpy(bk->stat_key, stat_key, stat_key_len);
    sc_statistics_add_empty(_t.stat, bk->stat_key);
#endif
  }

  bfamo_sub_kernel_t bsk = bfamo_sub_kernel_null();
  bsk.bk = bk;

  return bsk;
}

void bfamo_kernel_set_working_dims(bfamo_kernel_t *bk, int dim, occaDim items,
                                   occaDim groups)
{
  if (bk->dim != dim || bk->items.x != items.x || bk->items.y != items.y ||
      bk->items.z != items.z || bk->groups.x != groups.x ||
      bk->groups.y != groups.y || bk->groups.z != groups.z)
  {
    occaKernelSetWorkingDims(bk->kernel, dim, items, groups);
    bk->dim = dim;
    bk->items = items;
    bk->groups = groups;
  }
}

void bfamo_sub_kernel_set_working_dims(bfamo_sub_kernel_t *bsk, int dim,
                                       occaDim items, occaDim groups)
{
  bsk->dim = dim;
  bsk->items = items;
  bsk->groups = groups;
}

occaMemory bfamo_device_malloc(bfamo_device_t *bdevice, size_t bytecount,
                               void *src)
{
  occaDevice device = bdevice->device;
  double mem_fudge = bdevice->mem_fudge;
  bytecount = BFAM_MAX(1, bytecount);
  uintmax_t bytes = occaDeviceBytesAllocated(device);
  uintmax_t total_bytes = occaDeviceMemorySize(device);
  BFAM_ABORT_IF((double)(bytes + bytecount) > mem_fudge * (double)total_bytes,
                "Overmemory limit: "
                "current: allocated %ju (%.2f GiB) out of %ju (%.2f GiB) "
                "new val: allocated %ju (%.2f GiB) out of %ju (%.2f GiB) "
                "(fudge factor is: %.2f)",
                bytes, ((double)bytes) / GiB, total_bytes,
                ((double)total_bytes) / GiB, bytes + bytecount,
                ((double)(bytes + bytecount)) / GiB, total_bytes,
                ((double)total_bytes) / GiB, mem_fudge);
  if ((double)(bytes + bytecount) > 0.9 * mem_fudge * (double)total_bytes)
    BFAM_WARNING("At 90%% of memory limit: "
                 "current: allocated %ju (%.2f GiB) out of %ju (%.2f GiB) "
                 "new val: allocated %ju (%.2f GiB) out of %ju (%.2f GiB) "
                 "(fudge factor is: %.2f)",
                 bytes, ((double)bytes) / GiB, total_bytes,
                 ((double)total_bytes) / GiB, bytes + bytecount,
                 ((double)(bytes + bytecount)) / GiB, total_bytes,
                 ((double)total_bytes) / GiB, mem_fudge);

  return occaDeviceMalloc(device, bytecount, src);
}

occaKernelInfo bfamo_kernelinfo_new(int occa_mode)
{
  occaKernelInfo info = occaCreateKernelInfo();

  const bfamo_real_t p_PI = BFAMO_REAL(3.1415926535897932384626433832795);

  if (sizeof(bfamo_real_t) == sizeof(float))
  {
    occaKernelInfoAddDefine(info, "datafloat", occaString((char *)"float"));
    occaKernelInfoAddDefine(info, "datafloat4", occaString((char *)"float4"));
    occaKernelInfoAddDefine(info, "datafloatv", occaString((char *)"float4"));

    if (occa_mode == OPENCL)
      occaKernelInfoAddDefine(info, "DATAFLOAT_FMOD",
                              occaString((char *)"fmod"));
    else
      occaKernelInfoAddDefine(info, "DATAFLOAT_FMOD",
                              occaString((char *)"fmodf"));

    if (occa_mode == CUDA)
      occaKernelInfoAddDefine(info, "zero_datafloatv",
                              occaString((char *)"make_float4(0,0,0,0)"));
    else if (occa_mode == OPENMP)
      occaKernelInfoAddDefine(info, "zero_datafloatv",
                              occaString((char *)"float4(0,0,0,0)"));
    else
      occaKernelInfoAddDefine(info, "zero_datafloatv", occaString((char *)"0"));
  }
  else if (sizeof(bfamo_real_t) == sizeof(double))
  {
    occaKernelInfoAddDefine(info, "datafloat", occaString((char *)"double"));
    occaKernelInfoAddDefine(info, "datafloat4", occaString((char *)"double4"));
    occaKernelInfoAddDefine(info, "datafloatv", occaString((char *)"double4"));
    occaKernelInfoAddDefine(info, "DATAFLOAT_FMOD", occaString((char *)"fmod"));
    if (occa_mode == CUDA)
      occaKernelInfoAddDefine(info, "zero_datafloatv",
                              occaString((char *)"make_double4(0,0,0,0)"));
    else if (occa_mode == OPENMP)
      occaKernelInfoAddDefine(info, "zero_datafloatv",
                              occaString((char *)"double4(0,0,0,0)"));
    else
      occaKernelInfoAddDefine(info, "zero_datafloatv", occaString((char *)"0"));
  }
  else
  {
    BFAM_ABORT("Unknown floating point size");
  }
  if (sizeof(metric_real_t) == sizeof(double))
  {
    occaKernelInfoAddDefine(info, "metricfloat", occaString((char *)"double"));
  }
  else
  {
    BFAM_ABORT("Unknown metric floating point size");
  }

  occaKernelInfoAddDefine(info, "max_datafloat", occaReal(BFAMO_REAL_MAX));

  occaKernelInfoAddDefine(info, "DIM", occaInt(BFAMO_DIM));

  occaKernelInfoAddDefine(info, "ID_VGEO_JR0X0", occaInt(ID_VGEO_JR0X0));
  occaKernelInfoAddDefine(info, "ID_VGEO_JR0X1", occaInt(ID_VGEO_JR0X1));
#if BFAMO_DIM == 3
  occaKernelInfoAddDefine(info, "ID_VGEO_JR0X2", occaInt(ID_VGEO_JR0X2));
#endif
  occaKernelInfoAddDefine(info, "ID_VGEO_JR1X0", occaInt(ID_VGEO_JR1X0));
  occaKernelInfoAddDefine(info, "ID_VGEO_JR1X1", occaInt(ID_VGEO_JR1X1));
#if BFAMO_DIM == 3
  occaKernelInfoAddDefine(info, "ID_VGEO_JR1X2", occaInt(ID_VGEO_JR1X2));
  occaKernelInfoAddDefine(info, "ID_VGEO_JR2X0", occaInt(ID_VGEO_JR2X0));
  occaKernelInfoAddDefine(info, "ID_VGEO_JR2X1", occaInt(ID_VGEO_JR2X1));
  occaKernelInfoAddDefine(info, "ID_VGEO_JR2X2", occaInt(ID_VGEO_JR2X2));
#endif
  occaKernelInfoAddDefine(info, "ID_VGEO_JINV", occaInt(ID_VGEO_JINV));
  occaKernelInfoAddDefine(info, "ID_VGEO_W", occaInt(ID_VGEO_W));
  occaKernelInfoAddDefine(info, "ID_VGEO_X0", occaInt(ID_VGEO_X0));
  occaKernelInfoAddDefine(info, "ID_VGEO_X1", occaInt(ID_VGEO_X1));
#if BFAMO_DIM == 3
  occaKernelInfoAddDefine(info, "ID_VGEO_JINV2", occaInt(ID_VGEO_JINV2));
  occaKernelInfoAddDefine(info, "ID_VGEO_X2", occaInt(ID_VGEO_X2));
#endif
  occaKernelInfoAddDefine(info, "NVGEO", occaInt(NVGEO));

  occaKernelInfoAddDefine(info, "ID_SGEO_NX0", occaInt(ID_SGEO_NX0));
  occaKernelInfoAddDefine(info, "ID_SGEO_NX1", occaInt(ID_SGEO_NX1));
#if BFAMO_DIM == 3
  occaKernelInfoAddDefine(info, "ID_SGEO_NX2", occaInt(ID_SGEO_NX2));
#endif
  occaKernelInfoAddDefine(info, "ID_SGEO_SJWJ", occaInt(ID_SGEO_SJWJ));
  occaKernelInfoAddDefine(info, "NSGEO", occaInt(NSGEO));

  occaKernelInfoAddDefine(info, "ID_GGEO_NX0", occaInt(ID_GGEO_NX0));
  occaKernelInfoAddDefine(info, "ID_GGEO_NX1", occaInt(ID_GGEO_NX1));
#if BFAMO_DIM == 3
  occaKernelInfoAddDefine(info, "ID_GGEO_NX2", occaInt(ID_GGEO_NX2));
#endif

  occaKernelInfoAddDefine(info, "ID_GGEO_SJW", occaInt(ID_GGEO_SJW));

  occaKernelInfoAddDefine(info, "ID_GGEO_X0", occaInt(ID_GGEO_X0));
  occaKernelInfoAddDefine(info, "ID_GGEO_X1", occaInt(ID_GGEO_X1));
#if BFAMO_DIM == 3
  occaKernelInfoAddDefine(info, "ID_GGEO_X2", occaInt(ID_GGEO_X2));
#endif
  occaKernelInfoAddDefine(info, "NGGEO", occaInt(NGGEO));

  occaKernelInfoAddDefine(info, "P4EST_CHILDREN", occaInt(P4EST_CHILDREN));

  occaKernelInfoAddDefine(info, "ID_QUAD_TID", occaInt(ID_QUAD_TID));
  occaKernelInfoAddDefine(info, "ID_QUAD_LVL", occaInt(ID_QUAD_LVL));
  occaKernelInfoAddDefine(info, "ID_QUAD_X", occaInt(ID_QUAD_X));
  occaKernelInfoAddDefine(info, "ID_QUAD_Y", occaInt(ID_QUAD_Y));
  occaKernelInfoAddDefine(info, "ID_QUAD_Z", occaInt(ID_QUAD_Z));
  occaKernelInfoAddDefine(info, "NQUAD", occaInt(NQUAD));
  occaKernelInfoAddDefine(info, "QUAD_ROOT_LEN", occaInt(P4EST_ROOT_LEN));

  occaKernelInfoAddDefine(info, "p_PI", occaReal(p_PI));

  return info;
}

void bfamo_add_operator_metric_real(bfamo_device_t *bdevice,
                                    bfam_dictionary_t *c_ops,
                                    bfam_dictionary_t *d_ops,
                                    const char *prefix, const int N,
                                    const int Np)
{
  char name[BFAM_BUFSIZ];

  snprintf(name, BFAM_BUFSIZ, "%s_%d", prefix, N);
  if (!bfam_dictionary_contains(c_ops, name))
  {
    bfam_long_real_t *b_op = bfam_dictionary_get_value_ptr(d_ops, name);
    BFAM_ABORT_IF_NOT(b_op, name);
    metric_real_t *h_op = bfam_malloc_aligned(Np * sizeof(metric_real_t));

    for (int n = 0; n < Np; n++)
      h_op[n] = (metric_real_t)b_op[n];

    occaMemory c_op =
        bfamo_device_malloc(bdevice, Np * sizeof(metric_real_t), h_op);
    BFAMO_DICTIONARY_INSERT_PTR(c_ops, name, c_op);
    BFAMO_BFAM_FREE_ALIGNED(h_op);
  }
}

void bfamo_add_operator_bfamo_real(bfamo_device_t *bdevice,
                                   bfam_dictionary_t *c_ops,
                                   bfam_dictionary_t *d_ops, const char *prefix,
                                   const int N, const int Np)
{
  char name[BFAM_BUFSIZ];

  snprintf(name, BFAM_BUFSIZ, "%s_%d", prefix, N);
  if (!bfam_dictionary_contains(c_ops, name))
  {
    bfam_real_t *b_op = bfam_dictionary_get_value_ptr(d_ops, name);
    BFAM_ABORT_IF_NOT(b_op, name);
    bfamo_real_t *h_op = bfam_malloc_aligned(Np * sizeof(bfamo_real_t));

    for (int n = 0; n < Np; n++)
      h_op[n] = (bfamo_real_t)b_op[n];

    occaMemory c_op =
        bfamo_device_malloc(bdevice, Np * sizeof(bfamo_real_t), h_op);
    BFAMO_DICTIONARY_INSERT_PTR(c_ops, name, c_op);
    BFAMO_BFAM_FREE_ALIGNED(h_op);
  }
}

static void bfamo_add_projection_bfamo_real(bfamo_device_t *bdevice,
                                            bfam_dictionary_t *c_ops,
                                            bfam_dictionary_t *N2N,
                                            const int N_src, const int N_dst)
{
  char name[BFAM_BUFSIZ];

  snprintf(name, BFAM_BUFSIZ, "Pr_%d_%d", N_src, N_dst);
  if (!bfam_dictionary_contains(c_ops, name))
  {
    bfam_subdomain_dgx_interpolator_t *b_interpolator =
        bfam_subdomain_dgx_get_interpolator(N2N, N_src, N_dst, 0);

    BFAM_ASSERT(b_interpolator);
    BFAM_ASSERT(b_interpolator->N_src == N_src);
    BFAM_ASSERT(b_interpolator->N_dst == N_dst);
    bfam_real_t **b_prj = b_interpolator->prj;
#if BFAM_DEBUG
    bfam_real_t **b_wi_mass_prj = b_interpolator->wi_mass_prj;
    BFAM_ASSERT(b_prj);
    BFAM_ASSERT(b_wi_mass_prj);
#endif

    const int Nq_src = N_src + 1;
    const int Nq_dst = N_dst + 1;
    const int Np = Nq_src * Nq_dst;
    bfamo_real_t *h_op = bfam_malloc_aligned(Np * sizeof(bfamo_real_t));

    const int num_prj = b_interpolator->num_prj;
    BFAM_ASSERT(prj_base_names[num_prj] == NULL);
    for (int k = 0; k < num_prj; k++)
    {
      BFAM_ASSERT(prj_base_names[k]);
      snprintf(name, BFAM_BUFSIZ, "%s_%d_%d", prj_base_names[k], N_src, N_dst);
      bfam_real_t *b_op = b_prj[k];
      if (!b_op)
      {
        BFAM_ASSERT(N_src == N_dst);
        for (int i = 0; i < Nq_src; i++)
          for (int j = 0; j < Nq_src; j++)
          {
            const int n = i + j * Nq_src;
            if (i == j)
              h_op[n] = 1;
            else
              h_op[n] = 0;
          }
      }
      else
        for (int n = 0; n < Np; n++)
          h_op[n] = (bfamo_real_t)b_op[n];

      occaMemory c_op =
          bfamo_device_malloc(bdevice, Np * sizeof(bfamo_real_t), h_op);
      BFAMO_DICTIONARY_INSERT_PTR(c_ops, name, c_op);
    }

    BFAMO_BFAM_FREE_ALIGNED(h_op);
  }
}

void bfamo_add_proj_operators(bfamo_device_t *bdevice, bfam_dictionary_t *c_ops,
                              bfam_dictionary_t *N2N, const int N_src,
                              const int N_dst)
{
  bfamo_add_projection_bfamo_real(bdevice, c_ops, N2N, N_src, N_dst);
  bfamo_add_projection_bfamo_real(bdevice, c_ops, N2N, N_dst, N_src);
}

#ifdef COMPUTE_STATS
void bfamo_add_statisitics(bfam_domain_t *domain,
                           bfam_communicator_t *inter_comm, occaDevice device,
                           const char **tags, stats_t *tm)
{
  sc_statistics_add_empty(tm->stat, "compute time");
  sc_statistics_add_empty(tm->stat, "adapt time");
  sc_statistics_add_empty(tm->stat, "adapt mark time");
  sc_statistics_add_empty(tm->stat, "adapt coarsen time");
  sc_statistics_add_empty(tm->stat, "adapt refine time");
  sc_statistics_add_empty(tm->stat, "adapt grow time");
  sc_statistics_add_empty(tm->stat, "adapt repartition coarsen time");
  sc_statistics_add_empty(tm->stat, "adapt repartition refine time");

  bfam_subdomain_t **subdomains =
      bfam_malloc(domain->num_subdomains * sizeof(bfam_subdomain_t **));
  bfam_locidx_t num_subdomains = 0;
  bfam_domain_get_subdomains(domain, BFAM_DOMAIN_OR, tags,
                             domain->num_subdomains, subdomains,
                             &num_subdomains);

  uintmax_t num_elm = 0;
  uintmax_t num_pnt = 0;
  for (bfam_locidx_t s = 0; s < num_subdomains; ++s)
  {
    bfam_subdomain_dgx_t *sub = (bfam_subdomain_dgx_t *)subdomains[s];
    num_elm += sub->K;
    num_pnt += sub->K * sub->Np;
  }

  BFAMO_BFAM_FREE(subdomains);

  sc_statistics_add_empty(tm->stat, "number of elements");
  sc_statistics_accumulate(tm->stat, "number of elements", (double)num_elm);

  sc_statistics_add_empty(tm->stat, "number of grid points");
  sc_statistics_accumulate(tm->stat, "number of grid points", (double)num_pnt);

  if (inter_comm->num_procs > 0)
  {
    sc_statistics_add_empty(tm->stat, "number of neighgors");
    sc_statistics_accumulate(tm->stat, "number of neighgors",
                             (double)inter_comm->num_procs);

    sc_statistics_add_empty(tm->stat, "total send size");
    sc_statistics_accumulate(tm->stat, "total send size",
                             (double)inter_comm->send_sz);

    sc_statistics_add_empty(tm->stat, "total recv size");
    sc_statistics_accumulate(tm->stat, "total recv size",
                             (double)inter_comm->recv_sz);

    sc_statistics_add_empty(tm->stat, "average send size");
    sc_statistics_accumulate(tm->stat, "average send size",
                             ((double)inter_comm->send_sz) /
                                 inter_comm->num_procs);

    sc_statistics_add_empty(tm->stat, "average recv size");
    sc_statistics_accumulate(tm->stat, "average recv size",
                             ((double)inter_comm->recv_sz) /
                                 inter_comm->num_procs);
  }

  sc_statistics_add_empty(tm->stat, "device memory allocated");
  uintmax_t bytes = occaDeviceBytesAllocated(device);
  sc_statistics_accumulate(tm->stat, "device memory allocated",
                           ((double)bytes) / GiB);
}
#endif

/* TODO: Move all these calls to the communicator */
void bfamo_communicator_post_send(bfam_communicator_t *comm)
{
  for (int p = 0; p < comm->num_procs; p++)
    BFAM_MPI_CHECK(MPI_Isend(comm->proc_data[p].send_buf,
                             (int)comm->proc_data[p].send_sz, MPI_BYTE,
                             comm->proc_data[p].rank, comm->tag, comm->comm,
                             &comm->send_request[p]));
}

void bfamo_communicator_post_recv(bfam_communicator_t *comm)
{
  for (int p = 0; p < comm->num_procs; p++)
    BFAM_MPI_CHECK(MPI_Irecv(comm->proc_data[p].recv_buf,
                             (int)comm->proc_data[p].recv_sz, MPI_BYTE,
                             comm->proc_data[p].rank, comm->tag, comm->comm,
                             &comm->recv_request[p]));
}

void bfamo_communicator_send_wait(bfam_communicator_t *comm)
{
  BFAM_MPI_CHECK(
      MPI_Waitall(comm->num_procs, comm->send_request, comm->send_status));
}

void bfamo_communicator_recv_wait(bfam_communicator_t *comm)
{
  BFAM_MPI_CHECK(
      MPI_Waitall(comm->num_procs, comm->recv_request, comm->recv_status));
}
/* end TODO */

int bfamo_free_c_fields(const char *key, void *val, void *arg)
{
  const char **keep = (const char **)arg;
  if (keep)
    for (int n = 0; keep[n]; n++)
      if (0 == strcmp(key, keep[n]))
        return 0;
  BFAM_LDEBUG("Free occaMemory : %s (%p)", key, val);
  occaMemoryFree((occaMemory)val);
  return 0;
}

int bfamo_free_c_kernel_cache(const char *key, void *val, void *arg)
{
  bfamo_kernel_t *bk = (bfamo_kernel_t *)val;
  BFAM_LDEBUG("Free bfamo_kernel_free : %s %s -> (%p)", key, bk->key, val);
  BFAMO_BFAM_FREE(bk->key);
  BFAMO_BFAM_FREE(bk->stat_key);
  BFAMO_BFAM_FREE(bk->memread_key);
  BFAMO_BFAM_FREE(bk->flops_key);
  BFAMO_OCCA_KERNEL_FREE(bk->kernel);
  BFAMO_BFAM_FREE(bk);
  return 0;
}

void bfamo_prefs_print_file(const char *filename, int rank, int loglevel)
{
  size_t len;
  char *buffer = bfam_util_read_file(filename, &len);
  BFAM_ROOT_INFO("");
  BFAM_ROOT_INFO("----- Preferences File -----------------------------------");
  if (loglevel <= BFAM_LL_INFO && rank == 0)
    fwrite(buffer, sizeof(char), len, stdout);
  BFAM_ROOT_INFO("----------------------------------------------------------");
  BFAMO_BFAM_FREE(buffer);
}

int bfamo_weight_fn(p4est_t *p4est, p4est_topidx_t which_tree,
                    p4est_quadrant_t *quadrant)
{
  bfam_pxest_user_data_t *ud = quadrant->p.user_data;
  const int N = ud->N;

#if BFAMO_DIM == 2
  const int Np = (N + 1) * (N + 1);
#elif BFAMO_DIM == 3
  const int Np = (N + 1) * (N + 1) * (N + 1);
#else
#error "Bad Dimension"
#endif

  /* TODO may need to adjust weight for multiple physics */

  return Np;
}

void bfamo_lsrk_free(bfamo_lsrk_t *ts)
{
  ts->n_stages = 0;
  BFAMO_BFAM_FREE_ALIGNED(ts->A);
  BFAMO_BFAM_FREE_ALIGNED(ts->B);
  BFAMO_BFAM_FREE_ALIGNED(ts->C);
}

void bfamo_lsrk_init(bfamo_lsrk_t *ts, bfam_domain_t *dom,
                     const char *lsrk_prefix, bfamo_lsrk_method_t method,
                     bfam_domain_match_t subdom_match, const char **subdom_tags,
                     bfamo_aux_rates_t aux_rates, void *user_data)
{

  /*
   * get the subdomains and create rates we will need
   */
  bfam_subdomain_t *subs[dom->num_subdomains + 1];
  bfam_locidx_t num_subs = 0;
  bfam_domain_get_subdomains(dom, subdom_match, subdom_tags,
                             dom->num_subdomains, subs, &num_subs);
  for (int s = 0; s < num_subs; s++)
    aux_rates(subs[s], lsrk_prefix, user_data);

  switch (method)
  {
  default:
    BFAM_WARNING("Invalid LSRK scheme, using KC54");
  case BFAMO_TS_LSRK_KC54:
    ts->n_stages = 5;
    ts->A = bfam_malloc_aligned(ts->n_stages * sizeof(bfamo_real_t));
    ts->B = bfam_malloc_aligned(ts->n_stages * sizeof(bfamo_real_t));
    ts->C = bfam_malloc_aligned((ts->n_stages + 1) * sizeof(bfamo_real_t));

    ts->A[0] = (bfamo_real_t)(0.0L);
    ts->A[1] = (bfamo_real_t)(-567301805773.0L / 1357537059087.0L);
    ts->A[2] = (bfamo_real_t)(-2404267990393.0L / 2016746695238.0L);
    ts->A[3] = (bfamo_real_t)(-3550918686646.0L / 2091501179385.0L);
    ts->A[4] = (bfamo_real_t)(-1275806237668.0L / 842570457699.0L);

    ts->B[0] = (bfamo_real_t)(1432997174477.0L / 9575080441755.0L);
    ts->B[1] = (bfamo_real_t)(5161836677717.0L / 13612068292357.0L);
    ts->B[2] = (bfamo_real_t)(1720146321549.0L / 2090206949498.0L);
    ts->B[3] = (bfamo_real_t)(3134564353537.0L / 4481467310338.0L);
    ts->B[4] = (bfamo_real_t)(2277821191437.0L / 14882151754819.0L);

    ts->C[0] = (bfamo_real_t)(0.0L);
    ts->C[1] = (bfamo_real_t)(1432997174477.0L / 9575080441755.0L);
    ts->C[2] = (bfamo_real_t)(2526269341429.0L / 6820363962896.0L);
    ts->C[3] = (bfamo_real_t)(2006345519317.0L / 3224310063776.0L);
    ts->C[4] = (bfamo_real_t)(2802321613138.0L / 2924317926251.0L);
    ts->C[5] = (bfamo_real_t)(1.0L);
    break;
  case BFAMO_TS_LSRK_W33:
    ts->n_stages = 3;
    ts->A = bfam_malloc_aligned(ts->n_stages * sizeof(bfamo_real_t));
    ts->B = bfam_malloc_aligned(ts->n_stages * sizeof(bfamo_real_t));
    ts->C = bfam_malloc_aligned((ts->n_stages + 1) * sizeof(bfamo_real_t));

    ts->A[0] = (bfamo_real_t)(0.0L);
    ts->A[1] = (bfamo_real_t)(-5.0L / 9.0L);
    ts->A[2] = (bfamo_real_t)(-153.0L / 128.0L);

    ts->B[0] = (bfamo_real_t)(1.0L / 3.0L);
    ts->B[1] = (bfamo_real_t)(15.0L / 16.0L);
    ts->B[2] = (bfamo_real_t)(8.0L / 15.0L);

    ts->C[0] = (bfamo_real_t)(0.0L);
    ts->C[1] = (bfamo_real_t)(1.0L / 3.0L);
    ts->C[2] = (bfamo_real_t)(3.0L / 4.0L);
    ts->C[3] = (bfamo_real_t)(1.0L);
    break;
  case BFAMO_TS_LSRK_HEUN:
    ts->n_stages = 2;
    ts->A = bfam_malloc_aligned(ts->n_stages * sizeof(bfamo_real_t));
    ts->B = bfam_malloc_aligned(ts->n_stages * sizeof(bfamo_real_t));
    ts->C = bfam_malloc_aligned((ts->n_stages + 1) * sizeof(bfamo_real_t));

    ts->A[0] = (bfamo_real_t)(0.0L);
    ts->A[1] = -(bfamo_real_t)(1.0L);

    ts->B[0] = (bfamo_real_t)(1.0L);
    ts->B[1] = (bfamo_real_t)(1.0L / 2.0L);

    ts->C[0] = (bfamo_real_t)(0.0L);
    ts->C[1] = (bfamo_real_t)(1.0L);
    ts->C[2] = (bfamo_real_t)(1.0L);
    break;
  case BFAMO_TS_LSRK_FE:
    ts->n_stages = 1;
    ts->A = bfam_malloc_aligned(sizeof(bfamo_real_t));
    ts->B = bfam_malloc_aligned(sizeof(bfamo_real_t));
    ts->C = bfam_malloc_aligned((ts->n_stages + 1) * sizeof(bfamo_real_t));

    ts->A[0] = (bfamo_real_t)(0.0L);

    ts->B[0] = (bfamo_real_t)(1.0L);

    ts->C[0] = (bfamo_real_t)(0.0L);
    ts->C[1] = (bfamo_real_t)(1.0L);
    break;
  }
}
#endif
