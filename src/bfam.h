#ifndef BFAM_H
#define BFAM_H

#include <inttypes.h>
#include <mpi.h>
#include <stdint.h>
#include <stdio.h>

// {{{ Types:

// }}}

// {{{ base

#if defined __GNUC__ && !defined __GNUC_PREREQ
#ifndef __GNUC_MINOR__
#define __GNUC_PREREQ(maj, min) 0
#else
#define __GNUC_PREREQ(maj, min)                                                \
  ((__GNUC__ << 16) + __GNUC_MINOR__ >= ((maj) << 16) + (min))
#endif
#endif

#define BFAM_BUFSIZ 8192

#define BFAM_APPEND(x, y) x##y
#define BFAM_APPEND_EXPAND(x, y) BFAM_APPEND(x, y)

#define BFAM_MIN(a, b) (((a) < (b)) ? (a) : (b))
#define BFAM_MAX(a, b) (((a) > (b)) ? (a) : (b))

#define BFAM_NOOP()                                                            \
  do                                                                           \
  {                                                                            \
  } while (0)
#define BFAM_ABORT(...) bfam_abort_verbose(__FILE__, __LINE__, __VA_ARGS__)
#define BFAM_ABORT_IF(q, ...) ((q) ? BFAM_ABORT(__VA_ARGS__) : (void)0)
#define BFAM_ABORT_IF_NOT(q, ...) BFAM_ABORT_IF(!(q), __VA_ARGS__)

#ifdef BFAM_DEBUG
#define BFAM_ASSERT(expression)                                                \
  BFAM_ABORT_IF_NOT((expression), "Assert Failed: '" #expression "'")
#else
#define BFAM_ASSERT(expression) BFAM_NOOP()
#endif

#define BFAM_SYS_ERROR_CHECK(cond, msg)                                        \
  do                                                                           \
  {                                                                            \
    if (cond)                                                                  \
    {                                                                          \
      perror(msg);                                                             \
      BFAM_ABORT("perror");                                                    \
    }                                                                          \
  } while (0)

#define BFAM_MPI_CHECK(c) BFAM_ABORT_IF_NOT((c) == MPI_SUCCESS, "MPI Error")

#define BFAM_IS_ALIGNED(p, a) (((intptr_t)(p) & ((a)-1)) == 0)

/*
 * Switch for different compilers taken from this web page:
 *
 *     http://nadeausoftware.com/articles/2012/10/c_c_tip_how_detect_compiler_name_and_version_using_compiler_predefined_macros
 *
 */
#if defined(__clang__)
/* Clang/LLVM. ---------------------------------------------- */
#define BFAM_ASSUME_ALIGNED(lvalueptr, align) BFAM_NOOP()
#elif defined(__ICC) || defined(__INTEL_COMPILER)
/* Intel ICC/ICPC. ------------------------------------------ */
#define BFAM_ASSUME_ALIGNED(lvalueptr, align)                                  \
  __assume_aligned(lvalueptr, align);                                          \
  BFAM_ASSERT(BFAM_IS_ALIGNED(lvalueptr, align))

#elif defined(__GNUC__) || defined(__GNUG__)
/* GNU GCC/G++. --------------------------------------------- */
#if __GNUC_PREREQ(4, 7)
//      If  gcc_version >= 4.7
#define BFAM_ASSUME_ALIGNED(lvalueptr, align)                                  \
  lvalueptr = __builtin_assume_aligned(lvalueptr, align);                      \
  BFAM_ASSERT(BFAM_IS_ALIGNED(lvalueptr, align))
#else
//       Else
#define BFAM_ASSUME_ALIGNED(lvalueptr, align) BFAM_NOOP()
#endif

#elif defined(__HP_cc) || defined(__HP_aCC)
/* Hewlett-Packard C/aC++. ---------------------------------- */
#define BFAM_ASSUME_ALIGNED(lvalueptr, align) BFAM_NOOP()

#elif defined(__IBMC__) || defined(__IBMCPP__)
/* IBM XL C/C++. -------------------------------------------- */
#define BFAM_ASSUME_ALIGNED(lvalueptr, align) BFAM_NOOP()

#elif defined(_MSC_VER)
/* Microsoft Visual Studio. --------------------------------- */
#define BFAM_ASSUME_ALIGNED(lvalueptr, align) BFAM_NOOP()

#elif defined(__PGI)
/* Portland Group PGCC/PGCPP. ------------------------------- */
#define BFAM_ASSUME_ALIGNED(lvalueptr, align) BFAM_NOOP()

#elif defined(__SUNPRO_C) || defined(__SUNPRO_CC)
/* Oracle Solaris Studio. ----------------------------------- */
#define BFAM_ASSUME_ALIGNED(lvalueptr, align) BFAM_NOOP()

#endif

#if (defined __GNUC__) || (defined __PGI) || (defined __IBMC__)
#define BFAM_ALIGN(n) __attribute__((aligned(n)))
#elif (defined _MSC_VER)
#define BFAM_ALIGN(n) __declspec(align(n))
#else
#error Need equilvent of __attribute__((aligned(n))) for this compiler
#endif

#if defined __GNUC__
#define BFAM_ASM_COMMENT(X) __asm__("# " X)
#else
#define BFAM_ASM_COMMENT(X)
#endif

#ifndef BFAM_NORETURN
#if defined(__clang__)
#if __has_feature(attribute_analyzer_noreturn)
#define BFAM_NORETURN __attribute__((analyzer_noreturn))
#else
#define BFAM_NORETURN
#endif
#else
#define BFAM_NORETURN
#endif
#endif

/* Type for setup computations */
#ifdef BFAM_USE_LONG_DOUBLE
typedef long double bfam_long_real_t;
#define BFAM_LONG_REAL(x) BFAM_APPEND(x, L)
#define BFAM_LONG_REAL_PRIe "Le"
#define BFAM_LONG_REAL_PRIf "Lf"
#define BFAM_LONG_REAL_PRIg "Lg"

#define BFAM_LONG_REAL_ABS fabsl
#define BFAM_LONG_REAL_ATAN atanl
#define BFAM_LONG_REAL_COS cosl
#define BFAM_LONG_REAL_LGAMMA lgammal
#define BFAM_LONG_REAL_LOG logl
#define BFAM_LONG_REAL_EXP expl
#define BFAM_LONG_REAL_HYPOT hypotl
#define BFAM_LONG_REAL_SQRT sqrtl
#define BFAM_LONG_REAL_SIN sinl
#define BFAM_LONG_REAL_TAN tanl

#define BFAM_LONG_REAL_EPS LDBL_EPSILON
#define BFAM_LONG_REAL_MIN LDBL_MIN
#define BFAM_LONG_REAL_MAX LDBL_MAX
#else
typedef double bfam_long_real_t;
#define BFAM_LONG_REAL(x) x
#define BFAM_LONG_REAL_PRIe "e"
#define BFAM_LONG_REAL_PRIf "f"
#define BFAM_LONG_REAL_PRIg "g"

#define BFAM_LONG_REAL_ABS fabs
#define BFAM_LONG_REAL_ATAN atan
#define BFAM_LONG_REAL_COS cos
#define BFAM_LONG_REAL_LGAMMA lgamma
#define BFAM_LONG_REAL_LOG log
#define BFAM_LONG_REAL_EXP exp
#define BFAM_LONG_REAL_HYPOT hypot
#define BFAM_LONG_REAL_SQRT sqrt
#define BFAM_LONG_REAL_SIN sin
#define BFAM_LONG_REAL_TAN tan

#define BFAM_LONG_REAL_EPS DBL_EPSILON
#define BFAM_LONG_REAL_MIN DBL_MIN
#define BFAM_LONG_REAL_MAX DBL_MAX
#endif
#define BFAM_LONG_REAL_PI (4 * BFAM_LONG_REAL_ATAN(1))

#define BFAM_LONG_REAL_APPROX_EQ(x, y, K)                                      \
  BFAM_APPROX_EQ((x), (y), (K), BFAM_LONG_REAL_ABS, BFAM_LONG_REAL_EPS,        \
                 BFAM_LONG_REAL_EPS *BFAM_LONG_REAL_EPS)

/* Type for runtime computations */
typedef double bfam_real_t;
#define bfam_real_nan nan
#define BFAM_REAL_MPI MPI_DOUBLE
#define BFAM_REAL(x) BFAM_APPEND(x, )
#define BFAM_REAL_PRIe "e"
#define BFAM_REAL_PRIf "f"
#define BFAM_REAL_PRIg "g"

#define BFAM_REAL_NAN(x) nan(x)

#define BFAM_REAL_ABS fabs
#define BFAM_REAL_ISFINITE isfinite
#define BFAM_REAL_SQRT sqrt
#define BFAM_REAL_EXP exp
#define BFAM_REAL_COS cos
#define BFAM_REAL_SIN sin
#define BFAM_REAL_ASINH asinh
#define BFAM_REAL_HYPOT hypot
#define BFAM_REAL_HYPOT3(x, y, z) hypot((x), hypot((y), (z)))

#define BFAM_REAL_EPS DBL_EPSILON
#define BFAM_REAL_MIN DBL_MIN
#define BFAM_REAL_MAX DBL_MAX

#define BFAM_REAL_VTK "Float64"
#define BFAM_REAL_FMTe "24.16e"

#define BFAM_REAL_APPROX_EQ(x, y, K)                                           \
  BFAM_APPROX_EQ((x), (y), (K), BFAM_REAL_ABS, BFAM_REAL_EPS, BFAM_REAL_EPS)

/* Type for processor-local indexing */
typedef int32_t bfam_locidx_t;
#define BFAM_LOCIDX_MPI MPI_INT
#define BFAM_LOCIDX_VTK "Int32"
#define BFAM_LOCIDX_MIN INT32_MIN
#define BFAM_LOCIDX_MAX INT32_MAX
#define BFAM_LOCIDX_PRId PRId32
#define BFAM_LOCIDX_SCNd SCNd32

/* Type for globally unique indexing */
typedef int64_t bfam_gloidx_t;
#define BFAM_GLOIDX_MPI MPI_LONG_LONG_INT
#define BFAM_GLOIDX_VTK "Int64"
#define BFAM_GLOIDX_MIN INT64_MIN
#define BFAM_GLOIDX_MAX INT64_MAX
#define BFAM_GLOIDX_PRId PRId64
#define BFAM_GLOIDX_SCNd SCNd64

/* Flags for adaptation */
#define BFAM_FLAG_COARSEN (1 << 0)
#define BFAM_FLAG_REFINE (1 << 1)
#define BFAM_FLAG_ADAPTED (1 << 2)
#define BFAM_FLAG_SAME (0)

#ifndef BFAM_DGX_DIMENSION
#define BFAM_DGX_DIMENSION 3
#endif

/** Verbose abort function.
 *
 * Typically this function is not called directly but the \c BFAM_ABORT,
 * \c BFAM_ABORT_IF, and \c BFAM_ABORT_IF_NOT macros are used.
 *
 * \param[in] file file name the abort occurred in (can be obtained from
 *                 \c __FILE__)
 * \param[in] line line number the abort occurred at (can be obtained from
 *                 \c __LINE__)
 * \param[in] note note to print before aborting
 *
 * \return This function called \c bfam_abort() and will not return.
 */
void bfam_abort_verbose(const char *file, int line, ...) BFAM_NORETURN;

/** \c malloc wrapper.
 *
 * This wrapper calls \c malloc and aborts if there is a memory error.
 * The returned pointer needs to be freed by \c bfam_free();
 *
 * \param[in] size allocation size
 *
 * \return pointer to allocated memory.
 */
void *bfam_malloc(size_t size);

/** \c calloc wrapper.
 *
 * This wrapper calls \c calloc and aborts if there is a memory error.
 * The returned pointer needs to be freed by \c bfam_free();
 *
 * \param[in] nmemb number of elements
 * \param[in] size  size of each element
 *
 * \return pointer to allocated memory.
 */
void *bfam_calloc(size_t nmemb, size_t size);

/** \c realloc wrapper.
 *
 * This wrapper calls \c realloc and aborts if there is a memory error.
 * The returned pointer needs to be freed by \c bfam_free();
 *
 * \param[in] ptr  pointer to memory to reallocate
 * \param[in] size allocation size
 *
 * \return pointer to reallocated memory.
 */
void *bfam_realloc(void *ptr, size_t size);

/** \c free wrapper
 *
 * This function frees memory.
 *
 * \param[in,out] ptr pointer to memory to free.
 */
void bfam_free(void *ptr);

/** \c malloc wrapper for cache line aligned memory.
 *
 * This wrapper aligns memory to cache lines.  One needs to call \c
 * bfam_free_aligned() to free the allocated memory.
 *
 * \param[in] size allocation size
 *
 * \return pointer to cache line aligned allocated memory.
 */
void *bfam_malloc_aligned(size_t size);

/** \c free wrapper for cache line aligned memory.
 *
 * This function frees memory that was allocated with \c bfam_malloc_aligned().
 *
 * \param[in,out] ptr pointer to cache line aligned memory to free.
 */
void bfam_free_aligned(void *ptr);

/** Set a signal handler which prints stack traces on terminating signals.
 */
void bfam_signal_handler_set();

// }}}

// {{{ critbit

typedef struct
{
  void *root;
} bfam_critbit0_tree_t;

// }}}

// {{{ dictionary
typedef struct
{
  size_t num_entries;
  bfam_critbit0_tree_t t;
} bfam_dictionary_t;

/** Dictionary to initialize
 *
 * The following function takes a pointer to a dictionary, \a d, and initializes
 * it.
 *
 * \param [out] d pointer to dictionary
 */
void bfam_dictionary_init(bfam_dictionary_t *d);

/** Membership testing.
 *
 * The following function takes a dictionary, \a d, and a \c NULL terminated
 * string, \a u,  and returns non-zero iff \a u in \a d.
 *
 * \param [in] d dictionary
 * \param [in] u possible member
 * \returns non-zero iff \a u in \a d
 */
int bfam_dictionary_contains(bfam_dictionary_t *d, const char *u);

/** Inserting key and value pair where value is a pointer into a dictionary
 *
 * It takes a dictionary, \a d, and possibly mutates it such that a \c NULL
 * terminated strings, \a key and \a value, is a member on exit.
 *
 * \param [in,out] d dictionary
 * \param [in] key possible key
 * \param [in] val possible pointer value
 * \returns:
 *   $\cases{ 0 &if {\rm out of memory} \cr
 *            1 &if {\it key} {\rm was already a member} \cr
 *            2 &if {\it d} {\rm was mutated successfully}}$.
 */
int bfam_dictionary_insert_ptr(bfam_dictionary_t *d, const char *key,
                               const void *val);

/** Inserting key and value pair where value is a \c int into a
 * dictionary.
 *
 * It takes a dictionary, \a d, and possibly mutates it such that a \c NULL
 * terminated strings, \a key and \a value \c snprintf'd, is a member on exit.
 *
 * \param [in,out] d dictionary
 * \param [in] key possible key
 * \param [in] val possible \c int
 * \returns:
 *   $\cases{ 0 &if {\rm out of memory} \cr
 *            1 &if {\it key} {\rm was already a member} \cr
 *            2 &if {\it d} {\rm was mutated successfully}}$.
 */
int bfam_dictionary_insert_int(bfam_dictionary_t *d, const char *key,
                               const int val);

/** Inserting key and value pair where value is a \c bfam_locidx_t into a
 * dictionary.
 *
 * It takes a dictionary, \a d, and possibly mutates it such that a \c NULL
 * terminated strings, \a key and \a value \c snprintf'd, is a member on exit.
 *
 * \param [in,out] d dictionary
 * \param [in] key possible key
 * \param [in] val possible \c bfam_locidx_t
 * \returns:
 *   $\cases{ 0 &if {\rm out of memory} \cr
 *            1 &if {\it key} {\rm was already a member} \cr
 *            2 &if {\it d} {\rm was mutated successfully}}$.
 */
int bfam_dictionary_insert_locidx(bfam_dictionary_t *d, const char *key,
                                  const bfam_locidx_t val);

/** Return a value given a key assuming value is a pointer
 *
 * It takes a dictionary, \a d, returns a pointer that points to where the value
 * pointed associated with a \c NULL terminated \a key.
 *
 * \param [in] d dictionary
 * \param [in] key possible key
 * \returns:
 *   $\cases{ \c NULL & if {\it key} {\rm is not a member} \cr
 *          {\rm pointer to value} & if {\it key} {\rm is a member}$
 */
void *bfam_dictionary_get_value_ptr(bfam_dictionary_t *d, const char *key);

/** Return a value given a key assuming value is a \c bfam_locidx_t.
 *
 * It takes a dictionary, \a d, returns a pointer that points to where the
 * value pointed associated with a \c NULL terminated \a key.
 *
 * \param [in]  d dictionary
 * \param [in]  key possible key
 * \param [out] val value if function returned \c 1
 *
 * \returns:
 *   $\cases{ \c 0 & if {\it key} {\rm is not a member} \cr
 *            \c 1 & if {\it key} {\rm is a member}$
 */
int bfam_dictionary_get_value_locidx(bfam_dictionary_t *d, const char *key,
                                     bfam_locidx_t *val);

/** Clearing a dictionary.
 *
 * Clearing a dictionary (freeing all members) brings us our first code for
 * walking the whole dictionary rather than just tracing a path through it.
 *
 * So, the \c bfam_dictionary_clear function takes a dictionary, \a d, and frees
 * every member of it, mutating the dictionary such that it is empty on exit.
 *
 * \param [in,out] d dictionary
 */
void bfam_dictionary_clear(bfam_dictionary_t *d);

/** Fetching values with a given prefix.
 *
 * The following function takes a dictionary, \a d, and a \c NULL terminated
 * string, \a prefix. Let $S \subseteq d$ where $x \in S$ iff \a prefix is a
 * prefix of \c x, then $\forall x : S.$ \a handle is called with arguments \c x
 * and \c arg.
 * \returns:
 *   $\cases{ 0 &if {\it handle} {\rm returned 0} \cr
 *            1 &if {\rm successful} \cr
 *            2 &if {\it handle} {\rm returned a value} $\notin [0,1]$}$
 * \note (Note that, if |handle| returns 0, the iteration is aborted)
 */
int bfam_dictionary_allprefixed(bfam_dictionary_t *t, const char *prefix,
                                int (*handle)(const char *, const char *,
                                              void *),
                                void *arg);

/** Fetching pointer values with a given prefix.
 *
 * The following function takes a dictionary, \a d, and a \c NULL terminated
 * string, \a prefix. Let $S \subseteq d$ where $x \in S$ iff \a prefix is a
 * prefix of \c x, then $\forall x : S.$ \a handle is called with arguments \c x
 * and \c arg.
 * \returns:
 *   $\cases{ 0 &if {\it handle} {\rm returned 0} \cr
 *            1 &if {\rm successful} \cr
 *            2 &if {\it handle} {\rm returned a value} $\notin [0,1]$}$
 * \note (Note that, if |handle| returns 0, the iteration is aborted)
 *
 * \note The void * input to the handle is the pointer stored in the value, not
 * the pointer to the pointer
 */
int bfam_dictionary_allprefixed_ptr(bfam_dictionary_t *t, const char *prefix,
                                    int (*handle)(const char *, void *, void *),
                                    void *arg);

// }}}

// {{{ logging

#define BFAM_LC_ALL 0
#define BFAM_LC_ROOT 1

#define BFAM_LL_DEFAULT -1
#define BFAM_LL_ALWAYS 0
#define BFAM_LL_TRACE 1
#define BFAM_LL_DEBUG 2
#define BFAM_LL_VERBOSE 3
#define BFAM_LL_INFO 4
#define BFAM_LL_WARNING 5
#define BFAM_LL_ERROR 6
#define BFAM_LL_SILENT 7

/*
 * Setting a hard log threshold is set at compile time
 */
#ifdef BFAM_LOG_LEVEL
#define BFAM_LL_THRESHOLD BFAM_LOG_LEVEL
#else
#ifdef BFAM_DEBUG
#define BFAM_LL_THRESHOLD BFAM_LL_TRACE
#else
#define BFAM_LL_THRESHOLD BFAM_LL_VERBOSE
#endif
#endif

#define BFAM_ROOT_TRACE(...) BFAM_LOG(BFAM_LC_ROOT, BFAM_LL_TRACE, __VA_ARGS__)
#define BFAM_ROOT_LDEBUG(...) BFAM_LOG(BFAM_LC_ROOT, BFAM_LL_DEBUG, __VA_ARGS__)
#define BFAM_ROOT_VERBOSE(...)                                                 \
  BFAM_LOG(BFAM_LC_ROOT, BFAM_LL_VERBOSE, __VA_ARGS__)
#define BFAM_ROOT_INFO(...) BFAM_LOG(BFAM_LC_ROOT, BFAM_LL_INFO, __VA_ARGS__)
#define BFAM_ROOT_WARNING(...)                                                 \
  BFAM_LOG(BFAM_LC_ROOT, BFAM_LL_WARNING, __VA_ARGS__)
#define BFAM_ROOT_LERROR(...) BFAM_LOG(BFAM_LC_ROOT, BFAM_LL_ERROR, __VA_ARGS__)

#define BFAM_TRACE(...) BFAM_LOG(BFAM_LC_ALL, BFAM_LL_TRACE, __VA_ARGS__)
#define BFAM_LDEBUG(...) BFAM_LOG(BFAM_LC_ALL, BFAM_LL_DEBUG, __VA_ARGS__)
#define BFAM_VERBOSE(...) BFAM_LOG(BFAM_LC_ALL, BFAM_LL_VERBOSE, __VA_ARGS__)
#define BFAM_INFO(...) BFAM_LOG(BFAM_LC_ALL, BFAM_LL_INFO, __VA_ARGS__)
#define BFAM_WARNING(...) BFAM_LOG(BFAM_LC_ALL, BFAM_LL_WARNING, __VA_ARGS__)
#define BFAM_LERROR(...) BFAM_LOG(BFAM_LC_ALL, BFAM_LL_ERROR, __VA_ARGS__)

#define BFAM_LOG(category, level, ...)                                         \
  ((level) < BFAM_LL_THRESHOLD ? (void)0 : bfam_log_printf(                    \
                                               __FILE__, __LINE__, (category), \
                                               (level), __VA_ARGS__))

/** Initialization function for the logging system.
 *
 * \param[in] rank      the MPI rank of the current process
 * \param[in] stream    the stream to output the logging information
 * \param[in] threshold the threshold for the logging system (use
 *                      \c BFAM_LL_DEFAULT for the default value or
 *                      \c BFAM_LL_ALWAYS to print all log messages)
 *
 */
void bfam_log_init(int rank, FILE *stream, int threshold);

/** Logging function called by the logging macros.
 *
 * This function is typically not called directly in application code.  Instead
 * one of the logging macros (such as \c BFAM_INFO or \c BFAM_ROOT_INFO is
 * called) will call this function.
 *
 * Messages with \a category \c BFAM_LC_ROOT are only printed by the root rank
 * and \a category \c BFAM_LC_ALL are printed by all ranks.
 *
 * The parameter \a level is used to determine if the message will be printed
 * or not based on the threshold passed into \c bfam_log_init.
 *
 * \param[in] file     file name the abort occurred in (can be obtained from
 *                     \c __FILE__)
 * \param[in] line     line number the abort occurred at (can be obtained from
 *                     \c __LINE__)
 * \param[in] category logging category (e.g., \c BFAM_LC_ROOT or \c
 *                     BFAM_LC_ALL)
 * \param[in] level    logging level (e.g., \c BFAM_LL_TRACE, \c BFAM_LL_ERROR,
 *                     ...)
 * \param[in] ...      It is assumed that the first of the variable arguments
 *                     is a \c printf style format and its required arguments
 *                     follow.
 */
void bfam_log_printf(const char *file, int line, int category, int level, ...);

// }}}

// {{{ utilities

/** Read the contents of a file into a string.
 */
char *bfam_util_read_file(const char *filename, size_t *len);

/** Given a MPI Comm this returns the host rank.
 *
 * The host rank is the rank of the process in the list of processes that
 * have the same MPI processor name.
 */
int bfam_util_get_host_rank(MPI_Comm comm);

// }}}

// {{{ gopt

#define BFAM_GOPT_ONCE 0
#define BFAM_GOPT_REPEAT 1
#define BFAM_GOPT_NOARG 0
#define BFAM_GOPT_ARG 2

/*
 * Moved this from an anonymous struct inside of bfam_gopt_start because
 * clang was barfing.
 */
struct gopt_struct
{
  int k;
  int f;
  const char *s;
  const char *const *l;
};

#define bfam_gopt_start(...)                                                   \
  (const void *)(const struct gopt_struct[])                                   \
  {                                                                            \
    __VA_ARGS__, { 0, 0, NULL, NULL }                                          \
  }
#define bfam_gopt_option(k, f, s, l)                                           \
  {                                                                            \
    k, f, s, l                                                                 \
  }
#define bfam_gopt_shorts(...)                                                  \
  (const char *)(const char[]) { __VA_ARGS__, 0 }
#define bfam_gopt_longs(...)                                                   \
  (const char **)(const char *[]) { __VA_ARGS__, NULL }

/** Sorts options.
 *
 * This function prints to stderr and call exit() on error.
 *
 * \param[in] argc      number of arguments
 * \param[in] argv      array of arguments
 * \param[in] opt_specs options specification generated by \c bfam_gopt_start()
 *
 * \returns a pointer for use in the following calls
 */
void *bfam_gopt_sort(int *argc, const char **argv, const void *opt_specs);

/** Check if option associated with \a key is specified.
 *
 * \param[in] opts options returned from \c bfam_gopt_sort()
 * \param[in] key  option key
 *
 * \returns the number of times the option was specified
 *          which will be 0 or 1 unless \c BFAM_GOPT_REPEAT was used
 */
size_t bfam_gopt(const void *opts, int key);

#if 0
/** Get option associated with \a key and its argument.
 *
 * Writes a pointer to the option argument from the first (or only) occurrence
 * to \c *arg.
 *
 * \param[in]  opts options returned from \c bfam_gopt_sort()
 * \param[in]  key  option key
 * \param[out] arg  pointer to the value of the argument
 *
 * \returns the number of times the option was specified
 */
size_t bfam_gopt_arg(const void *opts, int key, const char **arg);
#endif

#if 0
/** Get \a i'th occurrence of an option.
 *
 * \param[in] opts options returned from \c bfam_gopt_sort()
 * \param[in] key  option key
 * \param[in] i    the occurrence number
 *
 * \returns a pointer to the \a i'th (starting at zero) occurrence of the
 *          option, or \c NULL if it was not specified that many times
 */
const char *bfam_gopt_arg_i(const void *opts, int key, size_t i);
#endif

#if 0
/** Get all the arguments to a repeated option.
 *
 * Writes pointers to the \a option arguments in the order of occurrence to
 * \a args[].  It writes at most \a args_len pointers.
 * If the return value is less than \a args_len, also writes a \c NULL pointer.
 *
 * \param[in]  opts options returned from \c bfam_gopt_sort()
 * \param[in]  key  option key
 * \param[out] args pointers to the values of the multiple arguments
 *
 * \returns the number of times the option was specified
 */
size_t bfam_gopt_args(const void *opts, int key, const char **args,
                      size_t args_len);
#endif

/** Release memory for options.
 *
 * Releases memory allocated in the corresponding call to \c bfam_gopt_sort()
 * \c opts can no longer be used.
 *
 * \param[in] opts options to free
 */
void bfam_gopt_free(void *opts);

// }}}

#ifdef BFAM_USE_LUA
// {{{ lua

#include <lauxlib.h>
#include <lua.h>
#include <lualib.h>

/*
 * Helper function for calling a lua function. Based on generic call function of
 * Listing 25.4-25.6 of
 * @book{Ierusalimschy2006Lua,
 *  author = {Ierusalimschy, Roberto},
 *  title = {Programming in Lua, Second Edition},
 *  year = {2006},
 *  isbn = {8590379825},
 *  publisher = {Lua.Org},
 * }
 */
int bfam_lua_global_function_call(lua_State *L, const char *name,
                                  const char *sig, ...);

/*
 * Some of the following functions were modified from the code at:
 *
 *    http://windrealm.org/tutorials/reading-a-lua-configuration-file-from-c.php
 */

/** Evaluates a Lua expression and returns the boolean result.
 *
 * If an error occurs or the result is not a boolean, def is returned.
 */
int bfam_lua_expr_boolean(lua_State *L, const char *expr, int def);

/** Evaluates a Lua expression and returns the integer result.
 *
 * If an error occurs or the result is not a integer, def is returned.
 */
lua_Integer bfam_lua_expr_integer(lua_State *L, const char *expr,
                                  lua_Integer def);

/** Evaluates a Lua expression and returns the number result.
 *
 * If an error occurs or the result is not a number, def is returned.
 */
lua_Number bfam_lua_expr_number(lua_State *L, const char *expr, lua_Number def);

/** Evaluates a Lua expression and returns the string result.
 *
 * If an error occurs or the result is not a string, a copy of def is returned.
 * The user of this function is responsible for cleaning up the memory.
 *
 */
char *bfam_lua_expr_string(lua_State *L, const char *expr, const char *def);

// }}} lua
#endif

// {{{ subdomain

struct bfam_subdomain;

/*
 * This is the field initialization function which will be called for each
 * subdomain.  This function fills the \a npoints values of \a field.
 *
 * \param [in]  npoints   this is the length of x, y, z, and field
 * \param [in]  x         x-coordinates of the points
 * \param [in]  y         y-coordinates of the points
 * \param [in]  z         z-coordinates of the points
 * \param [in]  s         pointer to the subdomain that the field is in
 * \param [in]  arg       user pointer
 * \param [out] field     the field values that need to be set
 *
 */
typedef void (*bfam_subdomain_init_field_t)(
    bfam_locidx_t npoints, const char *name, bfam_real_t time,
    bfam_real_t *restrict x, bfam_real_t *restrict y, bfam_real_t *restrict z,
    struct bfam_subdomain *s, void *arg, bfam_real_t *restrict field);

/**
 * base structure for to store glue data for a subdomain, i.e., plus and minus
 * sides stuff. It should also be included as first member with the name base
 * \code{.c}
 * typedef struct new_subdomain_glue_data_type
 * {
 *   bfam_subdomain_t base;
 *   ...
 * }
 */
typedef struct bfam_subdomain_glue_data
{
  bfam_locidx_t rank; /* Rank of the subdomain on this side */
  bfam_locidx_t id_s; /* Sort Id of the subdomain on this side */
  bfam_locidx_t id;   /* Id of the subdomain on this side */

  bfam_dictionary_t fields; /**< a dictionary storing glue fields */

  bfam_critbit0_tree_t tags; /**< critbit for tags for the glue */

  /* The following pointers should only be \ne NULL on the minus side */
  struct bfam_subdomain *sub_m; /* Local neighboring subdomain */
} bfam_subdomain_glue_data_t;

/**
 * base structure for all subdomains types. Any new subdomain should have this
 * as its first member with the name base, i.e.,
 * \code{.c}
 * typedef struct new_subdomain_type
 * {
 *   bfam_subdomain_t base;
 *   ...
 * }
 */
typedef struct bfam_subdomain
{
  bfam_locidx_t id;
  bfam_locidx_t uid;         /**< typically physics id */
  char *name;                /**< Name of the subdomain */
  bfam_critbit0_tree_t tags; /**< critbit for tags for the subdomain */
  bfam_dictionary_t fields;  /**< a dictionary storing pointers to fields */

  bfam_dictionary_t fields_face; /**< a dictionary storing face fields */

  /* glue quantities */
  bfam_subdomain_glue_data_t *glue_m;
  bfam_subdomain_glue_data_t *glue_p;

  /* user data */
  void *user_data;

  /* Function pointers that domain will need to call */
  void (*free)(struct bfam_subdomain *thisSubdomain);

  /**< Write a vtk vtu file */
  int (*vtk_write_vtu_piece)(struct bfam_subdomain *thisSubdomain, FILE *file,
                             bfam_real_t time, const char **scalars,
                             const char **vectors, const char **components,
                             int writeBinary, int writeCompressed, int rank,
                             bfam_locidx_t id, int Np_write);

  /**< Add a field to the subdomain */
  int (*field_add)(struct bfam_subdomain *thisSubdomain, const char *name);

  /**< Glue grid communication info */
  void (*glue_comm_info)(struct bfam_subdomain *thisSubdomain, int *rank,
                         bfam_locidx_t *s, int num_sort, size_t *send_sz,
                         size_t *recv_sz, void *args);

} bfam_subdomain_t;

/** Add a tag to the subdomain
 *
 * \param [in,out] thisSubdomain subdomain to andd the tag to
 * \param [in]     tag           tag of the domain (\0 terminated string)
 *
 */
void bfam_subdomain_add_tag(bfam_subdomain_t *thisSubdomain, const char *tag);

/** Check to see if a subdomain has a tag
 *
 * \param [in,out] thisSubdomain subdomain to search for the tag
 * \param [in]     tag           tag of the domain (\0 terminated string)
 *
 * \return nonzero iff \a thisSubdomain has the tag \a tag
 */
int bfam_subdomain_has_tag(bfam_subdomain_t *thisSubdomain, const char *tag);

// }}}

// {{{ domain

/**
 * structure containing a domain (which is a collection of subdomains!)
 */
typedef struct bfam_domain
{
  bfam_subdomain_t **subdomains; /**< array of pointers to subdomains */
  bfam_locidx_t num_subdomains;  /**< number of subdomains that are
                                     currently in the domain */
  bfam_locidx_t sizeSubdomains;  /**< total number of subdomains the domain
                                      can hold, i.e.  size of the array*/
  MPI_Comm comm;                 /**< communicator for the whole domain */
  bfam_dictionary_t name2num;    /**< dictionary map for convertings
                                      subdomain names to numbers */
} bfam_domain_t;

typedef enum bfam_domain_match {
  BFAM_DOMAIN_AND,
  BFAM_DOMAIN_OR,
} bfam_domain_match_t;

/** initializes a domain
 *
 * \param [in,out] domain pointer to the domain
 * \param [in]     domComm pointer to the communicator for the domain
 */
void bfam_domain_init(bfam_domain_t *domain, MPI_Comm domComm);

/** Get subdomain
 *
 * \param [in]     id Id of subdomain to get
 *
 * \return subdomain pointer
 */
bfam_subdomain_t *bfam_domain_get_subdomain_by_num(bfam_domain_t *thisDomain,
                                                   bfam_locidx_t id);

/** Get subdomains with tags passed in
 *
 * \param [in]  thisDomain    domain to search for subdomains in
 * \param [in]  matchType     type of match, \c BFAM_DOMAIN_OR will
 *                            match subdomains with any of the tags
 *                            and \c BFAM_DOMAIN_AND will match subdomains
 *                            with all of the tags.
 * \param [in]  tags          \c NULL terminated array of the tags to match
 * \param [in]  numEntries    number of entries in the \a subdomains array
 * \param [out] subdomains    array of pointers to be filled with matching
 *                            subdomains
 * \param [out] num_subdomains number of matching subdomains
 *
 */
void bfam_domain_get_subdomains(bfam_domain_t *thisDomain,
                                bfam_domain_match_t match, const char **tags,
                                bfam_locidx_t numEntries,
                                bfam_subdomain_t **subdomains,
                                bfam_locidx_t *num_subdomains);

/** Add tag to subdomains matching the tags passed in.
 *
 * \param [in,out]  thisDomain    domain to search for subdomains in
 * \param [in]      matchType     type of match, \c BFAM_DOMAIN_OR will
 *                                match subdomains with any of the tags
 *                                and \c BFAM_DOMAIN_AND will match subdomains
 *                                with all of the tags.
 * \param [in]      mtags         \c NULL terminated array of the tags to match
 * \param [in]      tag           tag to add to the subdomains
 *
 */
void bfam_domain_add_tag(bfam_domain_t *thisDomain, bfam_domain_match_t match,
                         const char **mtags, const char *tag);

/** Add fields to subdomains matching the tags passed in.
 *
 * \param [in,out]  thisDomain    domain to search for subdomains in
 * \param [in]      matchType     type of match, \c BFAM_DOMAIN_OR will
 *                                match subdomains with any of the tags
 *                                and \c BFAM_DOMAIN_AND will match subdomains
 *                                with all of the tags.
 * \param [in]      tags          \c NULL terminated array of the tags to match
 * \param [in]      fields        fields to add to the subdomains
 *
 */
void bfam_domain_add_fields(bfam_domain_t *thisDomain,
                            bfam_domain_match_t match, const char **tags,
                            const char **fields);
// }}}

// {{{ communicator
/**
 * structure for storing subdomain specific data for the communicator
 */
typedef struct bfam_comm_subdata
{
  bfam_subdomain_t *subdomain; /**< pointer to my local subdomain */

  size_t send_sz;     /**< amount of data this subdomain can send */
  void *send_buf;     /**< pointer for local send buffer */
  size_t send_offset; /**< offset into the global send buffer */

  size_t recv_sz;     /**< amount of data this subdomain should receive */
  void *recv_buf;     /**< pointer for local recv buffer */
  size_t recv_offset; /**< offset into the global receive buffer */
} bfam_comm_subdata_t;

/**
 * structure for storing processor specific data for the communicator
 */
typedef struct bfam_comm_procdata
{
  bfam_locidx_t rank; /**< neighboring processors rank */

  size_t send_sz; /**< amount to send */
  void *send_buf; /**< pointer to send buffer */

  size_t recv_sz; /**< amount to recv */
  void *recv_buf; /**< pointer to recv buffer */
} bfam_comm_procdata_t;

/**
 * structure for doing communication
 */
typedef struct bfam_communicator
{
  MPI_Comm comm;          /**< communicator used for communication */
  int tag;                /**< user specified tag for this communicator */
  bfam_locidx_t num_subs; /**< number of subdomains in the communicator */

  bfam_locidx_t num_procs; /**< number of processors in the communicator */

  MPI_Request *send_request; /**< send request */
  MPI_Request *recv_request; /**< recv request */

  MPI_Status *send_status; /**< send status */
  MPI_Status *recv_status; /**< recv status */

  void *send_buf; /**< full send buffer */
  void *recv_buf; /**< full recv buffer */

  size_t send_sz; /**< full send size */
  size_t recv_sz; /**< full recz size */

  bfam_comm_procdata_t *proc_data; /**< array of structure with neighboring
                                        processor data */
  bfam_comm_subdata_t *sub_data;   /**< array of structure with subdomains
                                        specific information */

  void *user_args; /**< user custom data to pass through */
} bfam_communicator_t;

/** create a communicator
 *
 * \param [in] domain     domain to output to communicate
 * \param [in] match      type of match, \c BFAM_DOMAIN_OR will
 *                        match subdomains with any of the tags
 *                        and \c BFAM_DOMAIN_AND will match subdomains
 *                        with all of the tags.
 * \param [in] tags       \c NULL terminated array of the tags to match
 *                        glue grids doing the communication
 * \param [in] comm       MPI communicator
 * \param [in] tag        user specified communicator tag
 * \param [in] userdata   user custom data to pass through
 *
 * \return the newly created communicator
 */
bfam_communicator_t *bfam_communicator_new(bfam_domain_t *domain,
                                           bfam_domain_match_t match,
                                           const char **tags, MPI_Comm comm,
                                           int tag, void *user_data);

/** Clean up communicator
 *
 * frees any memory allocated by the communicator
 *
 * \param [in,out] communicator communicator to clean up
 */
void bfam_communicator_free(bfam_communicator_t *communicator);

// }}}

// {{{ domain pxest

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wcast-align"
#endif
#if BFAM_DGX_DIMENSION == 2
#include <p4est.h>
#include <p4est_bits.h>
#include <p4est_connectivity.h>
#include <p4est_extended.h>
#include <p4est_ghost.h>
#include <p4est_iterate.h>
#include <p4est_lnodes.h>
#include <p4est_mesh.h>
#include <p4est_nodes.h>
#include <p4est_vtk.h>
#elif BFAM_DGX_DIMENSION == 3
#include <p4est_to_p8est.h>
#include <p8est.h>
#include <p8est_bits.h>
#include <p8est_connectivity.h>
#include <p8est_extended.h>
#include <p8est_ghost.h>
#include <p8est_iterate.h>
#include <p8est_lnodes.h>
#include <p8est_mesh.h>
#include <p8est_nodes.h>
#include <p8est_vtk.h>
#else
#error "bad dimension"
#endif
#if defined(__clang__)
#pragma clang diagnostic pop
#endif

/*  This is the user data each p4est qudrant will hold.
 */
typedef struct
{
  uint8_t flags;
  int8_t N;
  int8_t Nold;
  bfam_locidx_t subd_id;
  bfam_locidx_t elem_id;
  bfam_locidx_t root_id; /* used for element merging */
#if BFAM_DGX_DIMENSION == 2
  bfam_locidx_t glue_id[4]; /* used for element merging */
#elif BFAM_DGX_DIMENSION == 3
  bfam_locidx_t glue_id[6]; /* used for element merging */
#else
#error "bad dimension"
#endif
} bfam_pxest_user_data_t;

/**
 * structure containing a domain managed by p4est
 */
typedef struct bfam_domain_pxest
{
  bfam_domain_t base;         /** parent domain */
  bfam_dictionary_t *N2N;     /** Dictionary of projection operators */
  p4est_connectivity_t *conn; /** connectivity for p4est */
  p4est_t *pxest;             /** forest of quadtrees */
  bfam_dictionary_t *dgx_ops; /** Dictionary of dgx operators operators */
} bfam_domain_pxest_t;

typedef struct
{
  uint8_t *dst_to_adapt_flags;
  int8_t *dst_to_dst_chld_id;
  bfam_locidx_t *dst_to_src_subd_id;
  bfam_locidx_t *dst_to_src_elem_id;
  bfam_locidx_t num_dst;

  int8_t *coarse_dst_to_src_chld_id;
  bfam_locidx_t *coarse_dst_to_src_subd_id;
  bfam_locidx_t *coarse_dst_to_src_elem_id;
  bfam_locidx_t num_coarse_dst;
} bfam_domain_pxest_transfer_maps_t;

/** create a pxest managed domain
 *
 * \warning It is the callers responsibility to ensure that
 *          \a domComm and \a conn are freed after this domain is.
 *
 * \param [in]  domComm         pointer to the communicator for the domain
 * \param [in]  conn            pointer to the pxest connectivity for the domain
 * \param [in]  min_quadrants   Minimum initial quadrants per processor.
 *                              Makes the refinement pattern mpisize-specific.
 * \param [in]  min_level       The forest is refined at least to this level.
 *                              May be negative or 0, then it has no effect.
 * \param [in]  fill_uniform    If true, fill the forest with a uniform mesh
 *                              instead of the coarsest possible one.
 *                              The latter is partition-specific so that
 *                              is usually not a good idea.
*
 * \return the newly created pxest managed domain
 */
bfam_domain_pxest_t *bfam_domain_pxest_new_ext(MPI_Comm domComm,
                                               p4est_connectivity_t *conn,
                                               p4est_locidx_t min_quadrants,
                                               int min_level, int fill_uniform);

/** Clean up domain
 *
 * frees any memory allocated by the domain and calls free command on all
 * subdomains
 *
 * \param [in,out] domain domain to clean up
 */
void bfam_domain_pxest_free(bfam_domain_pxest_t *domain);

/** Fill a \c glueID based on tree ids.
 *
 * This fills a \c glueID array for the quadrants based on glue ids given for
 * the trees.
 *
 * \param [in]  pxest        pointer to a p4est
 * \param [in]  tree_to_glue an array indicating the glue ids for faces in the
 *                           connectivity structure of the p4est mesh.  The ids
 *                           for the tree faces are stored in -x +x -y +y -z +z
 *                           order for each tree with the face index running
 *                           the fastest.
 * \param [out] quad_to_glue an array indicating the glue ids for faces in the
 *                            p4est mesh.  The ids for the quadrant faces are
 *                            stored in -x +x -y +y -z +z order for each
 *                            quadrant with the face index running the fastest.
 */
void bfam_domain_pxest_quad_to_glueid(p4est_t *pxest,
                                      const bfam_locidx_t *tree_to_glueid,
                                      bfam_locidx_t *quad_to_glueid);

/** Takes an initialized domain and generates a DG hex mesh
 *
 * \param [in,out] domain        pointer to the initialized pxest managed
 *                               domain
 * \param [in]     num_subdomains number of volume subdomains to generate
 * \param [in]     subdomainID   array of length \c pxest->local_num_quadrants
 *                               which indicates the subdomain id for each
 *                               element
 * \param [in]     roots         array of roots for the volume subdomains
 * \param [in]     N             array of orders for the volume subdomains
 * \param [in]     glueID        array indicating what glue subdomains exist.
 *                               If the number is negative than no glue
 *                               subdomain will be created and if the number
 *                               is positive a subdomain will be created.
 *                               This is of length \c
 *                               pxest->local_num_quadrants*NumberOfFaces
 *                               if \c NULL it will be ignored.
 * \param [in] glue_order        user callback function to allow the user to
 *                               set the order of the glue grids
 * \param [in] go_user_args      user argument for glue_order
 */

typedef int (*bfam_glue_order_t)(const int N_m, const int N_p,
                                 const bfam_locidx_t uid, void *user_args);

void bfam_domain_pxest_split_dgx_subdomains(
    bfam_domain_pxest_t *domain, bfam_locidx_t num_subdomains,
    bfam_locidx_t *subdomainID, bfam_locidx_t *roots, int *N,
    bfam_locidx_t *glueID, bfam_glue_order_t glue_order, void *go_user_args);

/** Given a domain and a refined p4est generates a subdomain spliting
 *
 * \param [in]     pxest         pointer to a p4est
 * \param [in]     pflags        flags on what type of p adaptivity
 * \param [out]    num_subdomains number of volume subdomains to generate
 * \param [in]     subdomainID   array of length \c pxest->local_num_quadrants
 *                               which indicates the subdomain id for each
 *                               element
 * \param [in]     roots         array of roots for the volume subdomains
 * \param [in]     N             array of orders for the volume subdomains
 * \param [in]     glueID        array indicating what glue subdomains exist.
 *                               If the number is negative than no glue
 *                               subdomain will be created and if the number
 *                               is positive a subdomain will be created.
 *                               This is of length \c
 *                               pxest->local_num_quadrants*NumberOfFaces
 *                               if \c NULL it will be ignored.
 */
void bfam_domain_pxest_compute_split(p4est_t *pxest, uint8_t pflags,
                                     bfam_locidx_t *num_subdomains,
                                     bfam_locidx_t **subdomain_id,
                                     bfam_locidx_t **roots, int **N,
                                     bfam_locidx_t **glue_id);

/** Mark pxest quadrants with refinement info
 *
 * \param [in,out] domain        pointer to the initialized pxest managed
 *                               domain
 */
void bfam_domain_pxest_mark_elements(bfam_domain_pxest_t *domain);

/** Generate mesh transfer maps
 *
 * \param [in,out] maps     pointer to the transfer maps to be filled
 * \param [in] domain_dst   destination domain
 * \param [in] domain_src   source domain
 *
 */
void bfam_domain_pxest_transfer_maps_init(
    bfam_domain_pxest_transfer_maps_t *maps, bfam_domain_pxest_t *domain_dst,
    bfam_domain_pxest_t *domain_src);
/** Generate mesh transfer maps
 *
 * \param [in,out] maps     pointer to the transfer maps to be freed
 */
void bfam_domain_pxest_transfer_maps_free(
    bfam_domain_pxest_transfer_maps_t *maps);

/* Callbacks for pxest quadrants */
int bfam_domain_pxest_quadrant_coarsen(p4est_t *p4est,
                                       p4est_topidx_t which_tree,
                                       p4est_quadrant_t *quadrants[]);
int bfam_domain_pxest_quadrant_refine(p4est_t *p4est, p4est_topidx_t which_tree,
                                      p4est_quadrant_t *quadrant);
void bfam_domain_pxest_quadrant_init(p4est_t *p4est, p4est_topidx_t which_tree,
                                     p4est_quadrant_t *quadrant);
void bfam_domain_pxest_quadrant_replace(p4est_t *p4est,
                                        p4est_topidx_t which_tree,
                                        int num_outgoing,
                                        p4est_quadrant_t *outgoing[],
                                        int num_incoming,
                                        p4est_quadrant_t *incoming[]);

// }}}

// {{{ subdomain dgx

/*
 * Just like in p4est the quads are stored in z-order so the
 * face and corner orderings are:
 *
 *            2           3
 *            +-----------+
 *            |     3     |
 *            |           |
 *            | 0       1 |
 *            |           |
 *            |     2     |
 *            +-----------+
 *            0           1
 *
 * This subdomain is composed of nodal (tensor product LGL points)
 * The naming convention follows that of:
 *
 * @book{Hesthaven:2007:NDG:1557392,
 *   author = {Hesthaven, Jan S. and Warburton, Tim},
 *   title = {Nodal Discontinuous {Galerkin} Methods: Algorithms, Analysis, and
 *            Applications},
 *   year = {2007},
 *   isbn = {0387720650, 9780387720654},
 *   edition = {1st},
 *   publisher = {Springer Publishing Company, Incorporated},
 * }
 *
 * where it makes sense.
 *
 */

struct bfam_subdomain_dgx;

typedef struct bfam_subdomain_dgx_glue_data
{
  bfam_subdomain_glue_data_t base;

  int same_order; /* Boolean indicating if the elements connected
                   * through the glue are of the same order.
                   */

  bfam_locidx_t *EToEp;  /* Element     number on connected subdomain */
  int8_t *EToHp;         /* Hanging     number on connected subdomain */
  int8_t *EToOp;         /* Orientation number on connected subdomain */
  bfam_locidx_t *EToEm;  /* Element     number on local subdomain */
  int8_t *EToFm;         /* Face        number on local subdomain */
  int8_t *EToHm;         /* Hanging     number on local subdomain */
  bfam_locidx_t **mapOp; /* mapping the orientation */
  int num_orient;        /* number of orientations */

  /* The following pointers should only be \ne NULL on the minus side */
  int num_interp; /* number of interpolation operators */

  bfam_real_t **interpolation; /* array of interpolation operators;
                                * the first is for non-hanging faces
                                * the rest are for the hanging faces;
                                * if the operator is NULL it is assumed
                                * to be the identity operator.
                                */

  bfam_real_t **projection; /* array of projection operators;
                             * the first is for non-hanging faces
                             * the rest are for the hanging faces;
                             * if the operator is NULL it is assumed
                             * to be the identity operator.
                             */

  bfam_real_t **massprojection; /* array of mass projection operators;
                                 * the first is for non-hanging faces
                                 * the rest are for the hanging faces;
                                 * if the operator is NULL it is assumed
                                 * to be the identity operator.
                                 */

  bfam_real_t *exact_mass; /* exact mass matrix for this grid */

} bfam_subdomain_dgx_glue_data_t;

typedef struct bfam_subdomain_dgx
{
  bfam_subdomain_t base;

  int dim; /* dimensionality of this dgx */

  int N;    /* 1D Polynomial Order */
  int Np;   /* Number of points in the element */
  int *Ngp; /* Number of geometry points */

  int numg; /* number of geometry types */
  int *Ng;  /* geometry based quantities:
             * faces, edges, corners, etc.
             */
  // JK int              Nh;         /* Number of interpolations to glue */
  // JK int              No;         /* Number of orientations */

  /* These pointers are NOT `owned' by the subdomain, and thus not freed in
   * cleanup
   */
  bfam_real_t *r;        /* 1D LGL Nodal Point in [-1,1] */
  bfam_real_t *w;        /* 1D LGL Weights */
  bfam_real_t *wi;       /* inverse of 1D LGL Weights */
  bfam_real_t *Dr;       /* 1D LGL differentiation matrix */
  bfam_long_real_t *lr;  /* long format 1D LGL Nodal Point in [-1,1] */
  bfam_long_real_t *lw;  /* long format 1D LGL Nodal weights in [-1,1] */
  bfam_long_real_t *lDr; /* long format 1D LGL differentiation matrix */
  bfam_long_real_t *lV;  /* 1D Vandermonde matrix for this N */
  /* end of un-owned pointers */

  bfam_locidx_t K; /* Number of elements in the subdomain */

  bfam_locidx_t *vmapM; /* Mapping into the volume for the minus side of
                           the face mesh */
  bfam_locidx_t *vmapP; /* Mapping into the volume for the plus  side of
                           the face mesh */

  bfam_locidx_t *EToQ; /* Element to p4est local quadrant map */

  int ***gmask; /* geometry mask: same order as Ng */

  uint8_t *hadapt; /* length K where entries indicate h-adaptation
                    * which has the flags:
                    *   BFAM_FLAG_COARSEN
                    *   BFAM_FLAG_REFINE
                    *   BFAM_FLAG_ADAPTED
                    */

  int8_t *padapt; /* length K where entries indicate the desired
                   * order for each element after p-adaptation
                   */

  int8_t *lvl; /* lenght K where entries indicate current level */

  bfam_locidx_t *q_id; /* lenght K where entries are used by the subdomain
                          generator */
} bfam_subdomain_dgx_t;

// }}}

// {{{ vtk
/** Write out vtk files for each domain.
 *
 * This loops through the subdomains and writes out a vtk file for each one.
 *
 * \param [in] domain     domain to output to vtk files
 * \param [in] match      type of match, \c BFAM_DOMAIN_OR will
 *                        match subdomains with any of the tags
 *                        and \c BFAM_DOMAIN_AND will match subdomains
 *                        with all of the tags.
 * \param [in] tags       \c NULL terminated array of the tags to match
 * \param [in] prefix     prefix for the vtk files
 * \param [in] scalars    \c NULL terminated array of scalars to match
 * \param [in] vectors    \c NULL terminated array of vector that will
 *                        be outputted
 * \param [in] components \c NULL terminated array of vector components
 * \param [in] binary     boolean to indicate if the data should be
 *                        written in binary
 * \param [in] compress   boolean to indicate if the data should be compressed
 *                        when writing a binary data
 * \param [in] Np_write   Number of points to write out in each dimension of
 *                        element (poorman's high-order viz): Np_write = 0 means
 *                        use element Np
 *
 *
 * \note That a domain will only output scalars and vectors that it contains.
 *
 */

void bfam_vtk_write_file(bfam_domain_t *domain, bfam_domain_match_t match,
                         const char **tags, const char *directory,
                         const char *prefix, bfam_real_t time,
                         const char **scalars, const char **vectors,
                         const char **components, int binary, int compress,
                         int Np_write);

// }}}

// {{{ pcg32
// *Really* minimal PCG32 code / (c) 2014 M.E. O'Neill / pcg-random.org
// Licensed under Apache License 2.0 (NO WARRANTY, etc. see website)
// Copied from: http://www.pcg-random.org/
// Copied from: https://github.com/imneme/pcg-c-basic

typedef struct
{
  uint64_t state;
  uint64_t inc;
} bfam_pcg32_random_t;

void bfam_pcg32_srandom_r(bfam_pcg32_random_t *rng, uint64_t initstate,
                          uint64_t initseq);

uint32_t bfam_pcg32_boundedrand_r(bfam_pcg32_random_t *rng, uint32_t bound);
// }}}

// {{{ bfamo
#ifdef BFAM_USE_BFAMO
#include <bfam.h>
#include <occa_c.h>

#if defined(KERNEL_TIME) || defined(COMPUTE_STATS)
#define STATS
#endif

#ifdef STATS
#include <sc_statistics.h>
#endif

/* TODO: Add bfamo_locidx_t for occa */
/* TODO: Determine why SINGLE doesn't work */

#ifdef BFAMO_REAL_DOUBLE
typedef double metric_real_t;
#define occaMetric occaDouble
typedef double bfamo_real_t;
#define occaReal occaDouble
#define BFAMO_REAL_MAX DBL_MAX
#define BFAMO_REAL(x) (x)
#define BFAMO_REAL_FMTe "24.16e"
#define BFAMO_REAL_MPI MPI_DOUBLE
#define BFAMO_REAL_SQRT sqrt
#else
typedef float bfamo_real_t;
#define occaMetric occaDouble
typedef double metric_real_t;
#define occaReal occaFloat
#define BFAMO_REAL_MAX FLT_MAX
#define BFAMO_REAL(x) BFAM_APPEND(x, f)
#define BFAMO_REAL_FMTe "24.16e"
#define BFAMO_REAL_MPI MPI_FLOAT
#define BFAMO_REAL_SQRT sqrtf
#endif

#define MAX_COMMAND_LEN 2048
#if BFAM_DGX_DIMENSION == 2
#include <p4est_iterate.h>
#elif BFAM_DGX_DIMENSION == 3
#include <p8est_iterate.h>
#define P4EST_CONNECT_CORNER P8EST_CONNECT_CORNER
#else
#error "bad dimension"
#endif

#if BFAM_DGX_DIMENSION == 2
#define D3_AP(A1, A2) (A1)
#define D3_OP(A) BFAM_NOOP()

#define BFAMO_PXEST_CONNECT P4EST_CONNECT_FULL
#define BFAMO_DIM 2

#define ID_VGEO_JR0X0 0
#define ID_VGEO_JR0X1 1
#define ID_VGEO_JR1X0 2
#define ID_VGEO_JR1X1 3
#define ID_VGEO_JINV 4
#define ID_VGEO_W 5
#define ID_VGEO_X0 6
#define ID_VGEO_X1 7
#define NVGEO 8

#define ID_SGEO_NX0 0
#define ID_SGEO_NX1 1
#define ID_SGEO_SJWJ 2
#define NSGEO 3

#define ID_GGEO_NX0 0
#define ID_GGEO_NX1 1
#define ID_GGEO_SJW 2
#define ID_GGEO_X0 3
#define ID_GGEO_X1 4
#define NGGEO 5

#elif BFAM_DGX_DIMENSION == 3
#define D3_AP(A1, A2) (A1 A2)
#define D3_OP(A) A

#define BFAMO_PXEST_CONNECT P8EST_CONNECT_FULL
#define BFAMO_DIM 3

#define ID_VGEO_JR0X0 0
#define ID_VGEO_JR0X1 1
#define ID_VGEO_JR0X2 2
#define ID_VGEO_JR1X0 3
#define ID_VGEO_JR1X1 4
#define ID_VGEO_JR1X2 5
#define ID_VGEO_JINV 6
#define ID_VGEO_JR2X0 7
#define ID_VGEO_JR2X1 8
#define ID_VGEO_JR2X2 9
#define ID_VGEO_JINV2 10
#define ID_VGEO_W 11
#define ID_VGEO_X0 12
#define ID_VGEO_X1 13
#define ID_VGEO_X2 14
#define NVGEO 15

#define ID_SGEO_NX0 0
#define ID_SGEO_NX1 1
#define ID_SGEO_NX2 2
#define ID_SGEO_SJWJ 3
#define NSGEO 4

#define ID_GGEO_NX0 0
#define ID_GGEO_NX1 1
#define ID_GGEO_NX2 2
#define ID_GGEO_SJW 3
#define ID_GGEO_X0 4
#define ID_GGEO_X1 5
#define ID_GGEO_X2 6
#define NGGEO 7

#else
#error "bad dimension"
#endif

#define ID_QUAD_TID 0

#define ID_QUAD_LVL 1
#define ID_QUAD_X 2
#define ID_QUAD_Y 3
#define ID_QUAD_Z 4
#define NQUAD 5

/* OCCA modes */
#define SERIAL 0
#define OPENMP 1
#define OPENCL 2
#define CUDA 3
extern const char *occa_modes[5];

/* prefix for the projection dictionary */
extern const char *prj_base_names[6];

/* coarse, same, refine order */
#define BFAMO_ADAPT_PR 0
#define BFAMO_ADAPT_PS 1
#define BFAMO_ADAPT_PC 2
#define NBFAMO_ADAPT 3
#define BFAMO_ADAPT_P_CHANGE_MIN -1
extern const char *bfamo_adapt_postfix[4];

#define ID_FIELD(id, vec, n, e, Np, K)                                         \
  (((id) % (vec)) + (n) * (vec) + (e) * (Np) * (vec) +                         \
   ((id) / (vec)) * (K) * (Np) * (vec))
#define NFIELDS_VEC 4

#define BFAMO_DICTIONARY_INSERT_PTR(a, b, c)                                   \
  BFAM_ABORT_IF(bfam_dictionary_insert_ptr((a), (b), (c)) == 1,                \
                "Insert pointer error Error")

#define GB (1e9)
#define GiB (1024 * 1024 * 1024)

#ifdef COMPUTE_STATS
#define BFAMO_TIME_BLOCK(block, name)                                          \
  do                                                                           \
  {                                                                            \
    occaStream bfamo_xxx_curent_stream_xxx;                                    \
    bfamo_xxx_curent_stream_xxx = occaDeviceGetStream((_t.bdevice->device));   \
    occaDeviceSetStream((_t.bdevice->device), (_t.cmdx));                      \
    occaDeviceFinish((_t.bdevice->device));                                    \
    occaDeviceSetStream((_t.bdevice->device), (_t.copy));                      \
    occaDeviceFinish((_t.bdevice->device));                                    \
    occaDeviceSetStream((_t.bdevice->device), bfamo_xxx_curent_stream_xxx);    \
    occaGetStreamFree(bfamo_xxx_curent_stream_xxx);                            \
    double bfamo_xxx_time_xxx = MPI_Wtime();                                   \
    {                                                                          \
      block;                                                                   \
    }                                                                          \
    bfamo_xxx_curent_stream_xxx = occaDeviceGetStream((_t.bdevice->device));   \
    occaDeviceSetStream((_t.bdevice->device), (_t.cmdx));                      \
    occaDeviceFinish((_t.bdevice->device));                                    \
    occaDeviceSetStream((_t.bdevice->device), (_t.copy));                      \
    occaDeviceFinish((_t.bdevice->device));                                    \
    occaDeviceSetStream((_t.bdevice->device), bfamo_xxx_curent_stream_xxx);    \
    occaGetStreamFree(bfamo_xxx_curent_stream_xxx);                            \
    bfamo_xxx_time_xxx = MPI_Wtime() - bfamo_xxx_time_xxx;                     \
    sc_statistics_accumulate(_t.stat, (name), bfamo_xxx_time_xxx * 1e3);       \
  } while (0)

#else
#define BFAMO_TIME_BLOCK(block, name)                                          \
  do                                                                           \
  {                                                                            \
    block;                                                                     \
  } while (0)
#endif

#ifdef KERNEL_TIME

#define BFAMO_KERNEL_RUN(bsk, ...)                                             \
  do                                                                           \
  {                                                                            \
    occaDeviceFinish(_t.bdevice->device);                                      \
                                                                               \
    double bfamo_xxx_k_time_xxx = MPI_Wtime();                                 \
    bfamo_kernel_set_working_dims(bsk.bk, bsk.dim, bsk.items, bsk.groups);     \
    occaKernelRun(bsk.bk->kernel, __VA_ARGS__);                                \
    occaDeviceFinish(_t.bdevice->device);                                      \
    bfamo_xxx_k_time_xxx = MPI_Wtime() - bfamo_xxx_k_time_xxx;                 \
    sc_statistics_accumulate(_t.stat, bsk.bk->stat_key,                        \
                             bfamo_xxx_k_time_xxx * 1e3);                      \
                                                                               \
    int bfamo_xxx_flops1_xxx = bsk.bk->flops1;                                 \
    int bfamo_xxx_flops2_xxx = bsk.bk->flops2;                                 \
    if (bfamo_xxx_flops1_xxx + bfamo_xxx_flops2_xxx > 0)                       \
    {                                                                          \
      double bfamo_xxx_flops_xxx =                                             \
          (bfamo_xxx_flops1_xxx + bfamo_xxx_flops2_xxx * bsk.work2) *          \
          bsk.work1 / (bfamo_xxx_k_time_xxx * GB);                             \
      sc_statistics_accumulate(_t.stat, bsk.bk->flops_key,                     \
                               bfamo_xxx_flops_xxx);                           \
    }                                                                          \
                                                                               \
    int bfamo_xxx_memread_xxx = bsk.bk->memread;                               \
    if (bfamo_xxx_memread_xxx > 0)                                             \
    {                                                                          \
      double bfamo_xxx_bandwidth_xxx = bfamo_xxx_memread_xxx * bsk.work1 *     \
                                       ((double)sizeof(bfamo_real_t)) /        \
                                       (GB * bfamo_xxx_k_time_xxx);            \
      sc_statistics_accumulate(_t.stat, bsk.bk->memread_key,                   \
                               bfamo_xxx_bandwidth_xxx);                       \
    }                                                                          \
  } while (0)

#else

#define BFAMO_KERNEL_RUN(bsk, ...)                                             \
  do                                                                           \
  {                                                                            \
    bfamo_kernel_set_working_dims(bsk.bk, bsk.dim, bsk.items, bsk.groups);     \
    occaKernelRun(bsk.bk->kernel, __VA_ARGS__);                                \
  } while (0)

#endif

typedef struct
{
  occaDevice device;
  double mem_fudge;
} bfamo_device_t;

#ifdef STATS
typedef struct timing
{
  sc_statistics_t *stat;
  bfamo_device_t *bdevice;
  occaStream cmdx;
  occaStream copy;
  int sc_package_id;
} stats_t;

extern stats_t _t;
#endif

#define BFAMO_NULL_FREE(f, p)                                                  \
  do                                                                           \
  {                                                                            \
    if (p)                                                                     \
    {                                                                          \
      f(p);                                                                    \
      (p) = NULL;                                                              \
    }                                                                          \
  } while (0)

#define BFAMO_BFAM_FREE(p) BFAMO_NULL_FREE(bfam_free, p)
#define BFAMO_BFAM_FREE_ALIGNED(p) BFAMO_NULL_FREE(bfam_free_aligned, p)
#define BFAMO_OCCA_MEMORY_FREE(p) BFAMO_NULL_FREE(occaMemoryFree, p)
#define BFAMO_OCCA_KERNEL_FREE(p) BFAMO_NULL_FREE(occaKernelFree, p)
#define BFAMO_OCCA_KERNEL_INFO_FREE(p) BFAMO_NULL_FREE(occaKernelInfoFree, p)
#define BFAMO_OCCA_DEVICE_FREE(p) BFAMO_NULL_FREE(occaDeviceFree, p)
#define BFAMO_OCCA_STREAM_FREE(p) BFAMO_NULL_FREE(occaStreamFree, p)

typedef struct bfamo_kernel
{
  occaKernel kernel;
  int dim;
  occaDim items;
  occaDim groups;
  char *key;
  char *stat_key;
  char *flops_key;
  char *memread_key;

  int flops1;
  int flops2;
  int memread;
} bfamo_kernel_t;

typedef struct bfamo_sub_kernel
{
  bfamo_kernel_t *bk;
  int dim;
  occaDim items;
  occaDim groups;

  int work1;
  int work2;
} bfamo_sub_kernel_t;

bfamo_sub_kernel_t bfamo_sub_kernel_null();

void bfamo_sub_kernel_set_stats(bfamo_sub_kernel_t *bsk, int work1, int work2,
                                int flops1, int flops2, int memread);

bfamo_sub_kernel_t
bfamo_build_sub_kernel(bfam_dictionary_t *c_kernel_cache, const char *postfix,
                       bfamo_device_t *bdevice, const char *str,
                       const char *function_name, occaKernelInfo info);

void bfamo_kernel_set_working_dims(bfamo_kernel_t *bk, int dim, occaDim items,
                                   occaDim groups);

void bfamo_sub_kernel_set_working_dims(bfamo_sub_kernel_t *bsk, int dim,
                                       occaDim items, occaDim groups);

occaMemory bfamo_device_malloc(bfamo_device_t *bdevice, size_t bytecount,
                               void *src);

occaKernelInfo bfamo_kernelinfo_new(int occa_mode);

void bfamo_add_operator_metric_real(bfamo_device_t *bdevice,
                                    bfam_dictionary_t *c_ops,
                                    bfam_dictionary_t *d_ops,
                                    const char *prefix, const int N,
                                    const int Np);

void bfamo_add_operator_bfamo_real(bfamo_device_t *bdevice,
                                   bfam_dictionary_t *c_ops,
                                   bfam_dictionary_t *d_ops, const char *prefix,
                                   const int N, const int Np);

void bfamo_add_proj_operators(bfamo_device_t *bdevice, bfam_dictionary_t *c_ops,
                              bfam_dictionary_t *N2N, const int N_src,
                              const int N_dst);

#ifdef COMPUTE_STATS
void bfamo_add_statisitics(bfam_domain_t *domain,
                           bfam_communicator_t *inter_comm, occaDevice device,
                           const char **tags, stats_t *tm);
#endif

/* TODO: Move all these calls to the communicator */
void bfamo_communicator_post_send(bfam_communicator_t *comm);
void bfamo_communicator_post_recv(bfam_communicator_t *comm);
void bfamo_communicator_send_wait(bfam_communicator_t *comm);
void bfamo_communicator_recv_wait(bfam_communicator_t *comm);
/* end TODO */

int bfamo_free_c_fields(const char *key, void *val, void *arg);
int bfamo_free_c_kernel_cache(const char *key, void *val, void *arg);
void bfamo_prefs_print_file(const char *filename, int rank, int loglevel);
int bfamo_weight_fn(p4est_t *p4est, p4est_topidx_t which_tree,
                    p4est_quadrant_t *quadrant);

typedef struct
{
  int n_stages;
  bfamo_real_t *A;
  bfamo_real_t *B;
  bfamo_real_t *C;
} bfamo_lsrk_t;

typedef enum bfamo_lsrk_method {
  BFAMO_TS_LSRK_KC54,
  BFAMO_TS_LSRK_FE,
  BFAMO_TS_LSRK_HEUN,
  BFAMO_TS_LSRK_W33,
  BFAMO_TS_LSRK_NOOP,
} bfamo_lsrk_method_t;

typedef void (*bfamo_aux_rates_t)(bfam_subdomain_t *thisSubdomain,
                                  const char *prefix, void *user_data);

void bfamo_lsrk_free(bfamo_lsrk_t *ts);
void bfamo_lsrk_init(bfamo_lsrk_t *ts, bfam_domain_t *dom,
                     const char *lsrk_prefix, bfamo_lsrk_method_t method,
                     bfam_domain_match_t subdom_match, const char **subdom_tags,
                     bfamo_aux_rates_t aux_rates, void *user_data);
#endif

// }}}

#endif
