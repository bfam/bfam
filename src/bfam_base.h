#ifndef BFAM_BASE_H
#define BFAM_BASE_H

#if !defined(__APPLE__) && !defined(__FreeBSD__) && !defined(__USLC__) && \
    !defined(_M_UNIX) && !defined(__sgi) && !defined(__DragonFly__) && \
    !defined(__TANDEM)
#define _XOPEN_SOURCE 600
#endif

#include <bfam_config.h>

#include <stdarg.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <signal.h>

#include <sys/types.h>
#include <execinfo.h>
#include <errno.h>
#include <err.h>

#ifdef BFAM_HAVE_TIME_H
#include <time.h>
#endif

#ifdef BFAM_HAVE_SYS_TIME_H
#include <sys/time.h>
#endif

#ifdef BFAM_HAVE_MACH_MACH_TIME_H
#include <mach/mach_time.h>
#endif

#if defined (BFAM_HAVE_SYSEXITS_H)
#include <sysexits.h>
#elif defined (BFAM_HAVE_SYS_SYSEXITS_H)
#include <sys/sysexits.h>
#else
#define EX_OSERR EXIT_FAILURE
#define EX_USAGE EXIT_FAILURE
#endif

#if defined(__APPLE__)
#include <sys/sysctl.h>
#elif defined(_WIN32)
#include <windows.h>
#endif

#ifdef BFAM_HAVE_OPENCL
#  ifdef __APPLE__
#    include <OpenCL/cl.h>
#  else
#    include <CL/cl.h>
#  endif
#endif

#include <mpi.h>

#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wcast-align"
#include <p4est.h>
#include <p4est_bits.h>
#include <p4est_extended.h>
#include <p4est_ghost.h>
#include <p4est_lnodes.h>
#include <p4est_vtk.h>
#pragma clang diagnostic pop

#if defined __GNUC__ && !defined __GNUC_PREREQ
# ifndef __GNUC_MINOR__
#    define __GNUC_PREREQ(maj, min) 0
# else
#    define __GNUC_PREREQ(maj, min) \
         ((__GNUC__ << 16) + __GNUC_MINOR__ >= ((maj) << 16) + (min))
# endif
#endif

#define BFAM_BUFSIZ 8192

#define BFAM_APPEND(x, y) x ## y

#define BFAM_MIN(a,b) (((a)<(b))?(a):(b))
#define BFAM_MAX(a,b) (((a)>(b))?(a):(b))

#define BFAM_NOOP() do {} while(0)
#define BFAM_ABORT(...) bfam_abort_verbose(__FILE__, __LINE__, __VA_ARGS__)
#define BFAM_ABORT_IF(q,...) ((q) ? BFAM_ABORT(__VA_ARGS__) : (void) 0)
#define BFAM_ABORT_IF_NOT(q,...)  BFAM_ABORT_IF(!(q), __VA_ARGS__)

#ifdef BFAM_DEBUG
#define BFAM_ASSERT(expression)   \
  BFAM_ABORT_IF_NOT((expression), "Assert Failed: '" #expression "'")
#else
#define BFAM_ASSERT(expression) BFAM_NOOP()
#endif

#define BFAM_SYS_ERROR_CHECK(cond, msg) \
  do                                    \
  {                                     \
    if(cond)                            \
    {                                   \
      perror(msg);                      \
      BFAM_ABORT("perror");             \
    }                                   \
  } while (0)

#define BFAM_MPI_CHECK(c) BFAM_ABORT_IF_NOT((c) == MPI_SUCCESS, "MPI Error")

#define BFAM_IS_ALIGNED(p,a) (((intptr_t)(p) & ((a) - 1)) == 0)

/*
 * Switch for different compilers taken from this web page:
 *
 *     http://nadeausoftware.com/articles/2012/10/c_c_tip_how_detect_compiler_name_and_version_using_compiler_predefined_macros
 *
 */
#if defined(__clang__)
/* Clang/LLVM. ---------------------------------------------- */
#    define BFAM_ASSUME_ALIGNED(lvalueptr, align)  BFAM_NOOP()
#elif defined(__ICC) || defined(__INTEL_COMPILER)
/* Intel ICC/ICPC. ------------------------------------------ */
#    define BFAM_ASSUME_ALIGNED(lvalueptr, align) \
         __assume_aligned(lvalueptr, align)

#elif defined(__GNUC__) || defined(__GNUG__)
/* GNU GCC/G++. --------------------------------------------- */
#  if __GNUC_PREREQ(4,7)
//      If  gcc_version >= 4.7
#    define BFAM_ASSUME_ALIGNED(lvalueptr, align) \
         lvalueptr = __builtin_assume_aligned (lvalueptr, align)
#  else
//       Else
#    define BFAM_ASSUME_ALIGNED(lvalueptr, align)  BFAM_NOOP()
#  endif

#elif defined(__HP_cc) || defined(__HP_aCC)
/* Hewlett-Packard C/aC++. ---------------------------------- */
#    define BFAM_ASSUME_ALIGNED(lvalueptr, align)  BFAM_NOOP()

#elif defined(__IBMC__) || defined(__IBMCPP__)
/* IBM XL C/C++. -------------------------------------------- */
#    define BFAM_ASSUME_ALIGNED(lvalueptr, align)  BFAM_NOOP()

#elif defined(_MSC_VER)
/* Microsoft Visual Studio. --------------------------------- */
#    define BFAM_ASSUME_ALIGNED(lvalueptr, align)  BFAM_NOOP()

#elif defined(__PGI)
/* Portland Group PGCC/PGCPP. ------------------------------- */
#    define BFAM_ASSUME_ALIGNED(lvalueptr, align)  BFAM_NOOP()

#elif defined(__SUNPRO_C) || defined(__SUNPRO_CC)
/* Oracle Solaris Studio. ----------------------------------- */
#    define BFAM_ASSUME_ALIGNED(lvalueptr, align)  BFAM_NOOP()

#endif

#if (defined __GNUC__) || (defined __PGI) || (defined __IBMC__)
  #define BFAM_ALIGN(n) __attribute__((aligned(n)))
#elif (defined _MSC_VER)
  #define BFAM_ALIGN(n) __declspec(align(n))
#else
  #error Need equilvent of __attribute__((aligned(n))) for this compiler
#endif

#if defined __GNUC__
#define BFAM_ASM_COMMENT(X)  __asm__("# " X)
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

#define BFAM_APPROX_EQ(x, y, K, abs, eps, min) \
  ((abs)((x)-(y)) < (K) * (eps) * (abs)((x)+(y)) || (abs)((x)-(y)) < (min))

/* Type for setup computations */
typedef long double bfam_long_real_t;
#define BFAM_LONG_REAL(x) BFAM_APPEND(x, L)
#define BFAM_LONG_REAL_PRIe "Le"
#define BFAM_LONG_REAL_PRIf "Lf"
#define BFAM_LONG_REAL_PRIg "Lg"

#define BFAM_LONG_REAL_ABS      fabsl
#define BFAM_LONG_REAL_ATAN     atanl
#define BFAM_LONG_REAL_COS       cosl
#define BFAM_LONG_REAL_LGAMMA lgammal
#define BFAM_LONG_REAL_LOG       logl
#define BFAM_LONG_REAL_EXP       expl
#define BFAM_LONG_REAL_SQRT     sqrtl
#define BFAM_LONG_REAL_SIN       sinl
#define BFAM_LONG_REAL_TAN       tanl

#define BFAM_LONG_REAL_EPS LDBL_EPSILON
#define BFAM_LONG_REAL_MIN LDBL_MIN
#define BFAM_LONG_REAL_PI  (4*BFAM_LONG_REAL_ATAN(1))

#define BFAM_LONG_REAL_APPROX_EQ(x, y, K)                               \
  BFAM_APPROX_EQ((x), (y), (K), BFAM_LONG_REAL_ABS, BFAM_LONG_REAL_EPS, \
                 BFAM_LONG_REAL_MIN)

/* Type for runtime computations */
typedef double bfam_real_t;
#define BFAM_REAL(x) BFAM_APPEND(x, )
#define BFAM_REAL_PRIe "e"
#define BFAM_REAL_PRIf "f"
#define BFAM_REAL_PRIg "g"

#define BFAM_REAL_EPS DBL_EPSILON
#define BFAM_REAL_MIN DBL_MIN

#define BFAM_REAL_APPROX_EQ(x, y, K)                                         \
  BFAM_APPROX_EQ((x), (y), (K), BFAM_REAL_ABS, BFAM_REAL_EPS, BFAM_REAL_MIN)

/* Type for processor-local indexing */
typedef int32_t bfam_locidx_t;
#define BFAM_LOCIDX_MPI MPI_INT

/* Type for globally unique indexing */
typedef int64_t bfam_gloidx_t;
#define BFAM_GLOIDX_MPI MPI_LONG_LONG_INT

/** Abort function.
 *
 * This call will abort the program.
 *
 * This is a good function to set a breakpoint on when debugging.
 *
 * \return This function does not return.
 */
void bfam_abort() BFAM_NORETURN;


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
void * bfam_malloc(size_t size);

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
void * bfam_calloc(size_t nmemb, size_t size);

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
void * bfam_realloc(void *ptr, size_t size);

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

#endif
