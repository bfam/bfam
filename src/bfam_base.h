#ifndef BFAM_BASE_H
#define BFAM_BASE_H

#define _ISOC99_SOURCE

#include <bfam_config.h>

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#if defined (BFAM_HAVE_SYSEXITS_H)
#include <sysexits.h>
#elif defined (BFAM_HAVE_SYS_SYSEXITS_H)
#include <sys/sysexits.h>
#else
#define EX_OSERR EXIT_FAILURE
#define EX_USAGE EXIT_FAILURE
#endif

#include <mpi.h>

#define BFAM_NOOP() do {} while(0)
#define BFAM_ABORT(s) bfam_abort_verbose(__FILE__, __LINE__, (s))
#define BFAM_ABORT_IF(q,s) ((q) ? BFAM_ABORT(s) : (void) 0)
#define BFAM_ABORT_IF_NOT(q,s)  BFAM_ABORT_IF(!(q),(s))

#ifdef BFAM_DEBUG
#define BFAM_ASSERT(expression)   \
  BFAM_ABORT_IF_NOT((expression), "Assert Failed: '" #expression "'")
#else
#define BFAM_ASSERT(expression) BFAM_NOOP()
#endif

/** Abort function.
 *
 * This call will abort the program.
 *
 * This is a good function to set a breakpoint on when debugging.
 *
 * \return This function does not return.
 */
void bfam_abort();


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
void bfam_abort_verbose(const char *file, int line, const char *note);

#endif
