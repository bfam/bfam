#ifndef BFAM_LOG_H
#define BFAM_LOG_H

#include <bfam_base.h>

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

#endif
