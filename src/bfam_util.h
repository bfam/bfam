#ifndef BFAM_UTIL_H
#define BFAM_UTIL_H

/** Creates a comma separated list of strings.
 *
 * \param [out] str  comma separated list of strings
 * \param [in]  list list of strings to concatenate
 *
 * \note \a str must have enough spaces for the result.
 */
void bfam_util_strcsl(char *str, const char **list);

#endif
