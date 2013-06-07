#ifndef BFAM_DICTIONARY_H
#define BFAM_DICTIONARY_H

#include <bfam_critbit.h>

typedef struct {
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
int bfam_dictionary_insert(bfam_dictionary_t *d, const char *key,
                                  const char *val);

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
char* bfam_dictionary_get_value(bfam_dictionary_t *d, const char *key);

/** Clearing a dictionarty.
 *
 * Clearing a dictionarty (freeing all members) brings us our first code for
 * walking the whole dictionary rather than just tracing a path through it.
 *
 * So, the \c bfam_dictionary_clear function takes a dictionary, \a d, and frees
 * every member of it, mutating the dictionarty such that it is empty on exit.
 *
 * \param [in,out] d dictionary
 */
void bfam_dictionary_clear(bfam_dictionary_t *d);

#endif
