#ifndef BFAM_DICTIONARY_H
#define BFAM_DICTIONARY_H

#define BFAM_PTR_STR_LEN BFAM_BUFSIZ

#include <bfam_base.h>
#include <bfam_critbit.h>

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
int bfam_dictionary_delete(bfam_dictionary_t *d, const char *key);

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
char *bfam_dictionary_get_value(bfam_dictionary_t *d, const char *key);

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
int bfam_dictionary_get_value_int(bfam_dictionary_t *d, const char *key,
                                  int *val);

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

#endif
