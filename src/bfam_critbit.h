#ifndef BFAM_CRITBIT_H
#define BFAM_CRITBIT_H

typedef struct {
  void *root;
} bfam_critbit0_tree_t;

/** Membership testing.
 *
 * The following function takes a tree, \a t, and a \c NULL terminated string,
 * \a u,  and returns non-zero iff \a u in \a t.
 *
 * \param [in] t tree
 * \param [in] u possible member
 * \returns non-zero iff \a u in \a t
 */
int bfam_critbit0_contains(bfam_critbit0_tree_t *t, const char *u);

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
int bfam_critbit0_insert(bfam_critbit0_tree_t *t, const char *u);

/** Deleting elements.
 *
 * This function takes a tree, \a t, and a \c NULL terminated string,
 * \a u, and possibly mutates the tree such that $u \notin t$.
 *
 * \param [in,out] t tree
 * \param [in] u possible member to remove from \a t
 * \returns It returns 1 if the tree was mutated, 0 otherwise.
 */
int bfam_critbit0_delete(bfam_critbit0_tree_t *t, const char *u);

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
void bfam_critbit0_clear(bfam_critbit0_tree_t *t);

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
int bfam_critbit0_allprefixed(bfam_critbit0_tree_t *t, const char *prefix,
                              int (*handle) (const char *, void *), void *arg);

#endif
