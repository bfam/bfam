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

/** Matrix matrix multiplication.
 *
 * $C = A * B + C$
 *
 * \param [in]     m   The number of rows of the matrix \a A and \a C.
 * \param [in]     n   The number of columns of the matrix \a B and \a C.
 * \param [in]     k   The number of columns of the matrix \a A and
 *                     the number of rows    of the matrix \a B.
 * \param [in]     A   The matrix $A$ in column first storage.
 * \param [in]     lda The first dimension of \a A.
 * \param [in]     B   The matrix $B$ in column first storage.
 * \param [in]     ldb The first dimension of \a B.
 * \param [in,out] C   The matrix $C$ in column first storage.
 * \param [in]     ldc The first dimension of \a C.
 *
 * \note This code is not designed for speed.
 */
void bfam_util_mmmult(size_t m, size_t n, size_t k,
                      bfam_long_real_t *restrict A, size_t lda,
                      bfam_long_real_t *restrict B, size_t ldb,
                      bfam_long_real_t *restrict C, size_t ldc);

/** Matrix vector multiplication.
 *
 * $c = A * b + c$
 *
 * \param [in]     m   The number of rows of the matrix \a A and \a c.
 * \param [in]     n   The number of columns of the matrix \a A and
 *                     the number of rows    of the vector \a b.
 * \param [in]     A   The matrix $A$ in column first storage.
 * \param [in]     lda The first dimension of \a A.
 * \param [in]     b   The vector $b$.
 * \param [in,out] c   The matrix $c$.
 *
 * \note This code is not designed for speed.
 */
void bfam_util_mvmult(size_t m, size_t n,
                      bfam_long_real_t *restrict A, size_t lda,
                      bfam_long_real_t *restrict b,
                      bfam_long_real_t *restrict c);

/** LU Factorization using Gaussian Elimination with Complete Pivoting.
 *
 * This factors \a A in place returning pivot vector \a p and \a q.
 *
 * \param [in]      n The number of rows and columns of \a A.
 * \param [in,out]  A The matrix to be factorized in column first storage, the
 *                    matrix is overwritten with the factors upon return.
 * \param [out]     p Row    pivot vector of length \a n.
 * \param [out]     q Column pivot vector of length \a n.
 *
 * \note This code is not designed for speed.
 */
void bfam_util_lu_factor(size_t n, bfam_long_real_t *restrict A,
                         size_t *restrict p, size_t *restrict q);

/** LU Solve.
 *
 * This performs forward and backward substitution to solve a system
 * of equations given an LU factorization of a matrix.
 *
 * \param [in]     n    The number of rows and columns of \a LU.
 * \param [in]     LU   The LU factorization of a matrix.
 * \param [in,out] x    Upon entry it is the RHS of the system to solve
 *                      and upon exit it is the solution.
 * \param [in]     p    Row    pivot vector of length \a n.
 * \param [in]     q    Column pivot vector of length \a n.
 * \param [out]    work workspace.
 *
 * \note This code is not designed for speed.
 */
void bfam_util_lu_solve(size_t n, bfam_long_real_t *restrict LU,
                        bfam_long_real_t *restrict x, size_t *restrict p,
                        size_t *restrict q, bfam_long_real_t *restrict work);

#endif
