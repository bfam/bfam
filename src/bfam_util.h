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

/** In place matrix transposition.
 *
 * $A = A^T$
 *
 * \param [in]  m   The number of rows    of the matrix \a A.
 * \param [in]  n   The number of columns of the matrix \a A.
 * \param [in]  A   The matrix $A$ in column first storage.
 * \param [in]  lda The first dimension of \a A.
 * \param [out] B   The matrix $A^T$ in column first storage.
 * \param [in]  ldb The first dimension of \a B.
 *
 * \note This code is not designed for speed.
 */
void bfam_util_mtranspose(size_t m, size_t n,
                          bfam_long_real_t *restrict A, size_t lda,
                          bfam_long_real_t *restrict B, size_t ldb);

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

/** Forward slash.
 *
 * $C = A * B^{-1}$
 *
 * \param [in]  m   The number of rows of \a A and \a C.
 * \param [in]  n   The number of rows and columns of the matrix \a B and
 *                  the number of columns of \a A.
 * \param [in]  A   The matrix $A$ in column first storage.
 * \param [in]  B   The matrix $B$ in column first storage.
 * \param [out] C   The matrix $C$ in column first storage.
 *
 * \note This code is not designed for speed.
 */
void bfam_util_forwardslash(size_t m, size_t n, bfam_long_real_t *restrict A,
                            bfam_long_real_t *restrict B,
                            bfam_long_real_t *restrict C);

/** Backslash.
 *
 * $C = A^{-1} * B$
 *
 * \param [in]  m   The number of rows and columns of the matrix \a A and
 *                  the number of rows of \a C and \a B.
 * \param [in]  n   The number of columns of the matrix \a B and \a C.
 * \param [in]  A   The matrix $A$ in column first storage.
 * \param [in]  B   The matrix $B$ in column first storage.
 * \param [out] C   The matrix $C$ in column first storage.
 *
 * \note This code is not designed for speed.
 */
void bfam_util_backslash(size_t m, size_t n, bfam_long_real_t *restrict A,
                         bfam_long_real_t *restrict B,
                         bfam_long_real_t *restrict C);

/** Transfinite Interpolation
 *
 * \param [out] x   pointer to x dimension to initialize
 * \param [out] y   pointer to y dimension to initialize
 * \param [out] z   pointer to z dimension to initialize
 * \param [in]  N   global top index
 * \param [in]  Nl  local  top index
 * \param [in]  gx  my global starting index
 * \param [in]  xc  pointer to x values on corners
 * \param [in]  yc  pointer to y values on corners
 * \param [in]  zc  pointer to z values on corners
 * \param [in]  x0  pointer to x values on side 0
 * \param [in]  y0  pointer to y values on side 0
 * \param [in]  z0  pointer to z values on side 0
 * \param [in]  x1  pointer to x values on side 0
 * \param [in]  y1  pointer to y values on side 0
 * \param [in]  z1  pointer to z values on side 0
 * \param [in]  x2  pointer to x values on side 0
 * \param [in]  y2  pointer to y values on side 0
 * \param [in]  z2  pointer to z values on side 0
 * \param [in]  x3  pointer to x values on side 0
 * \param [in]  y3  pointer to y values on side 0
 * \param [in]  z3  pointer to z values on side 0
 *
 * \note If an out pointer is null, then this dimension is skipped
 *       If a side pointer is null, then linear interpolation between corners is
 *       used
 */
void bfam_util_transfinite(
          bfam_real_t   *x ,       bfam_real_t   *y ,       bfam_real_t   *z ,
    const bfam_gloidx_t *N , const bfam_locidx_t *Nl, const bfam_locidx_t *gx,
    const bfam_real_t   *xc, const bfam_real_t   *yc, const bfam_real_t   *zc,
    const bfam_real_t   *x0, const bfam_real_t   *y0, const bfam_real_t   *z0,
    const bfam_real_t   *x1, const bfam_real_t   *y1, const bfam_real_t   *z1,
    const bfam_real_t   *x2, const bfam_real_t   *y2, const bfam_real_t   *z2,
    const bfam_real_t   *x3, const bfam_real_t   *y3, const bfam_real_t   *z3);

#endif
