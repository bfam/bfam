#ifndef BFAM_JACOBI_H
#define BFAM_JACOBI_H


/** Compute normalized Jacobi Polynomials at \a x
 *
 * \param[in]  alpha Jacobi polynomial parameter
 * \param[in]  beta  Jacobi polynomial parameter
 * \param[in]  N     Jacobi polynomial order
 * \param[in]  nx    Number of locations to evaluate the Jacobi polynomial
 * \param[in]  x     Array of length \a nx containing the locations in [-1,1]
 *                   to evaluate the Jacobi polynomial
 * \param[out] P     An array of length \a nx containing the normalized Jacobi
 *                   polynomial $p^{(\alpha,\beta)}_N$ evaluated at \a x.
 */
void
bfam_jacobi_p(bfam_long_real_t alpha, bfam_long_real_t beta, int N,
    size_t nx, bfam_long_real_t *x, bfam_long_real_t *P);

/** Compute the derivative of the normalized Jacobi Polynomials at \a x
 *
 * \param[in]  alpha Jacobi polynomial parameter
 * \param[in]  beta  Jacobi polynomial parameter
 * \param[in]  N     Jacobi polynomial order
 * \param[in]  nx    Number of locations to evaluate the Jacobi polynomial
 * \param[in]  x     Array of length \a nx containing the locations in [-1,1]
 *                   to evaluate the Jacobi polynomial
 * \param[out] dP    An array of length \a nx containing the derivative of the
 *                   normalized Jacobi polynomial $p^{(\alpha,\beta)}_N$
 *                   evaluated at \a x.
 */
void
bfam_grad_jacobi_p(bfam_long_real_t alpha, bfam_long_real_t beta, int N,
    size_t nx, bfam_long_real_t *x, bfam_long_real_t *P);

/** Compute the nodes and weights of the Jacobi Gauss quadrature
 *
 * \param[in]  alpha Jacobi polynomial parameter
 * \param[in]  beta  Jacobi polynomial parameter
 * \param[in]  N     Jacobi polynomial order
 * \param[out] x     Nodes for the Jacobi Gauss quadrature
 * \param[out] w     Weights for the Jacobi Gauss quadrature
 */
void
bfam_jacobi_gauss_quadrature(bfam_long_real_t alpha, bfam_long_real_t beta,
    int N, bfam_long_real_t *x, bfam_long_real_t *w);

/** Compute the nodes and weights of the Jacobi Gauss Lobatto quadrature
 *
 * \param[in]  alpha Jacobi polynomial parameter
 * \param[in]  beta  Jacobi polynomial parameter
 * \param[in]  N     Jacobi polynomial order
 * \param[out] x     Nodes for the Jacobi Gauss Lobatto quadrature
 * \param[out] w     Weights for the Jacobi Gauss Lobatto quadrature
 */
void
bfam_jacobi_gauss_lobatto_quadrature(bfam_long_real_t alpha,
    bfam_long_real_t beta, int N, bfam_long_real_t *x, bfam_long_real_t *w);

/** Compute the Jacobi Vandermonde matrix.
 *
 * \param[in]  alpha Jacobi polynomial parameter
 * \param[in]  beta  Jacobi polynomial parameter
 * \param[in]  N     Jacobi polynomial order
 * \param[in]  nx    Number of locations to evaluate the Jacobi polynomial
 * \param[in]  x     Array of length \a nx containing the locations in [-1,1]
 *                   to evaluate the Jacobi polynomial
 * \param[out] V     A \a nx by \a N+1 matrix in column first order where
 *                   the $i,j$ entry contains the normalized Jacobi polynomials
 *                   $p^{(\alpha,\beta)}_j$ evaluated at \a x[i].
 */
void
bfam_jacobi_p_vandermonde(bfam_long_real_t alpha, bfam_long_real_t beta, int N,
    size_t nx, bfam_long_real_t *x, bfam_long_real_t *V);

/** Compute the Gradient Jacobi Vandermonde matrix.
 *
 * \param[in]  alpha Jacobi polynomial parameter
 * \param[in]  beta  Jacobi polynomial parameter
 * \param[in]  N     Jacobi polynomial order
 * \param[in]  nx    Number of locations to evaluate the Jacobi polynomial
 * \param[in]  x     Array of length \a nx containing the locations in [-1,1]
 *                   to evaluate the Jacobi polynomial
 * \param[out] V     A \a nx by \a N+1 matrix in column first order where
 *                   the $i,j$ entry contains the normalized Jacobi polynomials
 *                   $p^{(\alpha,\beta)}_j$ evaluated at \a x[i].
 */
void
bfam_grad_jacobi_p_vandermonde(bfam_long_real_t alpha, bfam_long_real_t beta,
    int N, size_t nx, bfam_long_real_t *x, bfam_long_real_t *V);

/** Compute the interpolation matrix using Jacobi polynomials.
 *
 * \param[in]  alpha Jacobi polynomial parameter
 * \param[in]  beta  Jacobi polynomial parameter
 * \param[in]  N     Jacobi polynomial order
 * \param[in]  nx    Number of locations to evaluate the Jacobi polynomial
 * \param[in]  x     Array of length \a nx containing the locations in [-1,1]
 *                   to evaluate the interpolating polynomial
 * \param[in]  V     A \a N+1 by \a N+1 matrix in column first order where
 *                   the $i,j$ entry contains the normalized Jacobi polynomials
 *                   $p^{(\alpha,\beta)}_j$ evaluated at the N+1 interpolating
 *                   points r[i] in [-1,1].
 * \param[out] I     Interpolation matrix.
 *
 */
void
bfam_jacobi_p_interpolation(bfam_long_real_t alpha, bfam_long_real_t beta,
    int N, size_t nx, bfam_long_real_t *x, bfam_long_real_t *V,
    bfam_long_real_t *I);

/** Compute the differentiation matrix using Jacobi polynomials.
 *
 * \param[in]  alpha Jacobi polynomial parameter
 * \param[in]  beta  Jacobi polynomial parameter
 * \param[in]  N     Jacobi polynomial order
 * \param[in]  nx    Number of locations to evaluate the Jacobi polynomial
 * \param[in]  x     Array of length \a nx containing the locations in [-1,1]
 *                   to evaluate the interpolating polynomial
 * \param[in]  V     A \a N+1 by \a N+1 matrix in column first order where
 *                   the $i,j$ entry contains the normalized Jacobi polynomials
 *                   $p^{(\alpha,\beta)}_j$ evaluated at the N+1 interpolating
 *                   points r[i] in [-1,1].
 * \param[out] D     Differentiation matrix.
 *
 */
void
bfam_jacobi_p_differentiation(bfam_long_real_t alpha, bfam_long_real_t beta,
    int N, size_t nx, bfam_long_real_t *x, bfam_long_real_t *V,
    bfam_long_real_t *D);

/** Compute the mass matrix using Jacobi polynomials.
 *
 * \param[in]  alpha Jacobi polynomial parameter
 * \param[in]  beta  Jacobi polynomial parameter
 * \param[in]  N     Jacobi polynomial order
 * \param[in]  V     A \a N+1 by \a N+1 matrix in column first order where
 *                   the $i,j$ entry contains the normalized Jacobi polynomials
 *                   $p^{(\alpha,\beta)}_j$ evaluated at the N+1 interpolating
 *                   points r[i] in [-1,1].
 * \param[out] M     mass matrix.
 *
 */
void
bfam_jacobi_p_mass(bfam_long_real_t alpha, bfam_long_real_t beta,
    int N, bfam_long_real_t *V, bfam_long_real_t *M);

#endif
