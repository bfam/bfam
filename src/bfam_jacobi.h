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

#endif
