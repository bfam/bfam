#include <bfam_base.h>
#include <bfam_jacobi.h>

/*
 * This function computes the normalization of the Jacobi polynomial
 * $\left\{h_N^{(\alpha,\beta)}\right\}^{-\frac12}$ where [see @Szego39
 * (4.3.3)]
 *
 * $$
 *  \begin{aligned}
 *    h_N^{(\alpha,\beta)} &=
 *      \int_{-1}^{+1} (1-x)^{\alpha} (1+x)^{\beta}
 *    \left\{P_N^{(\alpha,\beta)} (x)\right\}^2 \, dx \\
 *    &=
 *    \frac{2^{\alpha+\beta+1}}{2N+\alpha+\beta+1}
 *    \frac{\Gamma(N+\alpha+1)\Gamma(N+\beta+1)}
 *         {\Gamma(N+\alpha+\beta+1)\Gamma(N+1)}.
 *  \end{aligned}
 * $$
 */
static bfam_long_real_t
bfam_jacobi_h_inv_sqrt(bfam_long_real_t alpha, bfam_long_real_t beta, int N)
{
  BFAM_ASSERT(N>=0);
  BFAM_ASSERT(alpha>=BFAM_LONG_REAL(-1.0));
  BFAM_ASSERT(beta >=BFAM_LONG_REAL(-1.0));
  BFAM_ASSERT(!(BFAM_LONG_REAL_APPROX_EQ(alpha, BFAM_LONG_REAL(-0.5), 10) &&
                BFAM_LONG_REAL_APPROX_EQ(beta,  BFAM_LONG_REAL(-0.5), 10)));

  bfam_long_real_t lgn = -(alpha + beta + 1)*BFAM_LONG_REAL_LOG(2)
                         - BFAM_LONG_REAL_LGAMMA(N + alpha + 1)
                         - BFAM_LONG_REAL_LGAMMA(N + beta + 1)
                         + BFAM_LONG_REAL_LOG(2*N + alpha + beta + 1)
                         + BFAM_LONG_REAL_LGAMMA(N + 1)
                         + BFAM_LONG_REAL_LGAMMA(N + alpha + beta + 1);
  return BFAM_LONG_REAL_SQRT(BFAM_LONG_REAL_EXP(lgn));
}

/*
 * This function evaluates the orthonormal polynomial $p_N(x)$ associated with
 * the Jacobi polynomial where $p_N(x) =
 * \left\{h_N^{(\alpha,\beta)}\right\}^{-\frac12} P_N^{(\alpha,\beta)} (x)$.
 *
 * The Jacobi polynomials are a set of polynomial functions for $\alpha > -1$,
 * $\beta > -1$, and $N$ a non-negative integer.  The functions are defined
 * on $[-1, +1]$ and orthogonal with respect to the weight function
 * $$w(x)=(1-x)^\alpha(1+x)^\beta.$$  Here we use the same normalization as
 * @Szego39, i.e., $P_N^{(\alpha,\beta)}(1) = {n+\alpha\choose n}$.  Thus we
 * have
 * $$
 *   \int_{-1}^{+1} p_N(x) p_M(x) w(x) \, dx = \delta_{NM}.
 * $$
 *
 * The three term recurrence relation arrived at by rearranging @Szego39
 * [(4.5.1)] is
 * $$
 *   P_n^{(\alpha,\beta)}(x) = (ax-b) P_{n-1}^{(\alpha,\beta)}(x)
 *                             - c P_{n-2}^{(\alpha,\beta)}
 * $$
 * where
 * $$
 * \begin{aligned}
 *   a &= \frac{(2n + \alpha + \beta -1)(2n + \alpha + \beta)}
 *             {2n (n + \alpha + \beta)} \\
 *   b &= \frac{(\beta^2 - \alpha^2)(2n + \alpha + \beta - 1)}
 *             {2n(n + \alpha + \beta)(2n + \alpha + \beta - 2)} \\
 *   c &= \frac{(n + \alpha - 1)(n + \beta - 1)(2n + \alpha + \beta)}
 *             {n(n + \alpha + \beta)(2n + \alpha + \beta - 2)}
 * \end{aligned}
 * $$
 * with $P_0^{(\alpha,\beta)}(x) = 1$ and
 * $P_1^{(\alpha,\beta)}(x) =  \frac12(\alpha + \beta + 2)x
 *                           + \frac12(\alpha - \beta)$.
 */
void
bfam_jacobi_p(bfam_long_real_t alpha, bfam_long_real_t beta, int N,
   size_t nx, bfam_long_real_t *x, bfam_long_real_t *P)
{
  BFAM_ASSERT(N>=0);
  BFAM_ASSERT(alpha>=BFAM_LONG_REAL(-1.0));
  BFAM_ASSERT(beta >=BFAM_LONG_REAL(-1.0));
  BFAM_ASSERT(!(BFAM_LONG_REAL_APPROX_EQ(alpha, BFAM_LONG_REAL(-0.5), 10) &&
                BFAM_LONG_REAL_APPROX_EQ(beta,  BFAM_LONG_REAL(-0.5), 10)));

  for (size_t i=0; i < nx; ++i)
  {
    bfam_long_real_t P_n_2;
    bfam_long_real_t P_n_1 = 1;
    bfam_long_real_t P_n_0 = ((alpha + beta + 2)/2)*x[i] + (alpha - beta)/2;
    if (N==0)
    {
      P[i] = P_n_1;
    }
    else if (N==1)
    {
      P[i] = P_n_0;
    }
    else
    {
      for (int n=2; n < N+1; ++n)
      {
        bfam_long_real_t a = (2*n + alpha + beta - 1)*(2*n + alpha + beta)/
          (2*n*(n + alpha + beta));
        bfam_long_real_t b = (beta*beta - alpha*alpha)*
          (2*n + alpha + beta - 1)/
          (2*n*(n + alpha + beta)*(2*n + alpha + beta - 2));
        bfam_long_real_t c = (n + alpha - 1)*(n + beta - 1)*
          (2*n + alpha + beta)/(n*(n + alpha + beta)*(2*n + alpha + beta - 2));

        P_n_2 = P_n_1;
        P_n_1 = P_n_0;
        P_n_0 = (a*x[i] - b)*P_n_1 - c*P_n_2;
      }
      P[i] = P_n_0;
    }
  }


  /*
   * Normalize the Jacobi polynomials
   */
  bfam_long_real_t h_inv_sqrt = bfam_jacobi_h_inv_sqrt(alpha, beta, N);
  for(size_t i=0; i < nx; ++i)
    P[i] *= h_inv_sqrt;

  return;
}

/*
 * This function evaluates the derivative of the orthonormal polynomial
 * $p_N(x)$ associated with the Jacobi polynomial where
 * $p_N(x) = \left\{h_N^{(\alpha,\beta)}\right\}^{-\frac12}
 * P_N^{(\alpha,\beta)} (x)$.
 *
 * For the evaluation of the derivative we use the identity
 * $$
 *   \frac{d}{dx} P_N^{(\alpha,\beta)} (x) =
 *    \frac{N+\alpha+\beta+1}{2} P_{N-1}^{(\alpha+1,\beta+1)} (x)
 * $$
 * along with
 * $$
 *   h_N^{(\alpha,\beta)} =
 *     \frac{N+\alpha+\beta+1}{4N} h_{N-1}^{(\alpha+1,\beta+1)}
 * $$
 * to get
 * $$
 * \begin{aligned}
 *   \frac{d}{dx} p_N^{(\alpha,\beta)} (x)
 *     &= \left\{h_N^{(\alpha,\beta)}\right\}^{-\frac12}
 *        \frac{d}{dx} P_N^{(\alpha,\beta)} (x) \\
 *     &= \left\{h_N^{(\alpha,\beta)}\right\}^{-\frac12}
 *        \frac{N+\alpha+\beta+1}{2}
 *        P_{N-1}^{(\alpha+1,\beta+1)} (x) \\
 *     &= \left(\frac{4N}{N+\alpha+\beta+1}\right)^{\frac12}
 *        \frac{N+\alpha+\beta+1}{2}
 *        \left\{h_{N-1}^{(\alpha+1,\beta+1)}\right\}^{-\frac12}
 *        P_{N-1}^{(\alpha+1,\beta+1)} (x) \\
 *     &= \left(N(N+\alpha+\beta+1)\right)^{\frac12}
 *        p_{N-1}^{(\alpha+1,\beta+1)} (x).
 * \end{aligned}
 * $$
 */
void
bfam_grad_jacobi_p(bfam_long_real_t alpha, bfam_long_real_t beta, int N,
    size_t nx, bfam_long_real_t *x, bfam_long_real_t *dP)
{
  BFAM_ASSERT(N>=0);
  BFAM_ASSERT(alpha>=BFAM_LONG_REAL(-1.0));
  BFAM_ASSERT(beta >=BFAM_LONG_REAL(-1.0));

  if(N==0)
  {
    for (size_t i=0; i < nx; ++i)
    {
      dP[i] = BFAM_LONG_REAL(0.0);
    }
  }
  else
  {
    bfam_jacobi_p(alpha+1, beta+1, N-1, nx, x, dP);
    bfam_long_real_t scale = BFAM_LONG_REAL_SQRT(N*(N+alpha+beta+1));
    for (size_t i=0; i < nx; ++i)
    {
      dP[i] *= scale;
    }
  }

  return;
}
