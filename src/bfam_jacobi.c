#include <bfam_base.h>
#include <bfam_jacobi.h>
#include <bfam_util.h>

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

static void
bfam_jacobi_gauss_quadrature_half(bfam_long_real_t alpha,
    bfam_long_real_t beta, int N, int half,
    bfam_long_real_t *restrict x,
    bfam_long_real_t *restrict w)
{
  BFAM_ASSERT(N>=0);
  BFAM_ASSERT(alpha>=BFAM_LONG_REAL(-1.0));
  BFAM_ASSERT(beta >=BFAM_LONG_REAL(-1.0));
  BFAM_ASSERT(!(BFAM_LONG_REAL_APPROX_EQ(alpha, BFAM_LONG_REAL(-0.5), 10) &&
                BFAM_LONG_REAL_APPROX_EQ(beta,  BFAM_LONG_REAL(-0.5), 10)));

  const int MAX_ITERATIONS = 200;

  int nk = (half) ?
    (((N + 1)%2) ? (N + 1)/2 + 1 : (N + 1)/2) /* ceil((N + 1)/2) */
    : (N + 1)/2; /* floor((N + 1)/2) */

  if (nk == 0)
    return;

  bfam_long_real_t tworho = 2*(N+1)+alpha+beta+1;
  bfam_long_real_t * tmp;

  bfam_long_real_t *restrict theta0 =
    bfam_malloc_aligned(nk*sizeof(bfam_long_real_t));
  bfam_long_real_t *restrict theta1 =
    bfam_malloc_aligned(nk*sizeof(bfam_long_real_t));
  bfam_long_real_t *restrict p0 =
    bfam_malloc_aligned(nk*sizeof(bfam_long_real_t));
  bfam_long_real_t *restrict dp0 =
    bfam_malloc_aligned(nk*sizeof(bfam_long_real_t));

  BFAM_ASSUME_ALIGNED(theta0, 32);
  BFAM_ASSUME_ALIGNED(theta1, 32);
  BFAM_ASSUME_ALIGNED(p0, 32);
  BFAM_ASSUME_ALIGNED(dp0, 32);


  /*
   * Use Gatteschi and Pittaluga's approximation for the roots of the Jacobi
   * polynomials as an initial guess.  See equation (3.19) of
   * Nicholas Hale and Alex Townsend ``Fast and Accurate Computation of
   * Gauss–Legendre and Gauss–Jacobi Quadrature Nodes and Weights'' SIAM J.
   * SCI. COMPUT. Vol. 35, No. 2, pp. A652–A674.
   */
  for (int k=nk; k > 0; --k)
  {
    int khat = (half) ? nk - k : k - 1;

    bfam_long_real_t phik =
      (2*k + alpha - BFAM_LONG_REAL(0.5))*BFAM_LONG_REAL_PI / tworho;

    theta1[khat] = phik + 1/(tworho * tworho)
        * ((BFAM_LONG_REAL(0.25) - alpha*alpha)
           * 1/BFAM_LONG_REAL_TAN(BFAM_LONG_REAL(0.5)*phik) -
           (BFAM_LONG_REAL(0.25) - beta*beta)
           * BFAM_LONG_REAL_TAN(BFAM_LONG_REAL(0.5)*phik));
  }

  /*
   * Use Newton's method for finding the roots of the Jacobi polynomial.
   */
  int converged = 0;
  for (int i=0; i < MAX_ITERATIONS; ++i)
  {
    tmp    = theta0;
    theta0 = theta1;
    theta1 = tmp;

    for (int k=0; k < nk; ++k)
    {
      x[k] = BFAM_LONG_REAL_COS(theta0[k]);
    }

         bfam_jacobi_p(alpha, beta, N+1, nk, x,  p0);
    bfam_grad_jacobi_p(alpha, beta, N+1, nk, x, dp0);

    for (int k=0; k < nk; ++k)
    {
      theta1[k] = theta0[k] - p0[k]/(-BFAM_LONG_REAL_SIN(theta0[k])*dp0[k]);
    }

    int diff = 0;
    for (int k=0; k < nk; ++k)
    {
      diff += !BFAM_LONG_REAL_APPROX_EQ(theta0[k], theta1[k], 10);
    }
    if (!diff)
    {
      converged = 1;
      break;
    }
  }

  BFAM_ABORT_IF(!converged,
      "Newton's method does not converge when computing Jacobi Gauss points");
  /*
   * Nodes
   */
  for (int k=0; k < nk; ++k)
  {
    x[k] = BFAM_LONG_REAL_COS(theta1[k]);
  }

  /*
   * Weights
   */
  bfam_grad_jacobi_p(alpha, beta, N+1, nk, x, dp0);

  for (int k=0; k < nk; ++k)
  {
    bfam_long_real_t sint = BFAM_LONG_REAL_SIN(theta1[k]);
    w[k] = tworho/(sint*sint*dp0[k]*dp0[k]);
  }

  bfam_free_aligned(theta0);
  bfam_free_aligned(theta1);
  bfam_free_aligned(p0);
  bfam_free_aligned(dp0);

  return;
}

void
bfam_jacobi_gauss_quadrature(bfam_long_real_t alpha, bfam_long_real_t beta,
    int N, bfam_long_real_t *restrict x, bfam_long_real_t *restrict w)
{
  int nk_floor = (N+1)/2; /* floor((N + 1)/2) */

  bfam_jacobi_gauss_quadrature_half(alpha, beta, N, 1, x+nk_floor, w+nk_floor);
  bfam_jacobi_gauss_quadrature_half(beta, alpha, N, 0, x, w);

  for (int k=0; k < nk_floor; ++k)
    x[k] *= -1;


  return;
}

void
bfam_jacobi_gauss_lobatto_quadrature(bfam_long_real_t alpha,
    bfam_long_real_t beta, int N, bfam_long_real_t *restrict x,
    bfam_long_real_t *restrict w)
{
  BFAM_ASSERT(N>=1);

  x[0] = -1;
  x[N] =  1;

  if (N > 1)
  {
    bfam_jacobi_gauss_quadrature(alpha + 1, beta + 1, N-2, x+1, w+1);
  }

  bfam_jacobi_p(alpha, beta, N, N+1, x, w);
  bfam_long_real_t fac = (2*N + alpha + beta + 1) / (N*(N + alpha + beta + 1));
  for (int k=0; k < N+1; ++k)
  {
    w[k] = fac/(w[k]*w[k]);
  }

  w[0] *= (1 + beta);
  w[N] *= (1 + alpha);

  return;
}

void
bfam_jacobi_p_vandermonde(bfam_long_real_t alpha, bfam_long_real_t beta, int N,
    size_t nx, bfam_long_real_t *x, bfam_long_real_t *V)
{
  for (int j=0; j <= N; ++j)
    bfam_jacobi_p(alpha, beta, j, nx, x, V + j*nx);

  return;
}

void
bfam_grad_jacobi_p_vandermonde(bfam_long_real_t alpha, bfam_long_real_t beta,
    int N, size_t nx, bfam_long_real_t *x, bfam_long_real_t *V)
{
  for (int j=0; j <= N; ++j)
    bfam_grad_jacobi_p(alpha, beta, j, nx, x, V + j*nx);

  return;
}

void
bfam_jacobi_p_interpolation(bfam_long_real_t alpha, bfam_long_real_t beta,
    int N, size_t nx, bfam_long_real_t *x, bfam_long_real_t *V,
    bfam_long_real_t *I)
{

  bfam_long_real_t *Vx = bfam_malloc_aligned((N+1)*nx*sizeof(bfam_long_real_t));

  bfam_jacobi_p_vandermonde(alpha, beta, N, nx, x, Vx);

  bfam_util_forwardslash(nx, N+1, Vx, V, I);

  bfam_free_aligned(Vx);

  return;
}

void
bfam_jacobi_p_differentiation(bfam_long_real_t alpha, bfam_long_real_t beta,
    int N, size_t nx, bfam_long_real_t *x, bfam_long_real_t *V,
    bfam_long_real_t *D)
{

  bfam_long_real_t *Vx = bfam_malloc_aligned((N+1)*nx*sizeof(bfam_long_real_t));

  bfam_grad_jacobi_p_vandermonde(alpha, beta, N, nx, x, Vx);

  bfam_util_forwardslash(nx, N+1, Vx, V, D);

  bfam_free_aligned(Vx);

  return;
}

void
bfam_jacobi_p_mass(bfam_long_real_t alpha, bfam_long_real_t beta,
    int N, bfam_long_real_t *V, bfam_long_real_t *M)
{
  bfam_long_real_t *I =
    bfam_malloc_aligned((N+1)*(N+1)*sizeof(bfam_long_real_t));

  bfam_long_real_t *invV =
    bfam_malloc_aligned((N+1)*(N+1)*sizeof(bfam_long_real_t));

  bfam_long_real_t *invVT =
    bfam_malloc_aligned((N+1)*(N+1)*sizeof(bfam_long_real_t));

  for(int i = 0; i < (N+1)*(N+1); ++i)
    I[i] = 0;

  for(int i = 0; i < (N+1)*(N+1); ++i)
    M[i] = 0;

  for(int i = 0; i <= N; ++i)
    I[(N+1)*i + i] = 1;

  bfam_util_backslash(N+1, N+1, V, I, invV);

  bfam_util_mtranspose(N+1, N+1, invV, N+1, invVT, N+1);

  bfam_util_mmmult(N+1, N+1, N+1, invVT, N+1, invV, N+1, M, N+1);

  bfam_free_aligned(I);
  bfam_free_aligned(invV);
  bfam_free_aligned(invVT);

  return;
}
