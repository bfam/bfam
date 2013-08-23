#include <bfam_base.h>
#include <bfam_util.h>

void
bfam_util_strcsl(char *str, const char **list)
{
  str[0] = '\0';

  if(list != NULL)
  {
    for(size_t l = 0; list[l]; ++l)
    {
      strcat(str, list[l]);
      if(list[l+1])
        strcat(str, ",");
    }
  }
}


void
bfam_util_mtranspose(size_t m, size_t n, bfam_long_real_t *restrict A,
                     size_t lda, bfam_long_real_t *restrict B, size_t ldb)
{
  for(size_t i = 0; i < m; ++i)
    for(size_t j = 0; j < n; ++j)
      B[i * ldb + j] = A[j * lda + i];
}

void
bfam_util_mTmmult(size_t m, size_t n, size_t k,
                 bfam_long_real_t *restrict A, size_t lda,
                 bfam_long_real_t *restrict B, size_t ldb,
                 bfam_long_real_t *restrict C, size_t ldc)
{
  for(size_t i = 0; i < m; ++i)
    for(size_t j = 0; j < n; ++j)
      for(size_t p = 0; p < k; ++p)
        C[j * ldc +i ] += A[i * lda + p] * B[j * ldb + p];
}

void
bfam_util_mmmult(size_t m, size_t n, size_t k,
                 bfam_long_real_t *restrict A, size_t lda,
                 bfam_long_real_t *restrict B, size_t ldb,
                 bfam_long_real_t *restrict C, size_t ldc)
{
  for(size_t i = 0; i < m; ++i)
    for(size_t j = 0; j < n; ++j)
      for(size_t p = 0; p < k; ++p)
        C[j * ldc +i ] += A[p * lda + i] * B[j * ldb + p];
}

void
bfam_util_mvmult(size_t m, size_t n, bfam_long_real_t *restrict A, size_t lda,
                 bfam_long_real_t *restrict b, bfam_long_real_t *restrict c)
{
  for(size_t j = 0; j < n; ++j)
    for(size_t i = 0; i < m; ++i)
      c[i] += A[j * lda + i] * b[j];
}

void
bfam_util_lu_factor(size_t n, bfam_long_real_t *restrict A,
                    size_t *restrict p, size_t *restrict q)
{
  /*
   * Algorithm 3.4.2 (Gaussian Elimination with Complete Pivoting) from p. 118
   * of Matrix Computations, Third Edition, by Golub and van Loan.
   */

  /*
   * Initialize pivots
   */
  for(size_t k = 0; k < n; ++k)
  {
    p[k] = q[k] = k;
  }

  for(size_t k = 0; k < n - 1; ++k)
  {
    size_t mu = k, lambda = k;
    bfam_long_real_t a_max = 0.0;

    for(size_t j = k; j < n; ++j)
    {
      for(size_t i = k; i < n; ++i)
      {
        const bfam_long_real_t a_abs = BFAM_LONG_REAL_ABS(A[i + n * j]);
        if(a_abs > a_max)
        {
          a_max = a_abs;
          mu = i;
          lambda = j;
        }
      }
    }

    /*
     * Swap rows
     */
    const size_t ptmp = p[k];
    p[k] = p[mu];
    p[mu] = ptmp;
    for(size_t j = 0; j < n; ++j)
    {
      const bfam_long_real_t rtmp = A[k + n * j];
      A[k + n * j] = A[mu + n * j];
      A[mu + n * j] = rtmp;
    }

    /*
     * Swap columns
     */
    const size_t qtmp = q[k];
    q[k] = q[lambda];
    q[lambda] = qtmp;
    for(size_t i = 0; i < n; ++i)
    {
      const bfam_long_real_t rtmp = A[i + n * k];
      A[i + n * k] = A[i + n * lambda];
      A[i + n * lambda] = rtmp;
    }

    if(BFAM_LONG_REAL_ABS(A[k + n * k]) >
        (BFAM_LONG_REAL_EPS * BFAM_LONG_REAL_EPS))
    {
      for(size_t i = k + 1; i < n; ++i)
      {
        A[i + n * k] = A[i + n * k] / A[k + n * k];
      }

      for(size_t i = k + 1; i < n; ++i)
      {
        for(size_t j = k + 1; j < n; ++j)
        {
          A[i + n * j] =
            A[i + n * j] - A[i + n * k] * A[k + n * j];
        }
      }
    }
  }
}

void
bfam_util_lu_solve(size_t n, bfam_long_real_t *restrict LU,
                   bfam_long_real_t *restrict x, size_t *restrict p,
                   size_t *restrict q, bfam_long_real_t *restrict work)
{
  /*
   * Compute $Pb$
   */
  for(size_t k = 0; k < n; ++k)
  {
    work[k] = x[p[k]];
  }

  /*
   * Algorithm 3.1.3 (Forward Substitution: Column Version) from p. 90 of
   * Matrix Computations, Third Edition, by Golub and van Loan.
   *
   * Note: here we have L is unit lower triangular.
   *
   * Solve $Ly=Pb$.
   */
  for(size_t j = 0; j < n - 1; ++j)
  {
    /*
     * work[j] = work[j] / LU[j + n * j];
     */
    for(size_t i = j+1; i < n; ++i)
    {
      work[i] = work[i] - work[j] * LU[i + n * j];
    }
  }
  /*
   * work[n - 1] = work[n - 1] / LU[n - 1 + n * (n - 1)];
   */

  /*
   * Algorithm 3.1.4 (Back Substitution: Column Version) from p. 90 of
   * Matrix Computations, Third Edition, by Golub and van Loan.
   *
   * Solve $Uw=y$.
   */
  for(size_t j = n - 1; j > 0; --j)
  {
    work[j] = work[j] / LU[j + n * j];
    for(size_t i = 0; i < j; ++i)
    {
      work[i] = work[i] - work[j] * LU[i + n * j];
    }
  }
  work[0] = work[0] / LU[0 + n * 0];

  /*
   * Compute $Qw$
   */
  for(size_t k = 0; k < n; ++k)
  {
    x[q[k]] = work[k];
  }
}

void
bfam_util_forwardslash(size_t m, size_t n, bfam_long_real_t *restrict A,
                       bfam_long_real_t *restrict B,
                       bfam_long_real_t *restrict C)
{
  bfam_long_real_t *AT = bfam_malloc_aligned(m*m*sizeof(bfam_long_real_t));
  bfam_long_real_t *BT = bfam_malloc_aligned(n*m*sizeof(bfam_long_real_t));
  bfam_long_real_t *CT = bfam_malloc_aligned(n*m*sizeof(bfam_long_real_t));

  bfam_util_mtranspose(m, n, A, m, AT, n);
  bfam_util_mtranspose(n, n, B, n, BT, n);

  bfam_util_backslash(n, m, BT, AT, CT);

  bfam_util_mtranspose(n, m, CT, n, C, m);

  bfam_free_aligned(CT);
  bfam_free_aligned(BT);
  bfam_free_aligned(AT);
}

void
bfam_util_backslash(size_t m, size_t n, bfam_long_real_t *restrict A,
                    bfam_long_real_t *restrict B,
                    bfam_long_real_t *restrict C)
{
  bfam_long_real_t *LU   = bfam_malloc_aligned(m*m*sizeof(bfam_long_real_t));
  bfam_long_real_t *work = bfam_malloc_aligned(  m*sizeof(bfam_long_real_t));

  size_t *p = bfam_malloc_aligned(m*sizeof(size_t));
  size_t *q = bfam_malloc_aligned(m*sizeof(size_t));

  memcpy(LU, A, m*m*sizeof(bfam_long_real_t));
  memcpy( C, B, m*n*sizeof(bfam_long_real_t));

  bfam_util_lu_factor(m, LU, p, q);

  for(size_t j = 0; j < n; ++j)
    bfam_util_lu_solve(m, LU, C+j*m, p, q, work);

  bfam_free_aligned(q);
  bfam_free_aligned(p);
  bfam_free_aligned(work);
  bfam_free_aligned(LU);
}

void bfam_util_transfinite(
    bfam_real_t *x, bfam_real_t *y, bfam_real_t *z,
    const bfam_gloidx_t *N, const bfam_locidx_t *Nl, const bfam_gloidx_t *gx,
    const bfam_long_real_t *xc, const bfam_long_real_t *xe,
    const bfam_long_real_t *r,
    const bfam_long_real_t *yc, const bfam_long_real_t *ye,
    const bfam_long_real_t *s,
    const bfam_long_real_t *zc, const bfam_long_real_t *ze,
    const bfam_long_real_t *t)
{
}

/** Linear blending
 */
void bfam_util_linear_blend(
    bfam_real_t *restrict x, bfam_real_t *restrict y, bfam_real_t *restrict z,
    const int dim,
    const bfam_gloidx_t *N, const bfam_locidx_t *Nltmp,
    const bfam_gloidx_t *gxtmp, const bfam_long_real_t *xc,
    const bfam_long_real_t *yc, const bfam_long_real_t *zc)
{

  bfam_locidx_t Nl[3] = {0,0,0};
  if(Nltmp == NULL) for(int d = 0;d < dim; d++) Nl[d] = (bfam_locidx_t)N[d];
  else for(int d = 0;d < dim; d++) Nl[d] = Nltmp[d];

  bfam_gloidx_t gx[3] = {0,0,0};
  if(gxtmp != NULL) for(int d = 0;d < dim; d++) gx[d] = gxtmp[d];

  for(bfam_locidx_t k = 0; k <= Nl[2]; k++)
    for(bfam_locidx_t j = 0; j <= Nl[1]; j++)
      for(bfam_locidx_t i = 0; i <= Nl[0]; i++)
      {
        bfam_locidx_t ix = i + (Nl[0]+1)*(j + (Nl[1]+1)*k);

        bfam_long_real_t a = ((bfam_long_real_t) (gx[0]+i))/N[0];
        if(N[0] == 0) a = 0;

        if(x != NULL) x[ix] = (1-a)*xc[0]+a*xc[1];
        if(y != NULL) y[ix] = (1-a)*yc[0]+a*yc[1];
        if(z != NULL) z[ix] = (1-a)*zc[0]+a*zc[1];
        if(dim == 1) continue;

        bfam_long_real_t b = ((bfam_long_real_t) (gx[1]+j))/N[1];
        if(N[1] == 0) b = 0;

        if(x != NULL) x[ix] *= (1-b);
        if(x != NULL) x[ix] += b*((1-a)*xc[2]+a*xc[3]);
        if(y != NULL) y[ix] *= (1-b);
        if(y != NULL) y[ix] += b*((1-a)*yc[2]+a*yc[3]);
        if(z != NULL) z[ix] *= (1-b);
        if(z != NULL) z[ix] += b*((1-a)*zc[2]+a*zc[3]);

        if(dim == 2) continue;

        bfam_long_real_t c = ((bfam_long_real_t) (gx[2]+k))/N[2];
        if(N[2] == 0) c = 0;
        if(x != NULL) x[ix] *= (1-c);
        if(x != NULL) x[ix] += c*((1-b)*((1-a)*xc[4]+a*xc[5]) + b*((1-a)*xc[6]+a*xc[7]));
        if(y != NULL) y[ix] *= (1-c);
        if(y != NULL) y[ix] += c*((1-b)*((1-a)*yc[4]+a*yc[5]) + b*((1-a)*yc[6]+a*yc[7]));
        if(z != NULL) z[ix] *= (1-c);
        if(z != NULL) z[ix] += c*((1-b)*((1-a)*zc[4]+a*zc[5]) + b*((1-a)*zc[6]+a*zc[7]));
      }
}
