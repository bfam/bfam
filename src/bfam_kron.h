#ifndef BFAM_KRON_H
#define BFAM_KRON_H

/** $y = (I \otimes A) x$.
 *
 * \param [in]  N number of rows and columns of $A$
 * \param [in]  A column first representation of $A$
 * \param [in]  x vector $x$
 * \param [out] y vector $y$
 */
#define BFAM_KRON_IXA(N, A, x, y)                                  \
  do                                                               \
  {                                                                \
    for(int bfam_kron_n = 0; bfam_kron_n < (N)*(N); ++bfam_kron_n) \
      (y)[bfam_kron_n] = 0;                                        \
                                                                   \
    for(int bfam_kron_k = 0; bfam_kron_k < (N); ++bfam_kron_k)     \
      for(int bfam_kron_j = 0; bfam_kron_j < (N); ++bfam_kron_j)   \
        for(int bfam_kron_i = 0; bfam_kron_i < (N); ++bfam_kron_i) \
          (y)[(N) * bfam_kron_k + bfam_kron_j] +=                  \
            (A)[(N) * bfam_kron_i + bfam_kron_j] *                 \
            (x)[(N) * bfam_kron_k + bfam_kron_i];                  \
  } while (0)

/** $y = (I \otimes A^T) x$.
 *
 * \param [in]  N number of rows and columns of $A$
 * \param [in]  A column first representation of $A$
 * \param [in]  x vector $x$
 * \param [out] y vector $y$
 */
#define BFAM_KRON_IXAT(N, A, x, y)                                 \
  do                                                               \
  {                                                                \
    for(int bfam_kron_n = 0; bfam_kron_n < (N)*(N); ++bfam_kron_n) \
      (y)[bfam_kron_n] = 0;                                        \
                                                                   \
    for(int bfam_kron_k = 0; bfam_kron_k < (N); ++bfam_kron_k)     \
      for(int bfam_kron_j = 0; bfam_kron_j < (N); ++bfam_kron_j)   \
        for(int bfam_kron_i = 0; bfam_kron_i < (N); ++bfam_kron_i) \
          (y)[(N) * bfam_kron_k + bfam_kron_j] +=                  \
            (A)[(N) * bfam_kron_j + bfam_kron_i] *                 \
            (x)[(N) * bfam_kron_k + bfam_kron_i];                  \
  } while (0)

/** $y = (A \otimes I) x$.
 *
 * \param [in]  N number of rows and columns of $A$
 * \param [in]  A column first representation of $A$
 * \param [in]  x vector $x$
 * \param [out] y vector $y$
 */
#define BFAM_KRON_AXI(N, A, x, y)                                  \
  do                                                               \
  {                                                                \
    for(int bfam_kron_n = 0; bfam_kron_n < (N)*(N); ++bfam_kron_n) \
      (y)[bfam_kron_n] = 0;                                        \
                                                                   \
    for(int bfam_kron_i = 0; bfam_kron_i < (N); ++bfam_kron_i)     \
      for(int bfam_kron_k = 0; bfam_kron_k < (N); ++bfam_kron_k)   \
        for(int bfam_kron_j = 0; bfam_kron_j < (N); ++bfam_kron_j) \
          (y)[(N) * bfam_kron_k + bfam_kron_j] +=                  \
            (A)[(N) * bfam_kron_i + bfam_kron_k] *                 \
            (x)[(N) * bfam_kron_i + bfam_kron_j];                  \
  } while (0)

/** $y = (A^T \otimes I) x$.
 *
 * \param [in]  N number of rows and columns of $A$
 * \param [in]  A column first representation of $A$
 * \param [in]  x vector $x$
 * \param [out] y vector $y$
 */
#define BFAM_KRON_ATXI(N, A, x, y)                                 \
  do                                                               \
  {                                                                \
    for(int bfam_kron_n = 0; bfam_kron_n < (N)*(N); ++bfam_kron_n) \
      (y)[bfam_kron_n] = 0;                                        \
                                                                   \
    for(int bfam_kron_k = 0; bfam_kron_k < (N); ++bfam_kron_k)     \
      for(int bfam_kron_i = 0; bfam_kron_i < (N); ++bfam_kron_i)   \
        for(int bfam_kron_j = 0; bfam_kron_j < (N); ++bfam_kron_j) \
          (y)[(N) * bfam_kron_k + bfam_kron_j] +=                  \
            (A)[(N) * bfam_kron_k + bfam_kron_i] *                 \
            (x)[(N) * bfam_kron_i + bfam_kron_j];                  \
  } while (0)

#endif
