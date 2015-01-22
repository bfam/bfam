#include <bfam.h>

BFAM_ALIGN(32) int A[] = {
    35, 3,  31, 8,  30, 4,  1,  32, 9,  28, 5,  36, 6,  7,  2,  33, 34, 29,
    26, 21, 22, 17, 12, 13, 19, 23, 27, 10, 14, 18, 24, 25, 20, 15, 16, 11,
};

BFAM_ALIGN(32) int IXA_x[] = {
    398,  437,  398,  371,  356,  371,  1064, 1103, 1064, 1037, 1022, 1037,
    1730, 1769, 1730, 1703, 1688, 1703, 2396, 2435, 2396, 2369, 2354, 2369,
    3062, 3101, 3062, 3035, 3020, 3035, 3728, 3767, 3728, 3701, 3686, 3701,
};

BFAM_ALIGN(32) int IXAT_x[] = {
    340,  445,  502,  340,  364,  340,  1006, 1111, 1168, 1006, 1030, 1006,
    1672, 1777, 1834, 1672, 1696, 1672, 2338, 2443, 2500, 2338, 2362, 2338,
    3004, 3109, 3166, 3004, 3028, 3004, 3670, 3775, 3832, 3670, 3694, 3670,
};

BFAM_ALIGN(32) int AXI_x[] = {
    1833, 1944, 2055, 2166, 2277, 2388, 2067, 2178, 2289, 2400, 2511, 2622,
    1833, 1944, 2055, 2166, 2277, 2388, 1671, 1782, 1893, 2004, 2115, 2226,
    1581, 1692, 1803, 1914, 2025, 2136, 1671, 1782, 1893, 2004, 2115, 2226,
};

BFAM_ALIGN(32) int ATXI_x[] = {
    1485, 1596, 1707, 1818, 1929, 2040, 2115, 2226, 2337, 2448, 2559, 2670,
    2457, 2568, 2679, 2790, 2901, 3012, 1485, 1596, 1707, 1818, 1929, 2040,
    1629, 1740, 1851, 1962, 2073, 2184, 1485, 1596, 1707, 1818, 1929, 2040,
};

BFAM_ALIGN(32) int Ax[] = {
    35,  6,   93,  32,  150, 24,  7,   256, 81,  280, 55,  432,
    78,  98,  30,  528, 578, 522, 494, 420, 462, 374, 276, 312,
    475, 598, 729, 280, 406, 540, 744, 800, 660, 510, 560, 396,
};

int check(int N, int *x, int *y)
{
  int failures = 0;

  for (int n = 0; n < N; ++n)
    failures += (x[n] == y[n]) ? 0 : 1;

  return failures;
}

int main(int argc, char *argv[])
{
  int failures = 0;

  const int N = 6;

  int *x = bfam_malloc_aligned(N * N * sizeof(int));
  int *y = bfam_malloc_aligned(N * N * sizeof(int));

  for (int n = 0; n < N * N; ++n)
    x[n] = n + 1;

  BFAM_KRON_IXA(N, A, x, y);
  failures += check(N, y, IXA_x);

  BFAM_KRON_IXAT(N, A, x, y);
  failures += check(N, y, IXAT_x);

  BFAM_KRON_AXI(N, A, x, y);
  failures += check(N, y, AXI_x);

  BFAM_KRON_ATXI(N, A, x, y);
  failures += check(N, y, ATXI_x);

  BFAM_DOT_AX(N, A, x, y);
  failures += check(N, y, Ax);

  bfam_free_aligned(x);
  bfam_free_aligned(y);

  if (failures)
    return EXIT_FAILURE;
  else
    return EXIT_SUCCESS;
}
