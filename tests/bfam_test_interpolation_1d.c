#include <bfam.h>
#include <bfam_base.h>

#define REAL_APPROX_EQ(x, y, K)                                                \
  BFAM_APPROX_EQ((x), (y), (K), BFAM_REAL_ABS, BFAM_REAL_EPS, (K)*BFAM_REAL_EPS)

typedef bfam_real_t (*polyval_t)(bfam_long_real_t x, bfam_long_real_t y,
                                 bfam_long_real_t z);

static bfam_real_t poly1_field(bfam_long_real_t x, bfam_long_real_t y,
                               bfam_long_real_t z)
{
  return -x - y * x * y;
}

static bfam_real_t poly3_field(bfam_long_real_t x, bfam_long_real_t y,
                               bfam_long_real_t z)
{
  return -x + x * x - 0.5 * x * x * x;
}

static int test_1d(bfam_dictionary_t *N2N, int NC, int *NF, polyval_t poly)
{
  int failures = 0;

  /* Set up the coarse data */
  bfam_long_real_t lrC[NC + 1];
  bfam_long_real_t lwC[NC + 1];
  bfam_jacobi_gauss_lobatto_quadrature(0, 0, NC, lrC, lwC);
  bfam_real_t pC[(NC + 1)];
  bfam_real_t pC_trg[(NC + 1)];
  for (bfam_locidx_t ix = 0; ix <= NC; ix++)
  {
    pC[ix] = poly(lrC[ix], 0, 0);
    pC_trg[ix] = 0;
  }

  /* Set up the fine data */
  bfam_long_real_t *lrF[2];
  bfam_long_real_t *lwF[2];
  bfam_real_t *pF[2];
  bfam_real_t *pF_trg[2];
  const bfam_long_real_t x_off[2] = {BFAM_LONG_REAL(-1.0), BFAM_LONG_REAL(1.0)};
  const bfam_long_real_t HALF = BFAM_LONG_REAL(0.5);
  for (bfam_locidx_t k = 0; k < 2; k++)
  {
    pF[k] = bfam_malloc_aligned((NF[k] + 1) * sizeof(bfam_real_t));
    pF_trg[k] = bfam_malloc_aligned((NF[k] + 1) * sizeof(bfam_real_t));
    lrF[k] = bfam_malloc_aligned((NF[k] + 1) * sizeof(bfam_long_real_t));
    lwF[k] = bfam_malloc_aligned((NF[k] + 1) * sizeof(bfam_long_real_t));
    bfam_jacobi_gauss_lobatto_quadrature(0, 0, NF[k], lrF[k], lwF[k]);
    for (bfam_locidx_t ix = 0; ix <= NF[k]; ix++)
    {
      const bfam_long_real_t x = HALF * (lrF[k][ix] + x_off[k]);
      pF[k][ix] = poly(x, 0, 0);
      pF_trg[k][ix] = 0;
    }

    bfam_subdomain_dgx_interpolate_data_1(pF[k], k, NF[k], pC_trg, 0, NC,
                                          BFAM_FLAG_COARSEN, N2N, 1);
    bfam_subdomain_dgx_interpolate_data_1(pC, 0, NC, pF_trg[k], k, NF[k],
                                          BFAM_FLAG_REFINE, N2N, 1);

    for (bfam_locidx_t n = 0; n < (NF[k] + 1); n++)
      if (!REAL_APPROX_EQ(pF[k][n], pF_trg[k][n], 1000))
      {
        failures++;
        BFAM_VERBOSE("pF[%d][%d]: got %" BFAM_REAL_PRIe
                     " and expected %" BFAM_REAL_PRIe,
                     (int)k, (int)n, pF_trg[k][n], pF[k][n]);
      }
  }

  /* check the interpolation */
  for (bfam_locidx_t n = 0; n < (NC + 1); n++)
    if (!REAL_APPROX_EQ(pC[n], pC_trg[n], 1000))
    {
      failures++;
      BFAM_VERBOSE("pC[%d]: got %" BFAM_REAL_PRIe
                   " and expected %" BFAM_REAL_PRIe,
                   (int)n, pC_trg[n], pC[n]);
    }

  for (bfam_locidx_t k = 0; k < 2; k++)
  {
    bfam_free_aligned(pF[k]);
    bfam_free_aligned(pF_trg[k]);
    bfam_free_aligned(lrF[k]);
    bfam_free_aligned(lwF[k]);
  }

  return failures;
}

int main(int argc, char *argv[])
{
  int failures = 0;

  void *options = bfam_gopt_sort(
      &argc, (const char **)argv,
      bfam_gopt_start(bfam_gopt_option('h', 0, bfam_gopt_shorts('h', '?'),
                                       bfam_gopt_longs("help", "HELP")),
                      bfam_gopt_option('V', 0, bfam_gopt_shorts('V'),
                                       bfam_gopt_longs("version")),
                      bfam_gopt_option('v', BFAM_GOPT_REPEAT,
                                       bfam_gopt_shorts('v'),
                                       bfam_gopt_longs("verbose"))));

  const char *helpText =
      "\n"
      "\n"
      "  there are four possible options to this program, some of which "
      "have\n"
      "  multiple names:\n"
      "\n"
      "    -h -? --help --HELP\n"
      "    -V --version\n"
      "    -v --verbose  (which may be repeated for more verbosity)\n"
      "\n";

  if (bfam_gopt(options, 'h'))
  {
    /*
     * if any of the help options was specified
     */
    BFAM_ROOT_INFO(helpText);
    exit(EXIT_SUCCESS);
  }

  if (bfam_gopt(options, 'V'))
  {
    BFAM_ROOT_INFO("BFAM Version: %s", bfam_version_get());
    BFAM_ROOT_INFO("BFAM Compile Info:\n" BFAM_COMPILE_INFO);
    exit(EXIT_SUCCESS);
  }

  int verbosity = bfam_gopt(options, 'v');

  int logLevel = BFAM_MAX(BFAM_LL_INFO - verbosity, BFAM_LL_ALWAYS);

  bfam_log_init(0, stdout, logLevel);
  bfam_signal_handler_set();

  bfam_dictionary_t N2N;
  bfam_dictionary_init(&N2N);

  // int NF[] = {4, 5, 6, 3};
  // failures += test_coarsen_2d(&N2N, 4, NF, poly1_field);

  int NF[] = {5, 7};
  for (int k = 3; k < 9; k++)
  {
    failures += test_1d(&N2N, k, NF, poly1_field);
    failures += test_1d(&N2N, k, NF, poly3_field);
  }

  bfam_dictionary_allprefixed_ptr(
      &N2N, "", bfam_subdomain_dgx_clear_interpolation_dict_1, NULL);
  bfam_dictionary_clear(&N2N);

  bfam_gopt_free(options);

  if (failures > 0)
    BFAM_WARNING("FAIL! with failures %d", failures);
  return failures;
}
