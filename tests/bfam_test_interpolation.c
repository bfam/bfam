#include <bfam.h>
#include <bfam_base.h>

#define REAL_APPROX_EQ(x, y, K)                                                \
  BFAM_APPROX_EQ((x), (y), (K), BFAM_REAL_ABS, BFAM_REAL_EPS, (K)*BFAM_REAL_EPS)

typedef struct
{
  int N_src;
  int N_dst;
  bfam_locidx_t num_prj;
  bfam_real_t **prj; /* array of projection operators;
                      * no refinement (NULL is no change in order)
                      * coarsen from 0
                      * coarsen from 1
                      * refine  to   0
                      * refine  to   1
                      */
} interpolator_t;

typedef enum
{
  COARSEN,
  SAME,
  REFINE,
} crf_t;

static void init_interpolator(interpolator_t *interp_a2b, const int N_a,
                              const int N_b)
{
  const bfam_locidx_t num_prj = 5;
  const bfam_locidx_t Np_a = N_a + 1;
  const bfam_locidx_t Np_b = N_b + 1;

  interp_a2b->N_src = N_a;
  interp_a2b->N_dst = N_b;
  interp_a2b->num_prj = num_prj;
  interp_a2b->prj = bfam_malloc(num_prj * sizeof(bfam_real_t *));

  /* Storage for the projection operators */
  for (bfam_locidx_t k = 0; k < num_prj; k++)
  {
    if (k == 0 && N_a == N_b)
      interp_a2b->prj[k] = NULL;
    interp_a2b->prj[k] = bfam_malloc_aligned(Np_a * Np_b * sizeof(bfam_real_t));
  }
}

static void fill_grid_data(bfam_long_real_t N, bfam_long_real_t Np,
                           bfam_long_real_t *lr, bfam_long_real_t *lw,
                           bfam_long_real_t *V, bfam_long_real_t *M)
{
  bfam_jacobi_gauss_lobatto_quadrature(0, 0, Np, lr, lw);
  bfam_jacobi_p_vandermonde(0, 0, N, Np, lr, V);
  bfam_jacobi_p_mass(0, 0, N, V, M);
}

/* Build interpolation and projection between space a and g. Space a is the
 * lower order space and space g is the higher order space.
 */
static void fill_interp_proj_data(bfam_locidx_t N_a, bfam_locidx_t N_g,
                                  bfam_long_real_t *V_a, bfam_long_real_t *M_a,
                                  bfam_long_real_t *M_g, bfam_long_real_t *lr_g,
                                  bfam_long_real_t *I_a2g,
                                  bfam_long_real_t *P_g2a)
{
  BFAM_ASSERT(N_g >= N_a);
  bfam_locidx_t Np_a = N_a + 1;
  bfam_locidx_t Np_g = N_g + 1;

  bfam_jacobi_p_interpolation(0, 0, N_a, Np_g, lr_g, V_a, I_a2g);

  bfam_long_real_t MP_g2a[Np_g * Np_a];
  bfam_long_real_t Vt_MP_g2a[Np_g * Np_a];
  for (bfam_locidx_t n = 0; n < Np_g * Np_a; n++)
  {
    MP_g2a[n] = 0;
    Vt_MP_g2a[n] = 0;
  }

  /* MP_g2a = M_a * P_g2a = I_a2g^{T} * M_g */
  /* [Np_a X Np_g] = [Np_a X Np_g] [Np_g X Np_g] */
  bfam_util_mTmmult(Np_a, Np_g, Np_g, I_a2g, Np_g, M_g, Np_g, MP_g2a, Np_a);

  /* Vt_MP_g2a = V_a^{-1} * P_g2a = V_a^{T} MP_g2a */
  /* [Np_a X Np_g] = [Np_a X Np_a] [Np_a X Np_g] */
  bfam_util_mTmmult(Np_a, Np_g, Np_a, V_a, Np_a, MP_g2a, Np_g, Vt_MP_g2a, Np_a);

  /* P_g2a = V_a * V_a^{T} I_a2g^{T} * M_g */
  /* [Np_a X Np_g] = [Np_a X Np_a] [Np_a X Np_g] */
  bfam_util_mmmult(Np_a, Np_g, Np_a, V_a, Np_a, Vt_MP_g2a, Np_g, P_g2a, Np_a);
}

static void fill_hanging_data(bfam_locidx_t N, bfam_locidx_t Np,
                              bfam_long_real_t *lr_f, bfam_long_real_t *M,
                              bfam_long_real_t *V, bfam_long_real_t *I_f2t,
                              bfam_long_real_t *I_f2b, bfam_long_real_t *P_t2f,
                              bfam_long_real_t *P_b2f)
{
  const bfam_long_real_t HALF = BFAM_LONG_REAL(0.5);

  /* Grid for the top and bottom */
  bfam_long_real_t lr_t[Np];
  bfam_long_real_t lr_b[Np];

  for (bfam_locidx_t k = 0; k < Np; k++)
  {
    lr_t[k] = HALF * (lr_f[k] + 1);
    lr_b[k] = HALF * (lr_f[k] - 1);
  }

  /* Mass matrix for the top and bottom */
  bfam_long_real_t Mh[Np * Np];
  for (bfam_locidx_t k = 0; k < Np * Np; k++)
    Mh[k] = HALF * M[k];

  fill_interp_proj_data(N, N, V, M, Mh, lr_t, I_f2t, P_t2f);
  fill_interp_proj_data(N, N, V, M, Mh, lr_b, I_f2b, P_b2f);
}

static void create_interpolators(interpolator_t *interp_a2b,
                                 interpolator_t *interp_b2a, const int N_a,
                                 const int N_b)
{
  /* Figure out the glue space and the number of points */
  const int N_g = BFAM_MAX(N_a, N_b);
  const int Np_g = N_g + 1;
  const bfam_locidx_t Np_a = N_a + 1;
  const bfam_locidx_t Np_b = N_b + 1;

  /* initialize the interpolators */
  init_interpolator(interp_a2b, N_a, N_b);
  if (interp_b2a)
    init_interpolator(interp_b2a, N_b, N_a);

  /* Set up the reference grids for the side a */
  bfam_long_real_t lr_a[Np_a];
  bfam_long_real_t lw_a[Np_a];
  bfam_long_real_t V_a[Np_a * Np_a];
  bfam_long_real_t M_a[Np_a * Np_a];
  fill_grid_data(N_a, Np_a, lr_a, lw_a, V_a, M_a);

  /* Set up the reference grids for the side b */
  bfam_long_real_t lr_b[Np_b];
  bfam_long_real_t lw_b[Np_b];
  bfam_long_real_t V_b[Np_b * Np_b];
  bfam_long_real_t M_b[Np_b * Np_b];
  fill_grid_data(N_b, Np_b, lr_b, lw_b, V_b, M_b);

  /* Set up the reference grids for the glue space */
  bfam_long_real_t lr_g[Np_g];
  bfam_long_real_t lw_g[Np_g];
  bfam_long_real_t V_g[Np_g * Np_g];
  bfam_long_real_t M_g[Np_g * Np_g];
  fill_grid_data(N_g, Np_g, lr_g, lw_g, V_g, M_g);

  /* Interpolate to the hanging faces */
  bfam_long_real_t I_f2t[Np_g * Np_g]; /* full to top */
  bfam_long_real_t I_f2b[Np_g * Np_g]; /* full to btm */
  bfam_long_real_t P_t2f[Np_g * Np_g]; /* top  to full*/
  bfam_long_real_t P_b2f[Np_g * Np_g]; /* btm  to full*/
  fill_hanging_data(N_g, Np_g, lr_g, M_g, V_g, I_f2t, I_f2b, P_t2f, P_b2f);

  /* Interpolate to the intermediate space and projection back */
  bfam_long_real_t *I_a2g = NULL;
  bfam_long_real_t *P_g2a = NULL;
  if (N_a < N_g)
  {
    I_a2g = bfam_malloc_aligned(Np_g * Np_a * sizeof(bfam_long_real_t));
    P_g2a = bfam_malloc_aligned(Np_g * Np_a * sizeof(bfam_long_real_t));
    fill_interp_proj_data(N_a, N_g, V_a, M_a, M_g, lr_g, I_a2g, P_g2a);
  }

  bfam_long_real_t *I_b2g = NULL;
  bfam_long_real_t *P_g2b = NULL;
  if (N_b < N_g)
  {
    I_b2g = bfam_malloc_aligned(Np_g * Np_b * sizeof(bfam_long_real_t));
    P_g2b = bfam_malloc_aligned(Np_g * Np_b * sizeof(bfam_long_real_t));
    fill_interp_proj_data(N_b, N_g, V_b, M_b, M_g, lr_g, I_b2g, P_g2b);
  }

  /* free the temporary storage */
  if (I_a2g)
    bfam_free_aligned(I_a2g);
  if (P_g2a)
    bfam_free_aligned(P_g2a);
  if (I_b2g)
    bfam_free_aligned(I_b2g);
  if (P_g2b)
    bfam_free_aligned(P_g2b);
}

static void interp_2D(const bfam_real_t *src, const int N_src, bfam_real_t *dst,
                      const int N_dst, const crf_t crf,
                      const bfam_locidx_t child, bfam_dictionary_t *N2N)
{
  char str[BFAM_BUFSIZ];
  snprintf(str, BFAM_BUFSIZ, "%d_to_%d", N_src, N_dst);

  interpolator_t *interp =
      (interpolator_t *)bfam_dictionary_get_value_ptr(N2N, str);

  if (!interp)
  {
    interp = bfam_malloc(sizeof(interpolator_t));
    interpolator_t *interp2 = NULL;
    if (N_src != N_dst)
      interp2 = bfam_malloc(sizeof(interpolator_t));
    create_interpolators(interp, interp2, N_src, N_dst);
    int rval = bfam_dictionary_insert_ptr(N2N, str, interp);
    BFAM_ASSERT(rval != 1);

    if (interp2)
    {
      char str2[BFAM_BUFSIZ];
      snprintf(str2, BFAM_BUFSIZ, "%d_to_%d", N_dst, N_src);
      rval = bfam_dictionary_insert_ptr(N2N, str2, interp2);
      BFAM_ASSERT(rval != 1);
    }
  }

  switch (crf)
  {
  case COARSEN:
    break;
  case SAME:
    break;
  case REFINE:
    break;
  }
}

typedef bfam_real_t (*polyval_t)(bfam_long_real_t x, bfam_long_real_t y);

static bfam_real_t poly1_field(bfam_long_real_t x, bfam_long_real_t y)
{
  return -x - y * x * y;
}

static int test_coarsen_2d(bfam_dictionary_t *N2N, int NC, int *NF,
                           polyval_t poly)
{
  int failures = 0;

  /* Set up the coarse data */
  bfam_long_real_t lrC[NC + 1];
  bfam_long_real_t lwC[NC + 1];
  bfam_jacobi_gauss_lobatto_quadrature(0, 0, NC, lrC, lwC);
  bfam_real_t pC[(NC + 1) * (NC + 1)];
  bfam_real_t pC_trg[(NC + 1) * (NC + 1)];
  for (bfam_locidx_t iy = 0; iy <= NC; iy++)
    for (bfam_locidx_t ix = 0; ix <= NC; ix++)
    {
      pC[ix + iy * (NC + 1)] = poly(lrC[ix], lrC[iy]);
      pC_trg[ix + iy * (NC + 1)] = 0;
    }

  /* Set up the fine data */
  bfam_long_real_t *lrF[4];
  bfam_long_real_t *lwF[4];
  bfam_real_t *pF[4];
  const bfam_long_real_t x_off[4] = {
      BFAM_LONG_REAL(-1.0), BFAM_LONG_REAL(1.0), BFAM_LONG_REAL(-1.0),
      BFAM_LONG_REAL(1.0),
  };
  const bfam_long_real_t y_off[4] = {
      BFAM_LONG_REAL(-1.0), BFAM_LONG_REAL(-1.0), BFAM_LONG_REAL(1.0),
      BFAM_LONG_REAL(1.0),
  };
  const bfam_long_real_t HALF = BFAM_LONG_REAL(0.5);

  for (bfam_locidx_t k = 0; k < 4; k++)
  {
    pF[k] =
        bfam_malloc_aligned((NF[k] + 1) * (NF[k] + 1) * sizeof(bfam_real_t));
    lrF[k] = bfam_malloc_aligned((NF[k] + 1) * sizeof(bfam_long_real_t));
    lwF[k] = bfam_malloc_aligned((NF[k] + 1) * sizeof(bfam_long_real_t));
    bfam_jacobi_gauss_lobatto_quadrature(0, 0, NF[k], lrF[k], lwF[k]);
    for (bfam_locidx_t iy = 0; iy <= NF[k]; iy++)
      for (bfam_locidx_t ix = 0; ix <= NF[k]; ix++)
      {
        const bfam_long_real_t x = HALF * (lrF[k][ix] + x_off[k]);
        const bfam_long_real_t y = HALF * (lrF[k][iy] + y_off[k]);
        pF[k][ix + iy * (NF[k] + 1)] = poly(x, y);
      }

    interp_2D(pF[k], NF[k], pC_trg, NC, COARSEN, k, N2N);
  }

  /* check the interpolation */
  for (bfam_locidx_t k = 0; k < (NC + 1) * (NC + 1); k++)
    failures += !REAL_APPROX_EQ(pC[k], pC_trg[k], 1000);

  for (bfam_locidx_t k = 0; k < 4; k++)
  {
    bfam_free_aligned(pF[k]);
    bfam_free_aligned(lrF[k]);
    bfam_free_aligned(lwF[k]);
  }
  return failures;
}

static int clear_N2N_dict(const char *key, void *val, void *args)
{
  interpolator_t *interp = (interpolator_t *)val;
  for (bfam_locidx_t k = 0; k < interp->num_prj; k++)
    if (interp->prj[k])
      bfam_free_aligned(interp->prj[k]);
  bfam_free(interp->prj);
  bfam_free(val);
  return 1;
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

  int NF[] = {4, 5, 6, 3};
  bfam_dictionary_t N2N;
  bfam_dictionary_init(&N2N);
  failures += test_coarsen_2d(&N2N, 4, NF, poly1_field);

  bfam_dictionary_allprefixed_ptr(&N2N, "", clear_N2N_dict, NULL);
  bfam_dictionary_clear(&N2N);

  bfam_gopt_free(options);

  if (failures > 0)
    BFAM_WARNING("FAIL! with failures %d", failures);
  return failures;
}
