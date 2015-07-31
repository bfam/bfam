#include <bfam.h>
#include <bfam_domain_pxest_2.h>
#include <bfam_subdomain_dgx.h>

#define REAL_APPROX_EQ(x, y, K)                                                \
  BFAM_APPROX_EQ((x), (y), (K), BFAM_REAL_ABS, BFAM_REAL_EPS,                  \
                 10 * BFAM_REAL_EPS)

#define DIM 2

#define bfam_domain_pxest_t BFAM_APPEND_EXPAND(bfam_domain_pxest_t_, DIM)
#define bfam_domain_pxest_new BFAM_APPEND_EXPAND(bfam_domain_pxest_new_, DIM)
#define bfam_domain_pxest_init BFAM_APPEND_EXPAND(bfam_domain_pxest_init_, DIM)
#define bfam_domain_pxest_free BFAM_APPEND_EXPAND(bfam_domain_pxest_free_, DIM)
#define bfam_domain_pxest_split_dgx_subdomains                                 \
  BFAM_APPEND_EXPAND(bfam_domain_pxest_split_dgx_subdomains_, DIM)
#define bfam_domain_pxest_create_mesh                                          \
  BFAM_APPEND_EXPAND(bfam_domain_pxest_create_mesh_, DIM)
#define bfam_subdomain_dgx_init_grid                                           \
  BFAM_APPEND_EXPAND(bfam_subdomain_dgx_init_grid_, DIM)

bfam_real_t xex[] = {
    0.0,                1.000000000000000,  0.500000000000000,
    1.299038105676658,  0.0,                0.500000000000000,
    -0.500000000000000, 0.000000000000000,  -0.500000000000000,
    -1.299038105676658, 0.0,                -1.000000000000000,
    -1.299038105676658, -1.000000000000000, -0.500000000000000,
    0.0,                -0.500000000000000, 0.0,
    -0.000000000000000, 0.500000000000000,  1.299038105676658,
    1.000000000000000,  0.500000000000000,  0.0,
};

bfam_real_t yex[] = {
    0.0,                0.0,                0.866025403784439,
    0.750000000000000,  0.0,                0.866025403784439,
    0.866025403784439,  1.500000000000000,  0.866025403784439,
    0.750000000000000,  0.0,                0.000000000000000,
    -0.750000000000000, 0.000000000000000,  -0.866025403784438,
    0.0,                -0.866025403784438, 0.0,
    -1.500000000000000, -0.866025403784439, -0.750000000000001,
    0.0,                -0.866025403784439, 0.0,
};

bfam_real_t Jrxex[] = {
    0.216506350946110 * 2.000000000000000,
    0.187500000000000 * 2.000000000000000,
    0.187500000000000 * 2.309401076758503,
    0.158493649053890 * 2.366025403784439,
    0.216506350946110 * 2.000000000000000,
    0.187500000000000 * 1.690598923241497,
    0.187500000000000 * 2.309401076758504,
    0.158493649053890 * 2.000000000000000,
    0.187500000000000 * -2.309401076758502,
    0.158493649053890 * -2.366025403784438,
    0.216506350946110 * -2.000000000000000,
    0.187500000000000 * -2.000000000000000,
    -0.158493649053890 * 0.366025403784439,
    -0.187500000000000 * 0.000000000000000,
    -0.187500000000000 * 0.309401076758504,
    -0.216506350946110 * 0.000000000000000,
    -0.187500000000000 * 1.690598923241497,
    -0.216506350946110 * 1.999999999999999,
    -0.158493649053890 * 1.999999999999999,
    -0.187500000000000 * 2.309401076758502,
    0.158493649053890 * -0.366025403784436,
    0.187500000000000 * 0.0,
    0.187500000000000 * -0.309401076758501,
    0.216506350946110 * 0.0,
};

bfam_real_t Jryex[] = {
    0.216506350946110 * -1.154700538379252,
    0.187500000000000 * -0.797434948471088,
    0.187500000000000 * -1.333333333333334,
    0.158493649053890 * -0.943375672974065,
    0.216506350946110 * 1.154700538379251,
    0.187500000000000 * 1.333333333333333,
    0.187500000000000 * 1.333333333333333,
    0.158493649053890 * 1.577350269189626,
    0.187500000000000 * -1.333333333333332,
    0.158493649053890 * -0.943375672974064,
    0.216506350946110 * -1.154700538379251,
    0.187500000000000 * -0.797434948471088,
    -0.158493649053890 * 2.520725942163691,
    -0.187500000000000 * 2.666666666666668,
    -0.187500000000000 * 2.130768281804421,
    -0.216506350946110 * 2.309401076758503,
    -0.187500000000000 * 1.333333333333333,
    -0.216506350946110 * 1.154700538379251,
    -0.158493649053890 * 1.577350269189625,
    -0.187500000000000 * 1.333333333333333,
    0.158493649053890 * 2.520725942163689,
    0.187500000000000 * 2.666666666666664,
    0.187500000000000 * 2.130768281804422,
    0.216506350946110 * 2.309401076758503,
};

bfam_real_t Jsxex[] = {
    0.216506350946110 * 0.0,
    0.187500000000000 * 0.0,
    0.187500000000000 * 0.309401076758503,
    0.158493649053890 * 0.366025403784439,
    0.216506350946110 * -2.000000000000000,
    0.187500000000000 * -2.309401076758503,
    0.187500000000000 * -1.690598923241497,
    0.158493649053890 * -2.000000000000000,
    0.187500000000000 * 0.309401076758503,
    0.158493649053890 * 0.366025403784439,
    0.216506350946110 * -0.000000000000000,
    0.187500000000000 * -0.000000000000000,
    -0.158493649053890 * 2.366025403784439,
    -0.187500000000000 * 2.000000000000000,
    -0.187500000000000 * 2.309401076758504,
    -0.216506350946110 * 2.000000000000000,
    -0.187500000000000 * 2.309401076758501,
    -0.216506350946110 * 1.999999999999998,
    -0.158493649053890 * 1.999999999999999,
    -0.187500000000000 * 1.690598923241496,
    0.158493649053890 * -2.366025403784440,
    0.187500000000000 * -2.000000000000000,
    0.187500000000000 * -2.309401076758504,
    0.216506350946110 * -2.000000000000000,
};

bfam_real_t Jsyex[] = {
    0.216506350946110 * 2.309401076758503,
    0.187500000000000 * 2.666666666666667,
    0.187500000000000 * 2.130768281804421,
    0.158493649053890 * 2.520725942163691,
    0.216506350946110 * 1.154700538379252,
    0.187500000000000 * 1.333333333333334,
    0.187500000000000 * 1.333333333333333,
    0.158493649053890 * 1.577350269189626,
    0.187500000000000 * -2.130768281804421,
    0.158493649053890 * -2.520725942163691,
    0.216506350946110 * -2.309401076758503,
    0.187500000000000 * -2.666666666666667,
    -0.158493649053890 * -0.943375672974065,
    -0.187500000000000 * -0.797434948471088,
    -0.187500000000000 * -1.333333333333335,
    -0.216506350946110 * -1.154700538379253,
    -0.187500000000000 * -1.333333333333334,
    -0.216506350946110 * -1.154700538379252,
    -0.158493649053890 * -1.577350269189626,
    -0.187500000000000 * -1.333333333333334,
    0.158493649053890 * -0.943375672974064,
    0.187500000000000 * -0.797434948471087,
    0.187500000000000 * -1.333333333333334,
    0.216506350946110 * -1.154700538379252,
};

bfam_real_t Jex[] = {
    0.216506350946110,  0.187500000000000,  0.187500000000000,
    0.158493649053890,  0.216506350946110,  0.187500000000000,
    0.187500000000000,  0.158493649053890,  0.187500000000000,
    0.158493649053890,  0.216506350946110,  0.187500000000000,
    -0.158493649053890, -0.187500000000000, -0.187500000000000,
    -0.216506350946110, -0.187500000000000, -0.216506350946110,
    -0.158493649053890, -0.187500000000000, 0.158493649053890,
    0.187500000000000,  0.187500000000000,  0.216506350946110,
};

bfam_real_t nxex[] = {
    -0.866025403784439, -0.866025403784439, 0.928886923381516,
    0.928886923381516,  0,                  0,
    0.143699307140567,  0.143699307140567,  -0.866025403784439,
    -0.866025403784439, 0.785187616240949,  0.785187616240949,
    0.866025403784439,  0.866025403784439,  -0.785187616240949,
    -0.785187616240949, 0.866025403784439,  0.866025403784439,
    -0.928886923381516, -0.928886923381516, -0.143699307140567,
    -0.143699307140567, -0.000000000000000, -0.000000000000000,
    0.143699307140567,  0.143699307140567,  -0.000000000000000,
    -0.000000000000000, 0.928886923381516,  0.928886923381516,
    -0.866025403784438, -0.866025403784438, 0.785187616240949,
    0.785187616240949,  -0.866025403784439, -0.866025403784439,
    0.866025403784438,  0.866025403784438,  -0.785187616240948,
    -0.785187616240948, 0.143699307140566,  0.143699307140566,
    0,                  0,                  0.928886923381516,
    0.928886923381516,  -0.866025403784439, -0.866025403784439,
};

bfam_real_t nyex[] = {
    0.500000000000000,  0.500000000000000,  -0.370363447941103,
    -0.370363447941103, -1.000000000000000, -1.000000000000000,
    0.989621396862114,  0.989621396862114,  -0.500000000000000,
    -0.500000000000000, 0.619257948921010,  0.619257948921010,
    -0.500000000000000, -0.500000000000000, 0.619257948921010,
    0.619257948921010,  0.500000000000000,  0.500000000000000,
    -0.370363447941103, -0.370363447941103, 0.989621396862114,
    0.989621396862114,  -1.000000000000000, -1.000000000000000,
    0.989621396862114,  0.989621396862114,  -1.000000000000000,
    -1.000000000000000, -0.370363447941103, -0.370363447941103,
    0.500000000000000,  0.500000000000000,  0.619257948921010,
    0.619257948921010,  -0.500000000000000, -0.500000000000000,
    -0.500000000000000, -0.500000000000000, 0.619257948921011,
    0.619257948921011,  -0.989621396862114, -0.989621396862114,
    1.000000000000000,  1.000000000000000,  0.370363447941103,
    0.370363447941103,  -0.500000000000000, -0.500000000000000,
};

bfam_real_t sJex[] = {
    0.500000000000000, 0.500000000000000, 0.403708988210160, 0.403708988210160,
    0.500000000000000, 0.500000000000000, 0.403708988210160, 0.403708988210160,
    0.500000000000000, 0.500000000000000, 0.403708988210160, 0.403708988210160,
    0.500000000000000, 0.500000000000000, 0.403708988210160, 0.403708988210160,
    0.500000000000000, 0.500000000000000, 0.403708988210160, 0.403708988210160,
    0.403708988210160, 0.403708988210160, 0.500000000000000, 0.500000000000000,
    0.403708988210160, 0.403708988210160, 0.500000000000000, 0.500000000000000,
    0.403708988210160, 0.403708988210160, 0.500000000000000, 0.500000000000000,
    0.403708988210160, 0.403708988210160, 0.500000000000000, 0.500000000000000,
    0.500000000000000, 0.500000000000000, 0.403708988210160, 0.403708988210160,
    0.403708988210160, 0.403708988210160, 0.500000000000000, 0.500000000000000,
    0.403708988210160, 0.403708988210160, 0.500000000000000, 0.500000000000000,
};

/*
 * arbitrary dimension dgx test
 */
static int check_field_face_dgx(bfam_subdomain_t *s, const char *name,
                                bfam_real_t *ex)
{
  bfam_subdomain_dgx_t *sub = (bfam_subdomain_dgx_t *)s;
  int failures = 0;
  bfam_real_t *field =
      bfam_dictionary_get_value_ptr(&sub->base.fields_face, name);

  for (bfam_locidx_t i = 0; i < sub->K; ++i)
  {
    for (int j = 0; j < sub->Ngp[0] * sub->Ng[0]; ++j)
    {
      size_t idx = i * sub->Ngp[0] * sub->Ng[0] + j;

      int fail = !REAL_APPROX_EQ(field[idx], ex[idx], 10);

      if (fail)
        BFAM_INFO("Fail match (%s) %d %25.15" BFAM_REAL_PRIe
                  " %25.15" BFAM_REAL_PRIe " %d",
                  name, idx, field[idx], ex[idx],
                  BFAM_REAL_ABS(field[idx] - ex[idx]) < BFAM_REAL_MIN);

      failures += fail;
    }
  }

  return failures;
}

static int check_field_dgx(bfam_subdomain_t *s, const char *name,
                           bfam_real_t *ex)
{
  bfam_subdomain_dgx_t *sub = (bfam_subdomain_dgx_t *)s;
  int failures = 0;
  bfam_real_t *field = bfam_dictionary_get_value_ptr(&sub->base.fields, name);

  for (bfam_locidx_t i = 0; i < sub->K; ++i)
  {
    for (int j = 0; j < sub->Np; ++j)
    {
      size_t idx = i * sub->Np + j;

      int fail = !REAL_APPROX_EQ(field[idx], ex[idx], 10);

      if (fail)
        BFAM_INFO("Fail match (%s) %d %25.15" BFAM_REAL_PRIe
                  " %25.15" BFAM_REAL_PRIe " %d",
                  name, idx, field[idx], ex[idx],
                  BFAM_REAL_ABS(field[idx] - ex[idx]) < BFAM_REAL_MIN);

      failures += fail;
    }
  }

  return failures;
}

static int test_geo_dgx(MPI_Comm mpicomm)
{
  int failures = 0;
  int rank;
  BFAM_MPI_CHECK(MPI_Comm_rank(mpicomm, &rank));

  p4est_connectivity_t *conn = p4est_connectivity_new_star();

  for (p4est_topidx_t n = 0; n < conn->num_vertices; ++n)
    BFAM_INFO("%jd %25.15e %25.15e", (intmax_t)n + 1, conn->vertices[n * 3 + 0],
              conn->vertices[n * 3 + 1]);

  bfam_domain_pxest_t *domain = bfam_domain_pxest_new(mpicomm, conn);

  p4est_balance(domain->pxest, P4EST_CONNECT_CORNER, NULL);
  p4est_partition(domain->pxest, 1, NULL);

  bfam_locidx_t numSubdomains = 1;
  int N = 1;

  bfam_locidx_t *subdomainID =
      bfam_calloc(domain->pxest->local_num_quadrants, sizeof(bfam_locidx_t));

  bfam_domain_pxest_split_dgx_subdomains(domain, numSubdomains, subdomainID,
                                         NULL, &N, NULL, NULL, NULL);

  bfam_domain_pxest_create_mesh(domain, NULL, NULL);

  bfam_subdomain_t **subdomains =
      bfam_malloc(domain->base.numSubdomains * sizeof(bfam_subdomain_t **));

  const char *volume[] = {"_volume", NULL};
  bfam_domain_get_subdomains((bfam_domain_t *)domain, BFAM_DOMAIN_OR, volume,
                             domain->base.numSubdomains, subdomains,
                             &numSubdomains);

  BFAM_ABORT_IF_NOT(numSubdomains == 1, "We should only have one subdomain");

  failures += check_field_dgx(subdomains[0], "_grid_x0", xex);
  failures += check_field_dgx(subdomains[0], "_grid_x1", yex);

  failures += check_field_dgx(subdomains[0], "_grid_Jr0x0", Jrxex);
  failures += check_field_dgx(subdomains[0], "_grid_Jr0x1", Jryex);

  failures += check_field_dgx(subdomains[0], "_grid_Jr1x0", Jsxex);
  failures += check_field_dgx(subdomains[0], "_grid_Jr1x1", Jsyex);

  failures += check_field_dgx(subdomains[0], "_grid_J", Jex);

  failures += check_field_face_dgx(subdomains[0], "_grid_nx0", nxex);
  failures += check_field_face_dgx(subdomains[0], "_grid_nx1", nyex);

  bfam_free(subdomainID);
  bfam_free(subdomains);

  bfam_domain_pxest_free(domain);
  bfam_free(domain);

  p4est_connectivity_destroy(conn);
  return failures;
}

int main(int argc, char *argv[])
{
  int failures = 0;
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank;

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
      "  there are four possible options to this program, some of which have\n"
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

  BFAM_MPI_CHECK(MPI_Init(&argc, &argv));
  BFAM_MPI_CHECK(MPI_Comm_rank(comm, &rank));

  int logLevel = BFAM_MAX(BFAM_LL_INFO - verbosity, BFAM_LL_ALWAYS);

  bfam_log_init(rank, stdout, logLevel);
  bfam_signal_handler_set();

  int scLogPriorities = BFAM_MAX(SC_LP_STATISTICS - verbosity, SC_LP_ALWAYS);
  sc_init(comm, 0, 0, NULL, scLogPriorities);
  p4est_init(NULL, scLogPriorities);

  BFAM_MPI_CHECK(MPI_Barrier(comm));
  BFAM_ROOT_INFO("Test Geo DGX\n");
  failures += test_geo_dgx(comm);
  BFAM_MPI_CHECK(MPI_Barrier(comm));

  sc_finalize();
  BFAM_MPI_CHECK(MPI_Finalize());

  bfam_gopt_free(options);

  return failures;
}
