#include <bfam.h>

#define REAL_APPROX_EQ(x, y, K)                                              \
  BFAM_APPROX_EQ((x), (y), (K), BFAM_REAL_ABS, BFAM_REAL_EPS, 10*BFAM_REAL_EPS)

bfam_real_t xex[] = {
                 0.0,
   1.000000000000000,
   0.500000000000000,
   1.299038105676658,
                 0.0,
   0.500000000000000,
  -0.500000000000000,
   0.000000000000000,
  -0.500000000000000,
  -1.299038105676658,
                 0.0,
  -1.000000000000000,
  -1.299038105676658,
  -1.000000000000000,
  -0.500000000000000,
                 0.0,
  -0.500000000000000,
                 0.0,
  -0.000000000000000,
   0.500000000000000,
   1.299038105676658,
   1.000000000000000,
   0.500000000000000,
                 0.0,
};

bfam_real_t yex[] = {
                 0.0,
                 0.0,
   0.866025403784439,
   0.750000000000000,
                 0.0,
   0.866025403784439,
   0.866025403784439,
   1.500000000000000,
   0.866025403784439,
   0.750000000000000,
                 0.0,
   0.000000000000000,
  -0.750000000000000,
   0.000000000000000,
  -0.866025403784438,
                 0.0,
  -0.866025403784438,
                 0.0,
  -1.500000000000000,
  -0.866025403784439,
  -0.750000000000001,
                 0.0,
  -0.866025403784439,
                 0.0,
};

bfam_real_t rxex[] = {
   2.000000000000000,
   2.000000000000000,
   2.309401076758503,
   2.366025403784439,
   2.000000000000000,
   1.690598923241497,
   2.309401076758504,
   2.000000000000000,
  -2.309401076758502,
  -2.366025403784438,
  -2.000000000000000,
  -2.000000000000000,
   0.366025403784439,
   0.000000000000000,
   0.309401076758504,
   0.000000000000000,
   1.690598923241497,
   1.999999999999999,
   1.999999999999999,
   2.309401076758502,
  -0.366025403784436,
                 0.0,
  -0.309401076758501,
                 0.0,
};

bfam_real_t ryex[] = {
  -1.154700538379252,
  -0.797434948471088,
  -1.333333333333334,
  -0.943375672974065,
   1.154700538379251,
   1.333333333333333,
   1.333333333333333,
   1.577350269189626,
  -1.333333333333332,
  -0.943375672974064,
  -1.154700538379251,
  -0.797434948471088,
   2.520725942163691,
   2.666666666666668,
   2.130768281804421,
   2.309401076758503,
   1.333333333333333,
   1.154700538379251,
   1.577350269189625,
   1.333333333333333,
   2.520725942163689,
   2.666666666666664,
   2.130768281804422,
   2.309401076758503,
};

bfam_real_t sxex[] = {
                 0.0,
                 0.0,
   0.309401076758503,
   0.366025403784439,
  -2.000000000000000,
  -2.309401076758503,
  -1.690598923241497,
  -2.000000000000000,
   0.309401076758503,
   0.366025403784439,
  -0.000000000000000,
  -0.000000000000000,
   2.366025403784439,
   2.000000000000000,
   2.309401076758504,
   2.000000000000000,
   2.309401076758501,
   1.999999999999998,
   1.999999999999999,
   1.690598923241496,
  -2.366025403784440,
  -2.000000000000000,
  -2.309401076758504,
  -2.000000000000000,
};

bfam_real_t syex[] = {
   2.309401076758503,
   2.666666666666667,
   2.130768281804421,
   2.520725942163691,
   1.154700538379252,
   1.333333333333334,
   1.333333333333333,
   1.577350269189626,
  -2.130768281804421,
  -2.520725942163691,
  -2.309401076758503,
  -2.666666666666667,
  -0.943375672974065,
  -0.797434948471088,
  -1.333333333333335,
  -1.154700538379253,
  -1.333333333333334,
  -1.154700538379252,
  -1.577350269189626,
  -1.333333333333334,
  -0.943375672974064,
  -0.797434948471087,
  -1.333333333333334,
  -1.154700538379252,
};

bfam_real_t Jex[] = {
   0.216506350946110,
   0.187500000000000,
   0.187500000000000,
   0.158493649053890,
   0.216506350946110,
   0.187500000000000,
   0.187500000000000,
   0.158493649053890,
   0.187500000000000,
   0.158493649053890,
   0.216506350946110,
   0.187500000000000,
  -0.158493649053890,
  -0.187500000000000,
  -0.187500000000000,
  -0.216506350946110,
  -0.187500000000000,
  -0.216506350946110,
  -0.158493649053890,
  -0.187500000000000,
   0.158493649053890,
   0.187500000000000,
   0.187500000000000,
   0.216506350946110,
};

bfam_real_t nxex[] = {
  -0.866025403784439,
  -0.866025403784439,
   0.928886923381516,
   0.928886923381516,
                   0,
                   0,
   0.143699307140567,
   0.143699307140567,
  -0.866025403784439,
  -0.866025403784439,
   0.785187616240949,
   0.785187616240949,
   0.866025403784439,
   0.866025403784439,
  -0.785187616240949,
  -0.785187616240949,
   0.866025403784439,
   0.866025403784439,
  -0.928886923381516,
  -0.928886923381516,
  -0.143699307140567,
  -0.143699307140567,
  -0.000000000000000,
  -0.000000000000000,
   0.143699307140567,
   0.143699307140567,
  -0.000000000000000,
  -0.000000000000000,
   0.928886923381516,
   0.928886923381516,
  -0.866025403784438,
  -0.866025403784438,
   0.785187616240949,
   0.785187616240949,
  -0.866025403784439,
  -0.866025403784439,
   0.866025403784438,
   0.866025403784438,
  -0.785187616240948,
  -0.785187616240948,
   0.143699307140566,
   0.143699307140566,
                   0,
                   0,
   0.928886923381516,
   0.928886923381516,
  -0.866025403784439,
  -0.866025403784439,
};

bfam_real_t nyex[] = {
   0.500000000000000,
   0.500000000000000,
  -0.370363447941103,
  -0.370363447941103,
  -1.000000000000000,
  -1.000000000000000,
   0.989621396862114,
   0.989621396862114,
  -0.500000000000000,
  -0.500000000000000,
   0.619257948921010,
   0.619257948921010,
  -0.500000000000000,
  -0.500000000000000,
   0.619257948921010,
   0.619257948921010,
   0.500000000000000,
   0.500000000000000,
  -0.370363447941103,
  -0.370363447941103,
   0.989621396862114,
   0.989621396862114,
  -1.000000000000000,
  -1.000000000000000,
   0.989621396862114,
   0.989621396862114,
  -1.000000000000000,
  -1.000000000000000,
  -0.370363447941103,
  -0.370363447941103,
   0.500000000000000,
   0.500000000000000,
   0.619257948921010,
   0.619257948921010,
  -0.500000000000000,
  -0.500000000000000,
  -0.500000000000000,
  -0.500000000000000,
   0.619257948921011,
   0.619257948921011,
  -0.989621396862114,
  -0.989621396862114,
   1.000000000000000,
   1.000000000000000,
   0.370363447941103,
   0.370363447941103,
  -0.500000000000000,
  -0.500000000000000,
};

bfam_real_t sJex[] = {
   0.500000000000000,
   0.500000000000000,
   0.403708988210160,
   0.403708988210160,
   0.500000000000000,
   0.500000000000000,
   0.403708988210160,
   0.403708988210160,
   0.500000000000000,
   0.500000000000000,
   0.403708988210160,
   0.403708988210160,
   0.500000000000000,
   0.500000000000000,
   0.403708988210160,
   0.403708988210160,
   0.500000000000000,
   0.500000000000000,
   0.403708988210160,
   0.403708988210160,
   0.403708988210160,
   0.403708988210160,
   0.500000000000000,
   0.500000000000000,
   0.403708988210160,
   0.403708988210160,
   0.500000000000000,
   0.500000000000000,
   0.403708988210160,
   0.403708988210160,
   0.500000000000000,
   0.500000000000000,
   0.403708988210160,
   0.403708988210160,
   0.500000000000000,
   0.500000000000000,
   0.500000000000000,
   0.500000000000000,
   0.403708988210160,
   0.403708988210160,
   0.403708988210160,
   0.403708988210160,
   0.500000000000000,
   0.500000000000000,
   0.403708988210160,
   0.403708988210160,
   0.500000000000000,
   0.500000000000000,
};

static int
check_field_face(bfam_subdomain_t *s, const char *name, bfam_real_t *ex)
{

  bfam_subdomain_dgx_quad_t *sub = (bfam_subdomain_dgx_quad_t *) s;
  int failures = 0;
  bfam_real_t *field =
    bfam_dictionary_get_value_ptr(&sub->base.fields_face, name);

  for(bfam_locidx_t i = 0; i < sub->K; ++i)
  {
    for(int j=0; j<sub->Nfp*sub->Nfaces; ++j)
    {
      size_t idx = i*sub->Nfp*sub->Nfaces + j;

      int fail = !REAL_APPROX_EQ(field[idx], ex[idx], 10);

      if(fail)
        BFAM_LDEBUG("Fail match (%s) %25.15"BFAM_REAL_PRIe
            " %25.15"BFAM_REAL_PRIe " %d", name,
            field[idx], ex[idx],
            BFAM_REAL_ABS(field[idx]-ex[idx]) < BFAM_REAL_MIN);

      failures += fail;
    }
  }

  return failures;
}

static int
check_field(bfam_subdomain_t *s, const char *name, bfam_real_t *ex)
{

  bfam_subdomain_dgx_quad_t *sub = (bfam_subdomain_dgx_quad_t *) s;
  int failures = 0;
  bfam_real_t *field =
    bfam_dictionary_get_value_ptr(&sub->base.fields, name);

  for(bfam_locidx_t i = 0; i < sub->K; ++i)
  {
    for(int j=0; j<sub->Np; ++j)
    {
      size_t idx = i*sub->Np + j;

      int fail = !REAL_APPROX_EQ(field[idx], ex[idx], 10);

      if(fail)
        BFAM_LDEBUG("Fail match (%s) %25.15"BFAM_REAL_PRIe
            " %25.15"BFAM_REAL_PRIe " %d", name,
            field[idx], ex[idx],
            BFAM_REAL_ABS(field[idx]-ex[idx]) < BFAM_REAL_MIN);

      failures += fail;
    }
  }

  return failures;
}

static int
test_geo(MPI_Comm mpicomm)
{
  int failures = 0;
  int rank;
  BFAM_MPI_CHECK(MPI_Comm_rank(mpicomm, &rank));

  p4est_connectivity_t *conn = p4est_connectivity_new_star();

  for(p4est_topidx_t n = 0; n < conn->num_vertices; ++n)
    BFAM_INFO("%jd %25.15e %25.15e", (intmax_t) n+1,
       conn->vertices[n*3+0], conn->vertices[n*3+1]);

  bfam_domain_p4est_t* domain = bfam_domain_p4est_new(mpicomm, conn);

  p4est_balance(domain->p4est, P4EST_CONNECT_CORNER, NULL);
  p4est_partition(domain->p4est, NULL);
  // p4est_vtk_write_file(domain->p4est, NULL, "p4est_mesh");

  bfam_locidx_t numSubdomains = 1;
  int N = 1;

  bfam_locidx_t *subdomainID =
    bfam_calloc(domain->p4est->local_num_quadrants, sizeof(bfam_locidx_t));

  bfam_domain_p4est_split_dgx_quad_subdomains(domain, numSubdomains,
      subdomainID, &N);

  bfam_subdomain_t **subdomains =
    bfam_malloc(domain->d.numSubdomains*sizeof(bfam_subdomain_t**));

  const char* volume[] = {"_volume", NULL};
  bfam_domain_get_subdomains((bfam_domain_t*)domain, BFAM_DOMAIN_OR,
      volume, domain->d.numSubdomains, subdomains, &numSubdomains);

  BFAM_ABORT_IF_NOT(numSubdomains==1, "We should only have one subdomain");

  failures += check_field(subdomains[0], "_grid_x", xex);
  failures += check_field(subdomains[0], "_grid_y", yex);

  failures += check_field(subdomains[0], "_grid_rx", rxex);
  failures += check_field(subdomains[0], "_grid_ry", ryex);

  failures += check_field(subdomains[0], "_grid_sx", sxex);
  failures += check_field(subdomains[0], "_grid_sy", syex);

  failures += check_field(subdomains[0], "_grid_J", Jex);

  failures += check_field_face(subdomains[0], "_grid_nx", nxex);

  bfam_free(subdomainID);
  bfam_free(subdomains);

  bfam_domain_p4est_free(domain);
  bfam_free(domain);
  p4est_connectivity_destroy(conn);

  return failures;
}

int
main (int argc, char *argv[])
{
  int failures = 0;
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank;

  void *options= bfam_gopt_sort(&argc, (const char**)argv,
      bfam_gopt_start(
        bfam_gopt_option('h', 0,
                         bfam_gopt_shorts('h', '?'),
                         bfam_gopt_longs("help", "HELP")),
        bfam_gopt_option('V', 0,
                         bfam_gopt_shorts('V'),
                         bfam_gopt_longs("version")),
        bfam_gopt_option('v', BFAM_GOPT_REPEAT,
                         bfam_gopt_shorts('v'),
                         bfam_gopt_longs("verbose"))
        )
      );

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

  if(bfam_gopt(options, 'h'))
  {
    /*
     * if any of the help options was specified
     */
    BFAM_ROOT_INFO(helpText);
    exit(EXIT_SUCCESS);
  }

  if(bfam_gopt(options, 'V'))
  {
    BFAM_ROOT_INFO("BFAM Version: %s", bfam_version_get());
    BFAM_ROOT_INFO("BFAM Compile Info:\n" BFAM_COMPILE_INFO);
    exit( EXIT_SUCCESS );
  }

  int verbosity = bfam_gopt(options, 'v');

  BFAM_MPI_CHECK(MPI_Init(&argc,&argv));
  BFAM_MPI_CHECK(MPI_Comm_rank(comm, &rank));

  int logLevel = BFAM_MAX(BFAM_LL_INFO - verbosity, BFAM_LL_ALWAYS);

  bfam_log_init(rank, stdout, logLevel);
  bfam_signal_handler_set();

  int scLogPriorities = BFAM_MAX(SC_LP_STATISTICS - verbosity, SC_LP_ALWAYS);
  sc_init(comm, 0, 0, NULL, scLogPriorities);
  p4est_init(NULL, scLogPriorities);

  BFAM_MPI_CHECK(MPI_Barrier(comm));
  BFAM_ROOT_INFO("Test Geo\n");
  failures += test_geo(comm);


  sc_finalize();
  BFAM_MPI_CHECK(MPI_Finalize());

  bfam_gopt_free(options);

  return failures;
}

