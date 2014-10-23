#include <bfam.h>
#include <bfam_domain_pxest_2.h>

static int          refine_level = 0;

static int
refine_fn(p4est_t * p4est, p4est_topidx_t which_tree,
          p4est_quadrant_t * quadrant)
{
  if ((int) quadrant->level >= (refine_level - (int) (which_tree % 3)))
  {
    return 0;
  }
  if (quadrant->level == 1 && p4est_quadrant_child_id(quadrant) == 3)
  {
    return 1;
  }
  if (quadrant->x == P4EST_LAST_OFFSET (2) &&
      quadrant->y == P4EST_LAST_OFFSET (2))
  {
    return 1;
  }
  if (quadrant->x >= P4EST_QUADRANT_LEN (2))
  {
    return 0;
  }

  return 1;
}

static void
zero_field(bfam_locidx_t npoints, const char *name, bfam_real_t time,
    bfam_real_t *restrict x, bfam_real_t *restrict y, bfam_real_t *restrict z,
    struct bfam_subdomain *s, void *arg, bfam_real_t *restrict field)
{
  BFAM_ASSUME_ALIGNED(x, 32);
  BFAM_ASSUME_ALIGNED(y, 32);
  BFAM_ASSUME_ALIGNED(z, 32);
  BFAM_ASSUME_ALIGNED(field, 32);

  for(bfam_locidx_t n=0; n < npoints; ++n)
    field[n] = 0;
}

static void
x2_field(bfam_locidx_t npoints, const char *name, bfam_real_t time,
    bfam_real_t *restrict x, bfam_real_t *restrict y, bfam_real_t *restrict z,
    struct bfam_subdomain *s, void *arg, bfam_real_t *restrict field)
{
  BFAM_ASSUME_ALIGNED(x, 32);
  BFAM_ASSUME_ALIGNED(y, 32);
  BFAM_ASSUME_ALIGNED(z, 32);
  BFAM_ASSUME_ALIGNED(field, 32);

  for(bfam_locidx_t n=0; n < npoints; ++n)
    field[n] = x[n]*x[n];
}

static void
y2_field(bfam_locidx_t npoints, const char *name, bfam_real_t time,
    bfam_real_t *restrict x, bfam_real_t *restrict y, bfam_real_t *restrict z,
    struct bfam_subdomain *s, void *arg, bfam_real_t *restrict field)
{
  BFAM_ASSUME_ALIGNED(x, 32);
  BFAM_ASSUME_ALIGNED(y, 32);
  BFAM_ASSUME_ALIGNED(z, 32);
  BFAM_ASSUME_ALIGNED(field, 32);

  for(bfam_locidx_t n=0; n < npoints; ++n)
    field[n] = y[n]*y[n];
}

static void
z2_field(bfam_locidx_t npoints, const char *name, bfam_real_t time,
    bfam_real_t *restrict x, bfam_real_t *restrict y, bfam_real_t *restrict z,
    struct bfam_subdomain *s, void *arg, bfam_real_t *restrict field)
{
  BFAM_ASSUME_ALIGNED(x, 32);
  BFAM_ASSUME_ALIGNED(y, 32);
  BFAM_ASSUME_ALIGNED(z, 32);
  BFAM_ASSUME_ALIGNED(field, 32);

  for(bfam_locidx_t n=0; n < npoints; ++n)
    field[n] = x[n]*y[n];
}

static void
build_mesh(MPI_Comm mpicomm)
{
  int rank;
  BFAM_MPI_CHECK(MPI_Comm_rank(mpicomm, &rank));

  p4est_connectivity_t *conn = p4est_connectivity_new_corner();

  bfam_domain_pxest_t_2* domain = bfam_domain_pxest_new_2(mpicomm, conn);

  refine_level = 2;
  p4est_refine(domain->pxest, 1, refine_fn, NULL);
  p4est_balance(domain->pxest, P4EST_CONNECT_FACE, NULL);
  p4est_partition(domain->pxest, NULL);

  p4est_vtk_write_file(domain->pxest, NULL, "p4est_mesh");

  bfam_locidx_t numSubdomains = 4;
  bfam_locidx_t *subdomainID =
    bfam_malloc(domain->pxest->local_num_quadrants*sizeof(bfam_locidx_t));
  bfam_locidx_t *N = bfam_malloc(numSubdomains*sizeof(int));

  /*
   * Create an arbitrary splitting of the domain to test things.
   *
   * When use a subdomain id independent of MPI partition.  In practice
   * the subdomain id will be selected based on physics, element type, element
   * order, etc.
   *
   * For no particular reason increase element order with id
   */
  BFAM_ROOT_INFO("Splitting p4est into %jd DG Quad subdomains",
      (intmax_t) numSubdomains);
  for(bfam_locidx_t id = 0; id < numSubdomains; ++id)
  {
    N[id] = 3+id;

    p4est_gloidx_t first =
      p4est_partition_cut_gloidx(domain->pxest->global_num_quadrants,
          id, numSubdomains);

    p4est_gloidx_t last =
      p4est_partition_cut_gloidx(domain->pxest->global_num_quadrants,
          id + 1, numSubdomains) - 1;

    BFAM_ROOT_INFO("  id:%jd N:%d GIDs:%jd--%jd", (intmax_t) id, N[id],
        (intmax_t) first, (intmax_t) last);
  }

  p4est_gloidx_t gkOffset = domain->pxest->global_first_quadrant[rank];

  bfam_locidx_t idStart = 0;
  while(gkOffset >
      p4est_partition_cut_gloidx(domain->pxest->global_num_quadrants,
        idStart + 1, numSubdomains) - 1) ++idStart;

  for(p4est_locidx_t lk = 0, id = idStart;
      lk < domain->pxest->local_num_quadrants;
      ++lk)
  {
    p4est_gloidx_t gk = gkOffset + lk;

    if(gk > p4est_partition_cut_gloidx(domain->pxest->global_num_quadrants,
                                       id + 1, numSubdomains) - 1)
      ++id;

    BFAM_ASSERT(
      (gk >= p4est_partition_cut_gloidx(domain->pxest->global_num_quadrants,
                                   id, numSubdomains)) &&
      (gk < p4est_partition_cut_gloidx(domain->pxest->global_num_quadrants,
                                   id + 1, numSubdomains)));

    subdomainID[lk] = id;
  }

  bfam_domain_pxest_split_dgx_subdomains_2(domain, numSubdomains,
      subdomainID, N, NULL);

  const char *volume[] = {"_volume", NULL};
  bfam_vtk_write_file((bfam_domain_t*)domain, BFAM_DOMAIN_AND, volume,
                       NULL,"bfam_mesh",0, NULL, NULL, NULL, 1, 1,0);

  bfam_domain_add_field((bfam_domain_t*)domain, BFAM_DOMAIN_AND, volume, "v1");
  bfam_domain_add_field((bfam_domain_t*)domain, BFAM_DOMAIN_AND, volume, "v2");
  bfam_domain_add_field((bfam_domain_t*)domain, BFAM_DOMAIN_AND, volume, "v3");

  const char *glue[] = {"_glue", NULL};

  bfam_domain_add_minus_field((bfam_domain_t*)domain, BFAM_DOMAIN_AND, glue,
      "v1");
  bfam_domain_add_minus_field((bfam_domain_t*)domain, BFAM_DOMAIN_AND, glue,
      "v2");
  bfam_domain_add_minus_field((bfam_domain_t*)domain, BFAM_DOMAIN_AND, glue,
      "v3");

  bfam_domain_add_plus_field((bfam_domain_t*)domain, BFAM_DOMAIN_AND, glue,
      "v1");
  bfam_domain_add_plus_field((bfam_domain_t*)domain, BFAM_DOMAIN_AND, glue,
      "v2");
  bfam_domain_add_plus_field((bfam_domain_t*)domain, BFAM_DOMAIN_AND, glue,
      "v3");

  bfam_domain_init_field((bfam_domain_t*)domain, BFAM_DOMAIN_AND, volume, "v1",
      0, x2_field, NULL);
  bfam_domain_init_field((bfam_domain_t*)domain, BFAM_DOMAIN_AND, volume, "v2",
      0, y2_field, NULL);
  bfam_domain_init_field((bfam_domain_t*)domain, BFAM_DOMAIN_AND, volume, "v3",
      0, z2_field, NULL);

  const char *strain[] = {"E11", "E22", "E33", "E12", "E23", "E13", NULL};
  bfam_domain_add_fields((bfam_domain_t*)domain, BFAM_DOMAIN_AND, volume,
      strain);

  bfam_domain_add_minus_fields((bfam_domain_t*)domain, BFAM_DOMAIN_AND, glue,
      strain);
  bfam_domain_add_plus_fields( (bfam_domain_t*)domain, BFAM_DOMAIN_AND, glue,
      strain);

  bfam_domain_init_field((bfam_domain_t*)domain, BFAM_DOMAIN_AND, volume,
      "E11", 0, zero_field, NULL);
  bfam_domain_init_field((bfam_domain_t*)domain, BFAM_DOMAIN_AND, volume,
      "E22", 0, zero_field, NULL);
  bfam_domain_init_field((bfam_domain_t*)domain, BFAM_DOMAIN_AND, volume,
      "E33", 0, zero_field, NULL);
  bfam_domain_init_field((bfam_domain_t*)domain, BFAM_DOMAIN_AND, volume,
      "E12", 0, x2_field, NULL);
  bfam_domain_init_field((bfam_domain_t*)domain, BFAM_DOMAIN_AND, volume,
      "E23", 0, y2_field, NULL);
  bfam_domain_init_field((bfam_domain_t*)domain, BFAM_DOMAIN_AND, volume,
      "E13", 0, z2_field, NULL);

  const char *velocity[] = {"v", NULL};
  const char *velocityComp[] = {"v1", "v2", "v3", NULL};
  bfam_vtk_write_file((bfam_domain_t*)domain, BFAM_DOMAIN_AND, volume, NULL,
                      "bfam_fields",0, strain, velocity, velocityComp, 1, 1,0);

  const char* tags[] = {"_glue_parallel", "_glue_local", NULL};
  bfam_communicator_t* communicator =
    bfam_communicator_new((bfam_domain_t*)domain, BFAM_DOMAIN_AND, tags,
        mpicomm, 10, NULL);

  /* start recv_send */
  bfam_communicator_start(communicator);

  /* finish recv */
  bfam_communicator_finish(communicator);

  /* clean up */
  bfam_communicator_free(communicator);
  bfam_free(communicator);

  bfam_free(subdomainID);
  bfam_free(N);

  bfam_domain_pxest_free_2(domain);
  bfam_free(domain);
  p4est_connectivity_destroy(conn);
}

int
main (int argc, char *argv[])
{
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

  build_mesh(comm);

  sc_finalize();
  BFAM_MPI_CHECK(MPI_Finalize());

  bfam_gopt_free(options);

  return EXIT_SUCCESS;
}

