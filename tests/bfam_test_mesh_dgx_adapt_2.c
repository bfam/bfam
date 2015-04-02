#include <bfam.h>
#include <bfam_domain_pxest_2.h>

#define REAL_APPROX_EQ(x, y, K)                                                \
  BFAM_APPROX_EQ((x), (y), (K), BFAM_REAL_ABS, BFAM_REAL_EPS, (K)*BFAM_REAL_EPS)

static int refine_level = 0;

static int refine_fn(p4est_t *pxest, p4est_locidx_t which_tree,
                     p4est_quadrant_t *quadrant)
{
  if ((int)quadrant->level >= refine_level - (int)(1 - which_tree % 2))
    return 0;

  return 1;
}

static void poly0_field(bfam_locidx_t npoints, const char *name,
                        bfam_real_t time, bfam_real_t *restrict x,
                        bfam_real_t *restrict y, bfam_real_t *restrict z,
                        struct bfam_subdomain *s, void *arg,
                        bfam_real_t *restrict field)
{
  BFAM_ASSUME_ALIGNED(x, 32);
  BFAM_ASSUME_ALIGNED(y, 32);
  BFAM_ASSUME_ALIGNED(z, 32);
  BFAM_ASSUME_ALIGNED(field, 32);

  int *failures = (int *)arg;
  for (bfam_locidx_t n = 0; n < npoints; ++n)
  {
    const bfam_real_t val = 10;
    if (failures)
      (*failures) += !REAL_APPROX_EQ(field[n], val, 1000);
    else
      field[n] = val;
  }
}

static void poly1_field(bfam_locidx_t npoints, const char *name,
                        bfam_real_t time, bfam_real_t *restrict x,
                        bfam_real_t *restrict y, bfam_real_t *restrict z,
                        struct bfam_subdomain *s, void *arg,
                        bfam_real_t *restrict field)
{
  BFAM_ASSUME_ALIGNED(x, 32);
  BFAM_ASSUME_ALIGNED(y, 32);
  BFAM_ASSUME_ALIGNED(z, 32);
  BFAM_ASSUME_ALIGNED(field, 32);

  int *failures = (int *)arg;
  for (bfam_locidx_t n = 0; n < npoints; ++n)
  {
    const bfam_real_t val = x[n];
    if (failures)
      (*failures) += !REAL_APPROX_EQ(field[n], val, 1000);
    else
      field[n] = val;
  }
}

static void poly2_field(bfam_locidx_t npoints, const char *name,
                        bfam_real_t time, bfam_real_t *restrict x,
                        bfam_real_t *restrict y, bfam_real_t *restrict z,
                        struct bfam_subdomain *s, void *arg,
                        bfam_real_t *restrict field)
{
  BFAM_ASSUME_ALIGNED(x, 32);
  BFAM_ASSUME_ALIGNED(y, 32);
  BFAM_ASSUME_ALIGNED(z, 32);
  BFAM_ASSUME_ALIGNED(field, 32);

  int *failures = (int *)arg;
  for (bfam_locidx_t n = 0; n < npoints; ++n)
  {
    const bfam_real_t val = y[n];
    if (failures)
      (*failures) += !REAL_APPROX_EQ(field[n], val, 1000);
    else
      field[n] = val;
  }
}

static void poly3_field(bfam_locidx_t npoints, const char *name,
                        bfam_real_t time, bfam_real_t *restrict x,
                        bfam_real_t *restrict y, bfam_real_t *restrict z,
                        struct bfam_subdomain *s, void *arg,
                        bfam_real_t *restrict field)
{
  BFAM_ASSUME_ALIGNED(x, 32);
  BFAM_ASSUME_ALIGNED(y, 32);
  BFAM_ASSUME_ALIGNED(field, 32);

  int *failures = (int *)arg;
  for (bfam_locidx_t n = 0; n < npoints; ++n)
  {
    const bfam_real_t val = x[n] * y[n];
    if (failures)
      (*failures) += !REAL_APPROX_EQ(field[n], val, 1000);
    else
      field[n] = val;
  }
}

static void poly4_field(bfam_locidx_t npoints, const char *name,
                        bfam_real_t time, bfam_real_t *restrict x,
                        bfam_real_t *restrict y, bfam_real_t *restrict z,
                        struct bfam_subdomain *s, void *arg,
                        bfam_real_t *restrict field)
{
  BFAM_ASSUME_ALIGNED(x, 32);
  BFAM_ASSUME_ALIGNED(y, 32);
  BFAM_ASSUME_ALIGNED(z, 32);
  BFAM_ASSUME_ALIGNED(field, 32);

  int *failures = (int *)arg;
  for (bfam_locidx_t n = 0; n < npoints; ++n)
  {
    const bfam_real_t X = x[n] / 2.5;
    const bfam_real_t Y = y[n] / 2.5;
    const bfam_real_t Z = 1.0;
    const bfam_real_t val = X * Y * Z * Z + X * X * Y + Z * Z * Z;
    if (failures)
      (*failures) += !REAL_APPROX_EQ(field[n], val, 1000);
    else
      field[n] = val;
  }
}

static void mark_elements(bfam_domain_pxest_t_2 *domain)
{
  bfam_domain_t *dbase = &domain->base;
  bfam_subdomain_t **subdomains =
      bfam_malloc(dbase->numSubdomains * sizeof(bfam_subdomain_t *));

  bfam_locidx_t num_subdomains = 0;

  const char *volume[] = {"_volume", NULL};

  bfam_domain_get_subdomains(dbase, BFAM_DOMAIN_AND, volume,
                             dbase->numSubdomains, subdomains, &num_subdomains);

  for (bfam_locidx_t s = 0; s < num_subdomains / 2; ++s)
  {
    bfam_subdomain_dgx_t *sub = (bfam_subdomain_dgx_t *)subdomains[s];
    for (bfam_locidx_t k = 0; k < sub->K; ++k)
    {
      sub->padapt[k] = 8;
      sub->hadapt[k] = BFAM_FLAG_COARSEN;
    }
  }

  for (bfam_locidx_t s = num_subdomains / 2; s < num_subdomains; ++s)
  {
    bfam_subdomain_dgx_t *sub = (bfam_subdomain_dgx_t *)subdomains[s];
    for (bfam_locidx_t k = 0; k < sub->K; ++k)
    {
      sub->padapt[k] = 5;
      sub->hadapt[k] = BFAM_FLAG_REFINE;
    }
  }

  bfam_free(subdomains);
}

static int build_mesh(MPI_Comm mpicomm)
{
  int failures = 0;
  int rank;
  BFAM_MPI_CHECK(MPI_Comm_rank(mpicomm, &rank));

  p4est_connectivity_t *conn = p4est_connectivity_new_disk();

  bfam_domain_pxest_t_2 *domain = bfam_domain_pxest_new_2(mpicomm, conn);

  refine_level = 1;
  p4est_refine(domain->pxest, 2, refine_fn, bfam_domain_pxest_init_callback_2);
  p4est_balance(domain->pxest, P4EST_CONNECT_CORNER,
                bfam_domain_pxest_init_callback_2);
  p4est_partition(domain->pxest, 1, NULL);

  bfam_locidx_t numSubdomains = 11;
  bfam_locidx_t *subdomainID =
      bfam_malloc(domain->pxest->local_num_quadrants * sizeof(bfam_locidx_t));
  bfam_locidx_t *N = bfam_malloc(numSubdomains * sizeof(int));

  /*
   * Create an arbitrary splitting of the domain to test things.
   *
   * We use a subdomain id independent of MPI partition.  In practice the
   * subdomain id will be selected based on physics, element type, element
   * order, etc.
   *
   * For no particular reason increase element order with id
   */
  BFAM_ROOT_INFO("Splitting pxest into %jd DG Quad subdomains",
                 (intmax_t)numSubdomains);
  for (bfam_locidx_t id = 0; id < numSubdomains; ++id)
  {
    N[id] = 5 + id;

    p4est_gloidx_t first = p4est_partition_cut_gloidx(
        domain->pxest->global_num_quadrants, id, numSubdomains);

    p4est_gloidx_t last =
        p4est_partition_cut_gloidx(domain->pxest->global_num_quadrants, id + 1,
                                   numSubdomains) -
        1;

    BFAM_ROOT_INFO("  id:%jd N:%d GIDs:%jd--%jd", (intmax_t)id, N[id],
                   (intmax_t)first, (intmax_t)last);
  }

  p4est_gloidx_t gkOffset = domain->pxest->global_first_quadrant[rank];

  bfam_locidx_t idStart = 0;
  while (gkOffset >
         p4est_partition_cut_gloidx(domain->pxest->global_num_quadrants,
                                    idStart + 1, numSubdomains) -
             1)
    ++idStart;

  for (p4est_locidx_t lk = 0, id = idStart;
       lk < domain->pxest->local_num_quadrants; ++lk)
  {
    p4est_gloidx_t gk = gkOffset + lk;

    if (gk > p4est_partition_cut_gloidx(domain->pxest->global_num_quadrants,
                                        id + 1, numSubdomains) -
                 1)
      ++id;

    BFAM_ASSERT(
        (gk >= p4est_partition_cut_gloidx(domain->pxest->global_num_quadrants,
                                          id, numSubdomains)) &&
        (gk < p4est_partition_cut_gloidx(domain->pxest->global_num_quadrants,
                                         id + 1, numSubdomains)));

    subdomainID[lk] = id;
  }

  bfam_domain_pxest_split_dgx_subdomains_2(domain, numSubdomains, subdomainID,
                                           NULL, N, NULL, NULL, NULL, NULL,
                                           NULL);

  const char *volume[] = {"_volume", NULL};

  bfam_domain_add_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, volume, "p0");
  bfam_domain_add_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, volume, "p1");
  bfam_domain_add_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, volume, "p2");
  bfam_domain_add_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, volume, "p3");
  bfam_domain_add_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, volume, "p4");

  bfam_domain_init_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, volume, "p0",
                         0, poly0_field, NULL);
  bfam_domain_init_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, volume, "p1",
                         0, poly1_field, NULL);
  bfam_domain_init_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, volume, "p2",
                         0, poly2_field, NULL);
  bfam_domain_init_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, volume, "p3",
                         0, poly3_field, NULL);
  bfam_domain_init_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, volume, "p4",
                         0, poly4_field, NULL);

  const char *ps[] = {"p1", "p2", "p3", NULL};

  p4est_vtk_write_file(domain->pxest, NULL, "p4est_mesh_pre");
  bfam_vtk_write_file((bfam_domain_t *)domain, BFAM_DOMAIN_OR, volume, NULL,
                      "dgx_adapt_2_pre", 0, ps, NULL, NULL, 0, 0, 0);

  mark_elements(domain);
  bfam_domain_pxest_adapt_2(domain, NULL, NULL, NULL, NULL);

  p4est_vtk_write_file(domain->pxest, NULL, "p4est_mesh_post");
  bfam_vtk_write_file((bfam_domain_t *)domain, BFAM_DOMAIN_OR, volume, NULL,
                      "dgx_adapt_2_post", 0, ps, NULL, NULL, 0, 0, 0);

  bfam_domain_init_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, volume, "p0",
                         0, poly0_field, &failures);
  bfam_domain_init_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, volume, "p1",
                         0, poly1_field, &failures);
  bfam_domain_init_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, volume, "p2",
                         0, poly2_field, &failures);
  bfam_domain_init_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, volume, "p3",
                         0, poly3_field, &failures);
  bfam_domain_init_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, volume, "p4",
                         0, poly4_field, &failures);

  bfam_free(subdomainID);
  bfam_free(N);

  bfam_domain_pxest_free_2(domain);
  bfam_free(domain);
  p4est_connectivity_destroy(conn);

  return failures;
}

static int run_tests(MPI_Comm mpicomm)
{
  int failures = 0;

  failures += build_mesh(mpicomm);

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
      "  there are four possible options to this program, some of which \n"
      "  have multiple names:\n"
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

  failures += run_tests(comm);

  sc_finalize();
  BFAM_MPI_CHECK(MPI_Finalize());

  bfam_gopt_free(options);

  if (failures > 0)
    BFAM_WARNING("FAIL! with failures %d", failures);
  return failures;
}
