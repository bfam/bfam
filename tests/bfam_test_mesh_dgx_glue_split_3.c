#include <bfam.h>
#include <bfam_domain_pxest_3.h>

#define REAL_APPROX_EQ(x, y, K)                                              \
  BFAM_APPROX_EQ((x), (y), (K), BFAM_REAL_ABS, BFAM_REAL_EPS, (K)*BFAM_REAL_EPS)

static int          refine_level = 0;

static int
refine_fn(p8est_t* pxest, p4est_locidx_t which_tree, p8est_quadrant_t* quadrant)
{
  if ((int)quadrant->level >= refine_level - (int)(1 - which_tree % 2))
    return 0;

  return 1;
}

typedef struct
{
  p8est_connectivity_t *conn;
  bfam_domain_pxest_t_3 *domain;
} state_t;

static int
build_state(MPI_Comm mpicomm, state_t* state)
{
  int failures = 0;
  int rank;
  BFAM_MPI_CHECK(MPI_Comm_rank(mpicomm, &rank));

  state->conn = p8est_connectivity_new_twocubes();
  const bfam_locidx_t tree_to_glueid[2 * 6] = {
    -1,  0, -1, -1, -1, -1,
     0, -1, -1, -1, -1, -1,
  };
  state->domain = bfam_domain_pxest_new_3(mpicomm, state->conn);


  bfam_domain_pxest_t_3 *domain = state->domain;

  refine_level = 1;
  p8est_refine(domain->pxest, 2, refine_fn, NULL);
  p8est_balance(domain->pxest, P8EST_CONNECT_CORNER, NULL);
  p8est_partition(domain->pxest, NULL);

  p8est_vtk_write_file(domain->pxest, NULL, "p8est_mesh");

  bfam_locidx_t numSubdomains = 1;
  bfam_locidx_t *subdomainID =
    bfam_malloc(domain->pxest->local_num_quadrants*sizeof(bfam_locidx_t));
  bfam_locidx_t *N = bfam_malloc(numSubdomains*sizeof(int));

  bfam_locidx_t *glueID =
    bfam_malloc(domain->pxest->local_num_quadrants*P4EST_FACES*
        sizeof(bfam_locidx_t));

  /*
   * Create an arbitrary splitting of the domain to test things.
   *
   * When use a subdomain id independent of MPI partition.  In practice
   * the subdomain id will be selectmd based on physics, element type, element
   * order, etc.
   *
   * For no particular reason increase element order with id
   */
  BFAM_ROOT_INFO("Splitting pxest into %jd DG Quad subdomains",
      (intmax_t) numSubdomains);
  for(bfam_locidx_t id = 0; id < numSubdomains; ++id)
  {
    N[id] = 1 + id;

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

  bfam_domain_pxest_quad_to_glueid_3(domain->pxest, tree_to_glueid,
      glueID);

  for(p4est_locidx_t k = 0; k < domain->pxest->local_num_quadrants; ++k)
  {
    for(int f = 0; f < P4EST_FACES; ++f)
    {
      BFAM_LDEBUG("glueID[%3d][%d] = %jd", k, f,
          (intmax_t)glueID[P4EST_FACES*k+f]);
    }
  }


  bfam_domain_pxest_split_dgx_subdomains_3(domain, numSubdomains,
      subdomainID, N, glueID);

  const char *volume[] = {"_volume", NULL};
  const char *glue[]   = {"_glue_parallel", "_glue_local", NULL};

  bfam_domain_t *d = (bfam_domain_t*)domain;

  bfam_domain_add_field(d, BFAM_DOMAIN_OR, glue, "_grid_x0");
  bfam_domain_add_field(d, BFAM_DOMAIN_OR, glue, "_grid_x1");
  bfam_domain_add_field(d, BFAM_DOMAIN_OR, glue, "_grid_x2");

  bfam_domain_add_plus_field(d, BFAM_DOMAIN_OR, glue, "_grid_x0");
  bfam_domain_add_plus_field(d, BFAM_DOMAIN_OR, glue, "_grid_x1");
  bfam_domain_add_plus_field(d, BFAM_DOMAIN_OR, glue, "_grid_x2");

  bfam_domain_add_minus_field(d, BFAM_DOMAIN_OR, glue, "_grid_x0");
  bfam_domain_add_minus_field(d, BFAM_DOMAIN_OR, glue, "_grid_x1");
  bfam_domain_add_minus_field(d, BFAM_DOMAIN_OR, glue, "_grid_x2");

  bfam_subdomain_comm_args_t commargs;
  const char *comm_args_face_scalars[]      = {NULL};
  const char *comm_args_scalars[]           = {"_grid_x0",
                                               "_grid_x1",
                                               "_grid_x2",
                                               NULL};
  const char *comm_args_vectors[]           = {NULL};
  const char *comm_args_vector_components[] = {NULL};
  const char *comm_args_tensors[]           = {NULL};
  const char *comm_args_tensor_components[] = {NULL};

  commargs.scalars_m           = comm_args_scalars;
  commargs.scalars_p           = comm_args_scalars;

  commargs.vectors_m           = comm_args_vectors;
  commargs.vectors_p           = comm_args_vectors;
  commargs.vector_components_m = comm_args_vector_components;
  commargs.vector_components_p = comm_args_vector_components;

  commargs.tensors_m           = comm_args_tensors;
  commargs.tensors_p           = comm_args_tensors;
  commargs.tensor_components_m = comm_args_tensor_components;
  commargs.tensor_components_p = comm_args_tensor_components;

  commargs.face_scalars_m = comm_args_face_scalars;
  commargs.face_scalars_p = comm_args_face_scalars;

  bfam_communicator_t* communicator =
    bfam_communicator_new(d, BFAM_DOMAIN_OR, glue, mpicomm, 11, &commargs);
  bfam_communicator_start(communicator);
  bfam_communicator_finish(communicator);
  bfam_communicator_free(communicator);
  bfam_free(communicator);

  {
    bfam_locidx_t numElements = d->numSubdomains;
    bfam_subdomain_t **subdomains =
      bfam_malloc(numElements*sizeof(bfam_subdomain_t*));
    bfam_locidx_t numSubdomains;

    bfam_domain_get_subdomains(d, BFAM_DOMAIN_OR, glue,
        numElements, subdomains, &numSubdomains);

    for(bfam_locidx_t i = 0; i < numSubdomains; ++i)
    {
      bfam_subdomain_t *s = subdomains[i];
      bfam_subdomain_dgx_t *sub = (bfam_subdomain_dgx_t*) s;

      bfam_real_t *restrict x0 =
        bfam_dictionary_get_value_ptr(&s->fields, "_grid_x0");
      bfam_real_t *restrict x1 =
        bfam_dictionary_get_value_ptr(&s->fields, "_grid_x1");
      bfam_real_t *restrict x2 =
        bfam_dictionary_get_value_ptr(&s->fields, "_grid_x2");

      bfam_real_t *restrict x0_m =
        bfam_dictionary_get_value_ptr(&s->glue_m->fields, "_grid_x0");
      bfam_real_t *restrict x1_m =
        bfam_dictionary_get_value_ptr(&s->glue_m->fields, "_grid_x1");
      bfam_real_t *restrict x2_m =
        bfam_dictionary_get_value_ptr(&s->glue_m->fields, "_grid_x2");

      bfam_real_t *restrict x0_p =
        bfam_dictionary_get_value_ptr(&s->glue_p->fields, "_grid_x0");
      bfam_real_t *restrict x1_p =
        bfam_dictionary_get_value_ptr(&s->glue_p->fields, "_grid_x1");
      bfam_real_t *restrict x2_p =
        bfam_dictionary_get_value_ptr(&s->glue_p->fields, "_grid_x2");

      BFAM_ASSUME_ALIGNED(x0, 32);
      BFAM_ASSUME_ALIGNED(x1, 32);
      BFAM_ASSUME_ALIGNED(x2, 32);

      BFAM_ASSUME_ALIGNED(x0_m, 32);
      BFAM_ASSUME_ALIGNED(x1_m, 32);
      BFAM_ASSUME_ALIGNED(x2_m, 32);

      BFAM_ASSUME_ALIGNED(x0_p, 32);
      BFAM_ASSUME_ALIGNED(x1_p, 32);
      BFAM_ASSUME_ALIGNED(x2_p, 32);

      for(bfam_locidx_t j = 0; j < sub->Np * sub->K; ++j)
      {
        x0[j] = (x0_m[j] + x0_p[j])/2;
        x1[j] = (x1_m[j] + x1_p[j])/2;
        x2[j] = (x2_m[j] + x2_p[j])/2;
      }
    }

    bfam_free(subdomains);
  }

  bfam_vtk_write_file(d, BFAM_DOMAIN_OR, volume,
                      "","volume",0, NULL, NULL, NULL, 0, 0, 0);

  const char *glue_id_0[]   = {"_glue_id_0", NULL};

  bfam_vtk_write_file(d, BFAM_DOMAIN_OR, glue_id_0,
                      "","glue_id",0, NULL, NULL, NULL, 0, 0, 0);

  bfam_free(subdomainID);
  bfam_free(glueID);
  bfam_free(N);

  return failures;
}

static int
free_state(MPI_Comm mpicomm, state_t* state)
{
  int failures = 0;

  bfam_domain_pxest_free_3(state->domain);
  bfam_free(state->domain);
  p8est_connectivity_destroy(state->conn);

  return failures;
}

static int
run_tests(MPI_Comm mpicomm)
{
  int failures = 0;

  state_t astate; state_t *state = &astate;

  failures += build_state(mpicomm, state);

  failures += free_state(mpicomm, state);

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

  failures += run_tests(comm);

  sc_finalize();
  BFAM_MPI_CHECK(MPI_Finalize());

  bfam_gopt_free(options);

  if(failures > 0) BFAM_WARNING("FAIL! with failures %d",failures);
  return failures;
}

