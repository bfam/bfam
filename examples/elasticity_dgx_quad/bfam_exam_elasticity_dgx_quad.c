#include <bfam.h>
#include "bfam_exam_elasticity_dgx_quad_rhs.h"

static int refine_level = 0;

/*
 * Uniform refinement function
 */
static int
refine_fn(p4est_t * p4est, p4est_topidx_t which_tree,
          p4est_quadrant_t * quadrant)
{
  if((int)quadrant->level >= refine_level)
  {
    return 0;
  }
  else
  {
    return 1;
  }
}

#define REAL_APPROX_EQ(x, y, K)                                              \
  BFAM_APPROX_EQ((x), (y), (K), BFAM_REAL_ABS, BFAM_REAL_EPS, BFAM_REAL_EPS)

typedef struct prefs
{
  lua_State *L;

  int N;
  bfam_locidx_t num_subdomains;

  p4est_connectivity_t * (*conn_fn) (void);
} prefs_t;

typedef struct exam
{
  MPI_Comm mpicomm;
  int      mpirank;
  int      mpisize;

  p4est_connectivity_t *conn;
  bfam_domain_p4est_t  *domain;
} exam_t;

static void
init_mpi(exam_t *exam, MPI_Comm mpicomm)
{
  exam->mpicomm = mpicomm;
  BFAM_MPI_CHECK(MPI_Comm_rank(mpicomm, &exam->mpirank));
  BFAM_MPI_CHECK(MPI_Comm_size(mpicomm, &exam->mpisize));
}

static void
split_domain_arbitrary(exam_t *exam, int base_N, bfam_locidx_t num_subdomains)
{
  bfam_domain_p4est_t *domain = exam->domain;
  bfam_locidx_t *subdomain_id =
    bfam_malloc(domain->p4est->local_num_quadrants*sizeof(bfam_locidx_t));
  bfam_locidx_t *N = bfam_malloc(num_subdomains*sizeof(int));

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
      (intmax_t) num_subdomains);
  for(bfam_locidx_t id = 0; id < num_subdomains; ++id)
  {
    N[id] = base_N+id;

    p4est_gloidx_t first =
      p4est_partition_cut_gloidx(domain->p4est->global_num_quadrants,
          id, num_subdomains);

    p4est_gloidx_t last =
      p4est_partition_cut_gloidx(domain->p4est->global_num_quadrants,
          id + 1, num_subdomains) - 1;

    BFAM_ROOT_INFO("  id:%jd N:%d GIDs:%jd--%jd", (intmax_t) id, N[id],
        (intmax_t) first, (intmax_t) last);
  }

  p4est_gloidx_t gk_offset =
    domain->p4est->global_first_quadrant[exam->mpirank];

  bfam_locidx_t id_start = 0;
  while(gk_offset >
      p4est_partition_cut_gloidx(domain->p4est->global_num_quadrants,
        id_start + 1, num_subdomains) - 1) ++id_start;

  for(p4est_locidx_t lk = 0, id = id_start;
      lk < domain->p4est->local_num_quadrants;
      ++lk)
  {
    p4est_gloidx_t gk = gk_offset + lk;

    if(gk > p4est_partition_cut_gloidx(domain->p4est->global_num_quadrants,
                                       id + 1, num_subdomains) - 1)
      ++id;

    BFAM_ASSERT(
      (gk >= p4est_partition_cut_gloidx(domain->p4est->global_num_quadrants,
                                   id, num_subdomains)) &&
      (gk < p4est_partition_cut_gloidx(domain->p4est->global_num_quadrants,
                                   id + 1, num_subdomains)));

    subdomain_id[lk] = id;
  }

  bfam_domain_p4est_split_dgx_quad_subdomains(domain, num_subdomains,
      subdomain_id, N);

  bfam_free(subdomain_id);
  bfam_free(N);
}

static void
zero_field(bfam_locidx_t npoints, bfam_real_t time, bfam_real_t *restrict x,
    bfam_real_t *restrict y, bfam_real_t *restrict z, struct bfam_subdomain *s,
    void *arg, bfam_real_t *restrict field)
{
  BFAM_ASSUME_ALIGNED(x, 32);
  BFAM_ASSUME_ALIGNED(y, 32);
  BFAM_ASSUME_ALIGNED(z, 32);
  BFAM_ASSUME_ALIGNED(field, 32);

  for(bfam_locidx_t n=0; n < npoints; ++n)
    field[n] = 0;
}



static void
init_domain(exam_t *exam, prefs_t *prefs)
{
  exam->conn = prefs->conn_fn();

  exam->domain = bfam_domain_p4est_new(exam->mpicomm, exam->conn);

  p4est_refine(exam->domain->p4est, 1, refine_fn, NULL);
  p4est_balance(exam->domain->p4est, P4EST_CONNECT_CORNER, NULL);
  p4est_partition(exam->domain->p4est, NULL);

  p4est_vtk_write_file(exam->domain->p4est, NULL, "p4est_mesh");

  split_domain_arbitrary(exam, prefs->N, prefs->num_subdomains);

  const char *volume[] = {"_volume", NULL};
  const char *glue[]   = {"_glue_parallel", "_glue_local", NULL};

  bfam_domain_t *domain = (bfam_domain_t*)exam->domain;

  bfam_domain_add_field(domain, BFAM_DOMAIN_OR, volume, "rho");
  bfam_domain_add_field(domain, BFAM_DOMAIN_OR, volume, "lambda");
  bfam_domain_add_field(domain, BFAM_DOMAIN_OR, volume, "mu");

  bfam_domain_add_field(domain, BFAM_DOMAIN_OR, volume, "v1");
  bfam_domain_add_field(domain, BFAM_DOMAIN_OR, volume, "v2");
  bfam_domain_add_field(domain, BFAM_DOMAIN_OR, volume, "v3");

  bfam_domain_add_field(domain, BFAM_DOMAIN_OR, volume, "S11");
  bfam_domain_add_field(domain, BFAM_DOMAIN_OR, volume, "S22");
  bfam_domain_add_field(domain, BFAM_DOMAIN_OR, volume, "S33");
  bfam_domain_add_field(domain, BFAM_DOMAIN_OR, volume, "S12");
  bfam_domain_add_field(domain, BFAM_DOMAIN_OR, volume, "S13");
  bfam_domain_add_field(domain, BFAM_DOMAIN_OR, volume, "S23");

  bfam_domain_add_minus_field(domain, BFAM_DOMAIN_OR, glue, "v1");
  bfam_domain_add_minus_field(domain, BFAM_DOMAIN_OR, glue, "v2");
  bfam_domain_add_minus_field(domain, BFAM_DOMAIN_OR, glue, "v3");

  bfam_domain_add_plus_field( domain, BFAM_DOMAIN_OR, glue, "v1");
  bfam_domain_add_plus_field( domain, BFAM_DOMAIN_OR, glue, "v2");
  bfam_domain_add_plus_field( domain, BFAM_DOMAIN_OR, glue, "v3");

  bfam_domain_add_minus_field(domain, BFAM_DOMAIN_OR, glue, "S11");
  bfam_domain_add_minus_field(domain, BFAM_DOMAIN_OR, glue, "S22");
  bfam_domain_add_minus_field(domain, BFAM_DOMAIN_OR, glue, "S33");
  bfam_domain_add_minus_field(domain, BFAM_DOMAIN_OR, glue, "S12");
  bfam_domain_add_minus_field(domain, BFAM_DOMAIN_OR, glue, "S13");
  bfam_domain_add_minus_field(domain, BFAM_DOMAIN_OR, glue, "S23");

  bfam_domain_add_plus_field( domain, BFAM_DOMAIN_OR, glue, "S11");
  bfam_domain_add_plus_field( domain, BFAM_DOMAIN_OR, glue, "S22");
  bfam_domain_add_plus_field( domain, BFAM_DOMAIN_OR, glue, "S33");
  bfam_domain_add_plus_field( domain, BFAM_DOMAIN_OR, glue, "S12");
  bfam_domain_add_plus_field( domain, BFAM_DOMAIN_OR, glue, "S13");
  bfam_domain_add_plus_field( domain, BFAM_DOMAIN_OR, glue, "S23");

  bfam_domain_init_field(domain, BFAM_DOMAIN_OR, volume, "v1", 0, zero_field,
      NULL);
  bfam_domain_init_field(domain, BFAM_DOMAIN_OR, volume, "v2", 0, zero_field,
      NULL);
  bfam_domain_init_field(domain, BFAM_DOMAIN_OR, volume, "v3", 0, zero_field,
      NULL);
  bfam_domain_init_field(domain, BFAM_DOMAIN_OR, volume, "S11", 0, zero_field,
      NULL);
  bfam_domain_init_field(domain, BFAM_DOMAIN_OR, volume, "S22", 0, zero_field,
      NULL);
  bfam_domain_init_field(domain, BFAM_DOMAIN_OR, volume, "S33", 0, zero_field,
      NULL);
  bfam_domain_init_field(domain, BFAM_DOMAIN_OR, volume, "S12", 0, zero_field,
      NULL);
  bfam_domain_init_field(domain, BFAM_DOMAIN_OR, volume, "S13", 0, zero_field,
      NULL);
  bfam_domain_init_field(domain, BFAM_DOMAIN_OR, volume, "S23", 0, zero_field,
      NULL);
}

static void
free_exam(exam_t *exam)
{
  bfam_domain_p4est_free(exam->domain);
  bfam_free(exam->domain);
  p4est_connectivity_destroy(exam->conn);
}

static void
run(MPI_Comm mpicomm, prefs_t *prefs)
{
  exam_t exam;

  init_mpi(&exam, mpicomm);

  init_domain(&exam, prefs);

  free_exam(&exam);
}

static void
print_order(int N)
{
#define X(order) \
  case order: bfam_elasticity_dgx_quad_print_order_##order(N); break;

  switch(N)
  {
    BFAM_LIST_OF_DGX_QUAD_NORDERS
    default:
      bfam_elasticity_dgx_quad_print_order_(N);
      break;
  }
#undef X
}

struct conn_table {
  const char *name;
  p4est_connectivity_t * (*conn_fn) (void);
} conn_table[] = {
  {"unitsquare", &p4est_connectivity_new_unitsquare},
  {"periodic",   &p4est_connectivity_new_periodic},
  {"rotwrap",    &p4est_connectivity_new_rotwrap},
  {"corner",     &p4est_connectivity_new_corner},
  {"moebius",    &p4est_connectivity_new_moebius},
  {"star",       &p4est_connectivity_new_star},
  {NULL,         NULL}
};

static int
get_global_int(lua_State *L, const char *name, int def)
{
  lua_getglobal(L, name);
  int result = def;
  if(!lua_isnumber(L, -1))
    BFAM_ROOT_WARNING("`%s' not found, using default", name);
  else
    result = (int)lua_tonumber(L, -1);
  lua_pop(L, 1);
  return result;
}

static prefs_t *
new_prefs(const char *prefs_filename)
{
  prefs_t *prefs = bfam_malloc(sizeof(prefs_t));

  lua_State *L = luaL_newstate();
  luaL_openlibs(L);

  if(luaL_loadfile(L, prefs_filename) || lua_pcall(L, 0, 0, 0))
    BFAM_LERROR("cannot run configuration file: `%s'", lua_tostring(L, -1));

  prefs->L = L;

  BFAM_ASSERT(lua_gettop(L)==0);

  prefs->N = get_global_int(L, "N", 5);
  prefs->num_subdomains = get_global_int(L, "num_subdomains", 1);
  BFAM_ASSERT(lua_gettop(L)==0);

  lua_getglobal(L, "connectivity");
  prefs->conn_fn = conn_table[0].conn_fn;
  if(lua_isstring(L, -1))
  {
    int i;
    const char *conn_name = lua_tostring(L, -1);
    for(i = 0; conn_table[i].name != NULL; ++i)
    {
      if(strcmp(conn_name, conn_table[i].name) == 0)
        break;
    }

    if(conn_table[i].name == NULL)
      BFAM_LERROR("invalid connectivity name: `%s'; using default", conn_name);
    else
    {
      prefs->conn_fn = conn_table[i].conn_fn;
    }
  }
  else
  {
    BFAM_ROOT_WARNING("`connectivity' not found, using default");
  }
  BFAM_ASSERT(prefs->conn_fn != NULL);
  lua_pop(L, 1);

  BFAM_ASSERT(lua_gettop(L)==0);

  return prefs;
}

static void
print_prefs(prefs_t *prefs)
{
  BFAM_ROOT_INFO("----------Preferences----------");
  BFAM_ROOT_INFO("N=%d", prefs->N);
  BFAM_ROOT_INFO("num_subdomains=%"BFAM_LOCIDX_PRId, prefs->num_subdomains);
  for(int i=0; conn_table[i].name!=NULL; ++i)
    if(conn_table[i].conn_fn == prefs->conn_fn)
      BFAM_ROOT_INFO("conn_fn=`%s'", conn_table[i].name);
  BFAM_ROOT_INFO("-------------------------------");
}

static void
free_prefs(prefs_t *prefs)
{
  lua_close(prefs->L);
}

int
main(int argc, char *argv[])
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

  const char *help_text =
  "  %s [options] prefs_file\n"
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
    BFAM_ROOT_INFO(help_text, argv[0]);
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

  if(argc != 2)
  {
    BFAM_LERROR("Unexpected number of arguments.");
    BFAM_ROOT_INFO(help_text, argv[0]);
    exit(EXIT_FAILURE);
  }

  int logLevel = BFAM_MAX(BFAM_LL_INFO - verbosity, BFAM_LL_ALWAYS);

  bfam_log_init(rank, stdout, logLevel);
  bfam_signal_handler_set();

  int sc_log_priorities = BFAM_MAX(SC_LP_STATISTICS - verbosity, SC_LP_ALWAYS);
  sc_init(comm, 0, 0, NULL, sc_log_priorities);
  p4est_init(NULL, sc_log_priorities);

  prefs_t *prefs = new_prefs(argv[1]);
  print_prefs(prefs);

  print_order(prefs->N);

  run(comm, prefs);
  free_prefs(prefs);
  bfam_free(prefs);

  sc_finalize();
  BFAM_MPI_CHECK(MPI_Finalize());

  bfam_gopt_free(options);

  return EXIT_SUCCESS;
}
