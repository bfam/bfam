#include <bfam.h>
#include "beard_dgx_rhs.h"
#include <p4est_iterate.h>

/*
 * Lua helper functions
 */
static int
lua_get_global_int(lua_State *L, const char *name, int def, int warning)
{
  lua_getglobal(L, name);
  int result = def;
  if(!lua_isnumber(L, -1))
  { /* for some reason gcc combines the two if without these. Why? */
    if(warning) BFAM_ROOT_WARNING("`%s' not found, using default %d",
        name, def);
  }
  else
    result = (int)lua_tonumber(L, -1);
  lua_pop(L, 1);
  return result;
}

typedef struct prefs
{
  lua_State *L;

  int min_refine_level;
  int max_refine_level;
} prefs_t;

/*
 * Set up the beard preference file
 */
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

  prefs->min_refine_level = lua_get_global_int(prefs->L,"min_refine_level",0,1);
  prefs->max_refine_level = lua_get_global_int(prefs->L,"max_refine_level",0,1);

  return prefs;
}

static void
free_prefs(prefs_t *prefs)
{
  lua_close(prefs->L);
}

static void
print_prefs(prefs_t *prefs)
{
  BFAM_ROOT_INFO("----------Preferences----------");
  BFAM_ROOT_INFO(" minimum refinement level = %d",prefs->min_refine_level);
  BFAM_ROOT_INFO(" maximum refinement level = %d",prefs->max_refine_level);
  BFAM_ROOT_INFO("-------------------------------");
}

/*
 * run the beard
 */
static void
run(MPI_Comm mpicomm, prefs_t *prefs)
{
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

  run(comm, prefs);
  free_prefs(prefs);
  bfam_free(prefs);

  sc_finalize();
  BFAM_MPI_CHECK(MPI_Finalize());

  bfam_gopt_free(options);

  return EXIT_SUCCESS;
}
