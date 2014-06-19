#include <bfam.h>

#ifdef BLADE_DGX_DIMENSION

/* Handle the header files */
#if   BLADE_DGX_DIMENSION==2
#include <bfam_domain_pxest_2.h>
// #include "blade_dgx_rhs_2.h"
#include <p4est_iterate.h>
#elif BLADE_DGX_DIMENSION==3
#include <bfam_domain_pxest_3.h>
// #include "blade_dgx_rhs_3.h"
#include <p8est_iterate.h>
#define P4EST_CONNECT_CORNER P8EST_CONNECT_CORNER
#else
#error "bad dimension"
#endif

/* Some convience macros */
#if   BLADE_DGX_DIMENSION==2

#define DIM  2
#define BLADE_D3_AP(A1,A2) (A1)
#define BLADE_D3_OP(A) BFAM_NOOP()
#define bfam_domain_pxest_new bfam_domain_pxest_new_2
#define bfam_domain_pxest_quad_to_glueid \
  bfam_domain_pxest_quad_to_glueid_2
#define bfam_domain_pxest_split_dgx_subdomains \
  bfam_domain_pxest_split_dgx_subdomains_2
#define bfam_domain_pxest_free \
  bfam_domain_pxest_free_2
#define bfam_domain_pxest_t  bfam_domain_pxest_t_2

#elif BLADE_DGX_DIMENSION==3

#define DIM  3
#define BLADE_D3_AP(A1,A2) (A1 A2)
#define BLADE_D3_OP(A) A
#define bfam_domain_pxest_new bfam_domain_pxest_new_3
#define bfam_domain_pxest_quad_to_glueid \
  bfam_domain_pxest_quad_to_glueid_3
#define bfam_domain_pxest_split_dgx_subdomains \
  bfam_domain_pxest_split_dgx_subdomains_3
#define bfam_domain_pxest_free \
  bfam_domain_pxest_free_3
#define bfam_domain_pxest_t  bfam_domain_pxest_t_3

#else
#error "bad dimension"
#endif


#define BFAM_LOAD_FIELD_RESTRICT_ALIGNED(field,prefix,base,dictionary)         \
bfam_real_t *restrict field;                                                   \
{                                                                              \
  char bfam_load_field_name[BFAM_BUFSIZ];                                      \
  snprintf(bfam_load_field_name,BFAM_BUFSIZ,"%s%s",(prefix),(base));           \
  field = bfam_dictionary_get_value_ptr(dictionary, bfam_load_field_name);     \
  BFAM_ASSERT(field != NULL);                                                  \
}                                                                              \
BFAM_ASSUME_ALIGNED(field,32);
#define BFAM_LOAD_FIELD_ALIGNED(field,prefix,base,dictionary)                  \
bfam_real_t *field;                                                            \
{                                                                              \
  char bfam_load_field_name[BFAM_BUFSIZ];                                      \
  snprintf(bfam_load_field_name,BFAM_BUFSIZ,"%s%s",(prefix),(base));           \
  field = bfam_dictionary_get_value_ptr(dictionary, bfam_load_field_name);     \
  BFAM_ASSERT(field != NULL);                                                  \
}                                                                              \
BFAM_ASSUME_ALIGNED(field,32);

/* The communicator arguments */
const char *comm_args_face_scalars[]      = {NULL};
const char *comm_args_scalars[]           = {"q",NULL};
const char *comm_args_vectors[]           = {NULL};
const char *comm_args_vector_components[] = {NULL};
const char *comm_args_tensors[]           = {NULL};
const char *comm_args_tensor_components[] = {NULL};


/* The structure that we pass around to everything */
typedef struct blade
{
  MPI_Comm mpicomm;
  int      mpirank;
  int      mpisize;

  p4est_connectivity_t *conn;
  bfam_domain_pxest_t  *domain;

  bfam_ts_t                  *blade_ts;
  bfam_subdomain_comm_args_t *comm_args;
} blade_t;

/* Time stepping information */
struct lsrk_table {
  const char *name;
  bfam_ts_lsrk_method_t lsrk_method;
} lsrk_table[] = {
  {"KC54", BFAM_TS_LSRK_KC54},
  {"FE",   BFAM_TS_LSRK_FE},
  {"HEUN", BFAM_TS_LSRK_HEUN},
  {"W33",  BFAM_TS_LSRK_W33},
  {NULL,   BFAM_TS_LSRK_NOOP},
};

struct adams_table {
  const char *name;
  bfam_ts_adams_method_t adams_method;
} adams_table[] = {
  {"Adams 1", BFAM_TS_ADAMS_1},
  {"Adams 2", BFAM_TS_ADAMS_2},
  {"Adams 3", BFAM_TS_ADAMS_3},
  {"Adams 4", BFAM_TS_ADAMS_4},
  {NULL,      BFAM_TS_ADAMS_NOOP},
};

struct local_adams_table {
  const char *name;
  bfam_ts_local_adams_method_t local_adams_method;
} local_adams_table[] = {
  {"Adams 1", BFAM_TS_LOCAL_ADAMS_1},
  {"Adams 2", BFAM_TS_LOCAL_ADAMS_2},
  {"Adams 3", BFAM_TS_LOCAL_ADAMS_3},
  {"Adams 4", BFAM_TS_LOCAL_ADAMS_4},
  {NULL,      BFAM_TS_LOCAL_ADAMS_NOOP},
};

/* the preference structure. The brick_arges is data needed for the p4est brick
 * mesh type */
typedef struct brick_args
{
  int nx;
  int ny;
  int nz;
  int periodic_x;
  int periodic_y;
  int periodic_z;
} brick_args_t;

typedef struct prefs
{
  lua_State *L;

  int dimension;

  char output_prefix[BFAM_BUFSIZ];
  char data_directory[BFAM_BUFSIZ];
  int  vtk_binary;
  int  vtk_compress;

  brick_args_t* brick_args;

  bfam_ts_lsrk_method_t lsrk_method;
  char lsrk_name[BFAM_BUFSIZ];

  bfam_ts_adams_method_t adams_method;
  char adams_name[BFAM_BUFSIZ];

  bfam_ts_local_adams_method_t local_adams_method;
  char local_adams_name[BFAM_BUFSIZ];
} prefs_t;


/*
 * Lua helper functions
 */
/*
 * Helper function for calling a lua function. Based on generic call function of
 * Listing 25.4-25.6 of
 * @book{Ierusalimschy2006Lua,
 *  author = {Ierusalimschy, Roberto},
 *  title = {Programming in Lua, Second Edition},
 *  year = {2006},
 *  isbn = {8590379825},
 *  publisher = {Lua.Org},
 * }
 */
static int
lua_global_function_call(lua_State *L, const char *name, const char *sig, ...)
{
  va_list vl;
  int num_arg = 0;
  int num_res = 0;

  va_start(vl,sig);

  lua_getglobal(L,name);

  if(!lua_isfunction(L,-1))
  {
    BFAM_ROOT_WARNING("function `%s' not found in lua file", name);
    lua_pop(L,1);
    return 1;
  }

  for(num_arg = 0; sig[num_arg] && sig[num_arg] != '>'; num_arg++)
  {
    luaL_checkstack(L,1,"too many arguments");

    switch(sig[num_arg])
    {
      case 'd':
        lua_pushnumber(L,va_arg(vl,double));
        break;
      case 'r':
        lua_pushnumber(L,(double)va_arg(vl,bfam_real_t));
        break;
      case 'i':
        lua_pushinteger(L,va_arg(vl,int));
        break;
      case 's':
        lua_pushstring(L,va_arg(vl,char *));
        break;
      case '>':
        break;
      default:
        BFAM_ABORT("function '%s' invalid input argument (%c)",name,
            sig[num_arg]);
    }
  }

  BFAM_ABORT_IF_NOT(sig[num_arg] == '>',"arguments for '%s' does not contain "
      " a '>' character",name);

  num_res = strlen(sig) - num_arg - 1;

  BFAM_ABORT_IF_NOT(lua_pcall(L,num_arg,num_res,0) == 0,
      "error running function %s: %s", name,lua_tostring(L,-1));

  for(int n = 0; n < num_res;n++)
  {
    switch(sig[num_arg+1+n])
    {
      case 'r':
        BFAM_ABORT_IF_NOT(lua_isnumber(L,n-num_res),
            "for '%s' return %d expected number got '%s'",
            name, n, lua_tostring(L,n-num_res));
        *va_arg(vl, bfam_real_t*) = (bfam_real_t)lua_tonumber(L,n-num_res);
        break;
      case 'i':
        BFAM_ABORT_IF_NOT(lua_isnumber(L,n-num_res),
            "for '%s' return %d expected number got '%s'",
            name, n, lua_tostring(L,n-num_res));
        *va_arg(vl, int*) = lua_tointeger(L,n-num_res);
        break;
      case 's':
        BFAM_ABORT_IF_NOT(lua_isstring(L,n-num_res),
            "for '%s' return %d expected string got '%s'",
            name, n, lua_tostring(L,n-num_res));
        *va_arg(vl, const char **) = lua_tostring(L,n-num_res);
        break;
      default:
        BFAM_ABORT("function '%s' invalid output argument (%c)",name,
            sig[num_arg]);
    }
  }

  lua_pop(L,num_res);

  va_end(vl);
  return 0;
}

static int
lua_get_global_int(lua_State *L, const char *name, int def)
{
  lua_getglobal(L, name);
  int result = def;
  if(!lua_isnumber(L, -1))
    BFAM_ROOT_WARNING("`%s' not found, using default %d", name, def);
  else
    result = (int)lua_tonumber(L, -1);
  lua_pop(L, 1);
  return result;
}

static int
lua_get_table_int(lua_State *L, const char *table, const char *name, int def)
{
  int result = def;
  lua_getglobal(L, table);
  if(!lua_istable(L, -1))
    BFAM_ROOT_WARNING("table `%s' not found, using default %d", table, def);
  else
  {
    lua_pushstring(L,name);
    lua_gettable(L,-2);
    if(!lua_isnumber(L,-1))
      BFAM_ROOT_WARNING("table `%s' does not contain `%s', using default %d",
          table, name, def);
    else
      result = (int)lua_tonumber(L, -1);
    lua_pop(L, 1);
  }
  lua_pop(L, 1);
  return result;
}


static void
init_mpi(blade_t *blade, MPI_Comm mpicomm)
{
  blade->mpicomm = mpicomm;
  BFAM_MPI_CHECK(MPI_Comm_rank(mpicomm, &blade->mpirank));
  BFAM_MPI_CHECK(MPI_Comm_size(mpicomm, &blade->mpisize));
}

/*
 * Set up the blade preference file
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

  prefs->dimension = DIM;

  BFAM_ABORT_IF_NOT(prefs->dimension == 2 || prefs->dimension == 3,
      "blade not set up to handle dimension %d", prefs->dimension);

  /* get output prefix */
  lua_getglobal(L,"output_prefix");
  strncpy(prefs->output_prefix,"solution",BFAM_BUFSIZ);
  if(lua_isstring(L, -1))
  {
    const char *output_prefix = lua_tostring(L, -1);
    strncpy(prefs->output_prefix,output_prefix,BFAM_BUFSIZ);
  }
  else
    BFAM_ROOT_WARNING("using default 'output_prefix': %s",prefs->output_prefix);
  lua_pop(L, 1);

  /* get output prefix */
  lua_getglobal(L,"data_directory");
  strncpy(prefs->data_directory,".",BFAM_BUFSIZ);
  if(lua_isstring(L, -1))
  {
    const char *data_directory = lua_tostring(L, -1);
    strncpy(prefs->data_directory,data_directory,BFAM_BUFSIZ);
  }
  else
    BFAM_ROOT_WARNING("using default 'data_directory': %s",
        prefs->data_directory);
  lua_pop(L, 1);

  /* get vtk_binary */
  lua_getglobal(L,"vtk_binary");
  if(lua_isboolean(L, -1))
  {
    prefs->vtk_binary = lua_toboolean(L, -1);
  }
  else
  {
    prefs->vtk_binary = 1;
    BFAM_ROOT_WARNING("using default 'vtk_binary': %d",
        prefs->vtk_binary);
  }
  lua_pop(L, 1);

  /* get vtk_compress */
  lua_getglobal(L,"vtk_compress");
  if(lua_isboolean(L, -1))
  {
    prefs->vtk_compress = lua_toboolean(L, -1);
  }
  else
  {
    prefs->vtk_compress = 1;
    BFAM_ROOT_WARNING("using default 'vtk_compress': %d",
        prefs->vtk_compress);
  }
  lua_pop(L, 1);


  /* get the time stepper type: we look for both lsrk_method and adams_method */
  prefs->lsrk_method = BFAM_TS_LSRK_NOOP;
  lua_getglobal(L, "lsrk_method");
  if(lua_isstring(L, -1))
  {
    prefs->lsrk_method = lsrk_table[0].lsrk_method;
    strncpy(prefs->lsrk_name,lsrk_table[0].name,BFAM_BUFSIZ);
    if(lua_isstring(L, -1))
    {
      int i;
      const char *lsrk_name = lua_tostring(L, -1);
      for(i = 0; lsrk_table[i].name != NULL; ++i)
      {
        if(strcmp(lsrk_name, lsrk_table[i].name) == 0)
          break;
      }

      if(lsrk_table[i].name == NULL)
        BFAM_ROOT_WARNING("invalid lsrk method name: `%s'; using default %s",
            lsrk_name, prefs->lsrk_name);
      else
      {
        prefs->lsrk_method = lsrk_table[i].lsrk_method;
        strncpy(prefs->lsrk_name,lsrk_table[i].name,BFAM_BUFSIZ);
      }
    }
    else
    {
      BFAM_ROOT_WARNING("`lsrk method' not found, using default: %s",
          prefs->lsrk_name);
    }
    BFAM_ASSERT(prefs->lsrk_method != BFAM_TS_LSRK_NOOP);
  }
  lua_pop(L, 1);

  prefs->adams_method = BFAM_TS_ADAMS_NOOP;
  lua_getglobal(L, "adams_method");
  if(lua_isstring(L, -1))
  {
    prefs->adams_method = adams_table[0].adams_method;
    strncpy(prefs->adams_name,adams_table[0].name,BFAM_BUFSIZ);
    if(lua_isstring(L, -1))
    {
      int i;
      const char *adams_name = lua_tostring(L, -1);
      for(i = 0; adams_table[i].name != NULL; ++i)
      {
        if(strcmp(adams_name, adams_table[i].name) == 0)
          break;
      }

      if(adams_table[i].name == NULL)
        BFAM_ROOT_WARNING("invalid adams method name: `%s'; using default %s",
            adams_name, prefs->adams_name);
      else
      {
        prefs->adams_method = adams_table[i].adams_method;
        strncpy(prefs->adams_name,adams_table[i].name,BFAM_BUFSIZ);
      }
    }
    else
    {
      BFAM_ROOT_WARNING("`adams method' not found, using default: %s",
          prefs->adams_name);
    }
    BFAM_ASSERT(prefs->adams_method != BFAM_TS_ADAMS_NOOP);
  }
  lua_pop(L, 1);

  prefs->local_adams_method = BFAM_TS_LOCAL_ADAMS_NOOP;
  lua_getglobal(L, "local_adams_method");
  if(lua_isstring(L, -1))
  {
    prefs->local_adams_method = local_adams_table[0].local_adams_method;
    strncpy(prefs->local_adams_name,local_adams_table[0].name,BFAM_BUFSIZ);
    if(lua_isstring(L, -1))
    {
      int i;
      const char *local_adams_name = lua_tostring(L, -1);
      for(i = 0; local_adams_table[i].name != NULL; ++i)
      {
        if(strcmp(local_adams_name, local_adams_table[i].name) == 0)
          break;
      }

      if(local_adams_table[i].name == NULL)
        BFAM_ROOT_WARNING("invalid local adams method name: `%s';"
            " using default %s",
            local_adams_name, prefs->local_adams_name);
      else
      {
        prefs->local_adams_method =
          local_adams_table[i].local_adams_method;
        strncpy(prefs->local_adams_name,local_adams_table[i].name,BFAM_BUFSIZ);
      }
    }
    else
    {
      BFAM_ROOT_WARNING("`local_adams method' not found, using default: %s",
          prefs->local_adams_name);
    }
    BFAM_ASSERT(prefs->local_adams_method != BFAM_TS_LOCAL_ADAMS_NOOP);
  }
  lua_pop(L, 1);


  BFAM_ABORT_IF_NOT((prefs->lsrk_method        != BFAM_TS_LSRK_NOOP  &&
                     prefs->adams_method       == BFAM_TS_ADAMS_NOOP &&
                     prefs->local_adams_method == BFAM_TS_LOCAL_ADAMS_NOOP
                    )
                    ||
                    (
                     prefs->lsrk_method        == BFAM_TS_LSRK_NOOP  &&
                     prefs->adams_method       != BFAM_TS_ADAMS_NOOP &&
                     prefs->local_adams_method == BFAM_TS_LOCAL_ADAMS_NOOP
                    )
                    ||
                    (
                     prefs->lsrk_method        == BFAM_TS_LSRK_NOOP  &&
                     prefs->adams_method       == BFAM_TS_ADAMS_NOOP &&
                     prefs->local_adams_method != BFAM_TS_LOCAL_ADAMS_NOOP
                    )
                    ,"must have either LSRK, ADAMS, or LOCAL time stepper");

  /* get the connectivity type */
  prefs->brick_args     = bfam_malloc(sizeof(brick_args_t));

  prefs->brick_args->nx = lua_get_table_int(prefs->L, "brick", "nx",1);
  prefs->brick_args->ny = lua_get_table_int(prefs->L, "brick", "ny",1);

  prefs->brick_args->periodic_x = lua_get_table_int(prefs->L, "brick",
      "periodic_x",1);
  prefs->brick_args->periodic_y = lua_get_table_int(prefs->L, "brick",
      "periodic_y",1);

#if   DIM==2
#elif DIM==3
  prefs->brick_args->nz = lua_get_table_int(prefs->L, "brick", "nz",1);
  prefs->brick_args->periodic_z =
    lua_get_table_int(prefs->L, "brick", "periodic_z",1);
#else
#error "Bad dimension"
#endif

  return prefs;
}

static void
print_prefs(prefs_t *prefs)
{
  BFAM_ROOT_INFO("----------Preferences----------");
  BFAM_ROOT_INFO(" blade Dimension          = %d",prefs->dimension);
  BFAM_ROOT_INFO(" Output prefix            = %s",prefs->output_prefix);
  if(prefs->lsrk_method != BFAM_TS_LSRK_NOOP)
    BFAM_ROOT_INFO(" Low Storage Time Stepper = %s",prefs->lsrk_name);
  if(prefs->adams_method != BFAM_TS_ADAMS_NOOP)
    BFAM_ROOT_INFO(" Adams-Bashforth Scheme   = %s",prefs->adams_name);
  if(prefs->local_adams_method != BFAM_TS_LOCAL_ADAMS_NOOP)
    BFAM_ROOT_INFO(" Local Adams-Bashforth Scheme   = %s",
        prefs->local_adams_name);
  if(prefs->brick_args != NULL)
  {
    BFAM_ROOT_INFO(" brick arguments");
    BFAM_ROOT_INFO("  nx         = %d", prefs->brick_args->nx);
    BFAM_ROOT_INFO("  ny         = %d", prefs->brick_args->ny);
#if   DIM==2
#elif DIM==3
    BFAM_ROOT_INFO("  nz         = %d", prefs->brick_args->nz);
#else
#error "Bad dimension"
#endif
    BFAM_ROOT_INFO("  periodic_x = %d", prefs->brick_args->periodic_x);
    BFAM_ROOT_INFO("  periodic_y = %d", prefs->brick_args->periodic_y);
#if   DIM==2
#elif DIM==3
    BFAM_ROOT_INFO("  periodic_z = %d", prefs->brick_args->periodic_z);
#else
#error "Bad dimension"
#endif
  }
  BFAM_ROOT_INFO("-------------------------------");
}


static void
free_prefs(prefs_t *prefs)
{
  if(prefs->brick_args != NULL) bfam_free(prefs->brick_args);
  lua_close(prefs->L);
}


static int
static_refine_fn(p4est_t * p4est, p4est_topidx_t which_tree,
    p4est_quadrant_t * quadrant)
{
  return 1;
}

static int
refine_fn(p4est_t * p4est, p4est_topidx_t which_tree,
    p4est_quadrant_t * quadrant)
{
  lua_State *L = (lua_State*)p4est->user_pointer;

  double vxyz[P4EST_CHILDREN*3];

#if DIM==2
  for(int iy = 0; iy < 2; iy++)
    for(int ix = 0; ix < 2; ix++)
    {
      int ox = ix*(1 << (P4EST_MAXLEVEL-quadrant->level));
      int oy = iy*(1 << (P4EST_MAXLEVEL-quadrant->level));
      p4est_qcoord_to_vertex (p4est->connectivity, which_tree,
          quadrant->x+ox, quadrant->y+oy,&vxyz[3*(ix + iy*2)]);
    }

  int val = 0;
  int result = lua_global_function_call(L,"refinement_function",
      "ddd" "ddd" "ddd" "ddd" "ii" ">" "i",
      vxyz[0], vxyz[1], vxyz[2], vxyz[3], vxyz[ 4], vxyz[ 5],
      vxyz[6], vxyz[7], vxyz[8], vxyz[9], vxyz[10], vxyz[11],
      quadrant->level, which_tree, &val);
  BFAM_ABORT_IF_NOT(result == 0, "Expected 0 got %d", result);
#elif DIM==3
  for(int iz = 0; iz < 2; iz++)
    for(int iy = 0; iy < 2; iy++)
      for(int ix = 0; ix < 2; ix++)
    {
      int ox = ix*(1 << (P4EST_MAXLEVEL-quadrant->level));
      int oy = iy*(1 << (P4EST_MAXLEVEL-quadrant->level));
      int oz = iz*(1 << (P4EST_MAXLEVEL-quadrant->level));
      p4est_qcoord_to_vertex (p4est->connectivity, which_tree,
          quadrant->x+ox, quadrant->y+oy, quadrant->z+oz,
          &vxyz[3*(ix + 2*(iy + 2*iz))]);
    }

  int val = 0;
  int result = lua_global_function_call(L,"refinement_function",
      "ddd" "ddd" "ddd" "ddd" "ddd" "ddd" "ddd" "ddd" "ii" ">" "i",
      vxyz[ 0], vxyz[ 1], vxyz[ 2], vxyz[ 3], vxyz[ 4], vxyz[ 5],
      vxyz[ 6], vxyz[ 7], vxyz[ 8], vxyz[ 9], vxyz[10], vxyz[11],
      vxyz[12], vxyz[13], vxyz[14], vxyz[15], vxyz[16], vxyz[17],
      vxyz[18], vxyz[19], vxyz[20], vxyz[21], vxyz[22], vxyz[23],
      quadrant->level, which_tree, &val);
  BFAM_ABORT_IF_NOT(result == 0, "Expected 0 got %d", result);
#else
#error "Bad Dimension"
#endif

  return val;
}

/*
 * data structure used to split the domain into subdomains. Namely we loop
 * through all the local elements and assign a subdomain ID. Right now the
 * subdomain IDs correspond to the element orders. This could be further split
 * in the future to allow for further splitting.
 */
typedef struct split_iter_data
{
  lua_State *L; /* Lua file */
  int loc;      /* counter for which element I am processing */
  bfam_locidx_t *subdomain_id; /*array to state the subdomain ID of each cell */
  bfam_locidx_t  max_N;        /* largest polynomial order found */
} split_iter_data_t;

static void
get_element_order(p4est_iter_volume_info_t *info, void *arg)
{
  split_iter_data_t *data = (split_iter_data_t*) arg;

  double vxyz[P4EST_CHILDREN*3];

#if   DIM==2
  for(int ix = 0; ix < 2; ix++)
    for(int iy = 0; iy < 2; iy++)
    {
      int ox = ix*(1 << (P4EST_MAXLEVEL-info->quad->level));
      int oy = iy*(1 << (P4EST_MAXLEVEL-info->quad->level));
      p4est_qcoord_to_vertex (info->p4est->connectivity, info->treeid,
          info->quad->x+ox, info->quad->y+oy,&vxyz[3*(ix + 2*iy)]);
    }

  int N = 0;
  int result = lua_global_function_call(data->L,"element_order",
      "ddd" "ddd" "ddd" "ddd" "ii" ">" "i",
      vxyz[0], vxyz[1], vxyz[2], vxyz[3], vxyz[ 4], vxyz[ 5],
      vxyz[6], vxyz[7], vxyz[8], vxyz[9], vxyz[10], vxyz[11],
      info->quad->level, info->treeid, &N);
#elif DIM==3
  for(int iz = 0; iz < 2; iz++)
    for(int iy = 0; iy < 2; iy++)
      for(int ix = 0; ix < 2; ix++)
    {
      int ox = ix*(1 << (P4EST_MAXLEVEL-info->quad->level));
      int oy = iy*(1 << (P4EST_MAXLEVEL-info->quad->level));
      int oz = iz*(1 << (P4EST_MAXLEVEL-info->quad->level));
      p4est_qcoord_to_vertex (info->p4est->connectivity, info->treeid,
          info->quad->x+ox, info->quad->y+oy, info->quad->z+oz,
          &vxyz[3*(ix + 2*(iy + 2*iz))]);
    }

  int N = 0;
  int result = lua_global_function_call(data->L,"element_order",
      "ddd" "ddd" "ddd" "ddd" "ddd" "ddd" "ddd" "ddd" "ii" ">" "i",
      vxyz[ 0], vxyz[ 1], vxyz[ 2], vxyz[ 3], vxyz[ 4], vxyz[ 5],
      vxyz[ 6], vxyz[ 7], vxyz[ 8], vxyz[ 9], vxyz[10], vxyz[11],
      vxyz[12], vxyz[13], vxyz[14], vxyz[15], vxyz[16], vxyz[17],
      vxyz[18], vxyz[19], vxyz[20], vxyz[21], vxyz[22], vxyz[23],
      info->quad->level, info->treeid, &N);
#else
#error "Bad Dimension"
#endif

  BFAM_ABORT_IF_NOT(result == 0, "no 'element_order' function");
  BFAM_ABORT_IF_NOT(N > 0, "N must be greater than 0");
  data->subdomain_id[data->loc] = N-1;
  data->max_N = BFAM_MAX(N,data->max_N);
  data->loc++;
}

/*
 * function to figure out what the order of that is on this processor is.
 */
static void
split_domain(blade_t *blade, prefs_t *prefs)
{
  bfam_domain_pxest_t *domain = blade->domain;
  if(domain->pxest->local_num_quadrants == 0) BFAM_ABORT("No quadrants");
  bfam_locidx_t *sub_ids =
    bfam_malloc(domain->pxest->local_num_quadrants*sizeof(bfam_locidx_t));

  split_iter_data_t data = {prefs->L,0,sub_ids,0};
  p4est_iterate (domain->pxest, NULL, &data, get_element_order, NULL, NULL
#if DIM==3
      ,NULL
#endif
      );

  bfam_locidx_t *N = bfam_malloc(data.max_N*sizeof(bfam_locidx_t));
  for(int n = 0; n < data.max_N;n++) N[n] = n+1;

  /* Last argument is NULL since we don't need user specified glues for blade */
  bfam_domain_pxest_split_dgx_subdomains(domain, data.max_N, sub_ids, N, NULL);

  bfam_free(sub_ids);
  bfam_free(N);
}

static void
field_set_val(bfam_locidx_t npoints, const char *name, bfam_real_t t,
    bfam_real_t *restrict x, bfam_real_t *restrict y, bfam_real_t *restrict z,
    struct bfam_subdomain *s, void *arg, bfam_real_t *restrict field)
{
  BFAM_ASSUME_ALIGNED(x, 32);
  BFAM_ASSUME_ALIGNED(y, 32);
  BFAM_ASSUME_ALIGNED(z, 32);
  BFAM_ASSUME_ALIGNED(field, 32);
  bfam_real_t val = 0;

  lua_State *L = (lua_State*) arg;
  lua_getglobal(L,name);

  if(lua_isfunction(L,-1))
  {
    BFAM_ROOT_INFO("field_set_val: using '%s' as lua callback function",name);
    lua_pop(L,1);
#if   DIM==2
    bfam_real_t tmpz = 0;
#endif
    for(bfam_locidx_t n=0; n < npoints; ++n)
#if   DIM==2
      lua_global_function_call(L, name, "rrrr>r", x[n],y[n],tmpz,t,&field[n]);
#elif DIM==3
      lua_global_function_call(L, name, "rrrr>r", x[n],y[n],z[n],t,&field[n]);
#else
#error "Bad Dimension"
#endif
    return;
  }
  else if(lua_isnumber(L,-1)) val = (bfam_real_t)lua_tonumber(L,-1);
  else BFAM_WARNING("Did not find '%s' in lua as a function or number using 0: "
                    "lua message: '%s'", name,lua_tostring(L,-1));
  lua_pop(L,1);

  BFAM_ROOT_INFO("field_set_val: using '%s' with value %"BFAM_REAL_FMTe,
      name,val);

  for(bfam_locidx_t n=0; n < npoints; ++n)
    field[n] = val;
}

static void
domain_add_fields(blade_t *blade, prefs_t *prefs)
{
  const char *volume[] = {"_volume",NULL};
  const char *fields[] = {"u", "q", NULL};
  bfam_domain_t *domain = (bfam_domain_t*)blade->domain;
  for(int f = 0; fields[f] != NULL; f++)
  {
    bfam_domain_add_field (domain, BFAM_DOMAIN_OR, volume, fields[f]);
    bfam_domain_init_field(domain, BFAM_DOMAIN_OR, volume, fields[f], 0,
        field_set_val, prefs->L);
  }

  /* add glue fields */
  const char *glue[]   = {"_glue_parallel", "_glue_local", NULL};
  for(int f = 0 ; comm_args_scalars[f] != NULL; f++)
  {
    bfam_domain_add_minus_field(domain, BFAM_DOMAIN_OR, glue,
        comm_args_scalars[f]);
    bfam_domain_add_plus_field( domain, BFAM_DOMAIN_OR, glue,
        comm_args_scalars[f]);
  }

//JK   /* exchange material properties to glue */
//JK   const char *glue_mat[] = {"u","_grid_x0","_grid_x1",
//JK #if DIM==3
//JK     "_grid_x2",
//JK #endif
//JK     NULL};
//JK   for(bfam_locidx_t g = 0; glue_mat[g] != NULL; g++)
//JK   {
//JK     bfam_domain_add_minus_field(domain, BFAM_DOMAIN_OR, glue, glue_mat[g]);
//JK     bfam_domain_add_plus_field( domain, BFAM_DOMAIN_OR, glue, glue_mat[g]);
//JK   }
//JK   bfam_domain_add_field(domain, BFAM_DOMAIN_OR, glue, "_grid_x0");
//JK   bfam_domain_add_field(domain, BFAM_DOMAIN_OR, glue, "_grid_x1");
//JK #if DIM==3
//JK   bfam_domain_add_field(domain, BFAM_DOMAIN_OR, glue, "_grid_x2");
//JK #endif
//JK 
//JK   const char *glue_face_scalar[] = {"_grid_nx0","_grid_nx1",
//JK #if DIM==3
//JK     "_grid_nx2",
//JK #endif
//JK     NULL};
//JK   for(bfam_locidx_t g = 0; glue_face_scalar[g] != NULL; g++)
//JK   {
//JK     bfam_domain_add_minus_field(domain, BFAM_DOMAIN_OR, glue,
//JK         glue_face_scalar[g]);
//JK     bfam_domain_add_plus_field( domain, BFAM_DOMAIN_OR, glue,
//JK         glue_face_scalar[g]);
//JK   }
//JK 
//JK   bfam_communicator_t material_comm;
//JK 
//JK   bfam_subdomain_comm_args_t mat_args;
//JK   const char * mat_NULL[]      = {NULL};
//JK   mat_args.scalars_m           = glue_mat;
//JK   mat_args.vectors_m           = mat_NULL;
//JK   mat_args.vector_components_m = mat_NULL;
//JK   mat_args.tensors_m           = mat_NULL;
//JK   mat_args.tensor_components_m = mat_NULL;
//JK   mat_args.face_scalars_m      = glue_face_scalar;
//JK 
//JK   mat_args.scalars_p           = glue_mat;
//JK   mat_args.vectors_p           = mat_NULL;
//JK   mat_args.vector_components_p = mat_NULL;
//JK   mat_args.tensors_p           = mat_NULL;
//JK   mat_args.tensor_components_p = mat_NULL;
//JK   mat_args.face_scalars_p      = glue_face_scalar;
//JK 
//JK   mat_args.user_comm_info       = NULL;
//JK   mat_args.user_put_send_buffer = NULL;
//JK   mat_args.user_get_recv_buffer = NULL;
//JK   mat_args.user_data            = NULL;
//JK 
//JK   mat_args.user_prefix_function = NULL;
//JK 
//JK 
//JK   bfam_communicator_init(&material_comm,domain,BFAM_DOMAIN_OR,glue,
//JK       blade->mpicomm,10,&mat_args);
//JK   bfam_communicator_start( &material_comm);
//JK   bfam_communicator_finish(&material_comm);
//JK   bfam_communicator_free(  &material_comm);
//JK 
//JK 
//JK   /* we can trick init fields into handling locations to glue */
//JK   bfam_domain_init_field(domain, BFAM_DOMAIN_OR, glue, "_grid_x0", 0,
//JK       blade_grid_glue, NULL);
//JK 
//JK   const char *boundary[] = {"_glue_boundary",NULL};
//JK   const char *boundary_fields[] = {"_grid_x0","_grid_x1",
//JK #if DIM==3
//JK     "_grid_x2",
//JK #endif
//JK     NULL};
//JK   for(bfam_locidx_t g = 0; boundary_fields[g] != NULL; g++)
//JK     bfam_domain_add_field(domain, BFAM_DOMAIN_OR, boundary, boundary_fields[g]);
//JK 
//JK   bfam_domain_init_field(domain, BFAM_DOMAIN_OR, boundary, "_grid_x1", 0,
//JK       blade_grid_boundary, NULL);
//JK 
//JK   /*
//JK    * Loop through glue info
//JK    */
//JK 
//JK   lua_State *L = prefs->L;
//JK #ifdef BFAM_DEBUG
//JK   int top = lua_gettop(L);
//JK #endif
//JK 
//JK   lua_getglobal(L,"glue_info");
//JK   if(!lua_istable(L, -1))
//JK     BFAM_ROOT_WARNING("table `%s' not found", "glue_info");
//JK   else
//JK   {
//JK     luaL_checktype(L, -1, LUA_TTABLE);
//JK 
//JK     int N_glueids = luaL_getn(L, -1);
//JK     BFAM_LDEBUG("N_glueids: %3d", N_glueids);
//JK 
//JK     BFAM_ROOT_VERBOSE("Reading glue info from lua");
//JK     for(int i = 1; i <= N_glueids; ++i)
//JK     {
//JK       BFAM_ROOT_VERBOSE("glue ID: %d", i);
//JK       /*
//JK        * get glue_id tag
//JK        */
//JK       const char *this_glue[2];
//JK       char this_glue_tag[BFAM_BUFSIZ];
//JK       this_glue[0] = this_glue_tag;
//JK       this_glue[1] = NULL;
//JK       snprintf(this_glue_tag,BFAM_BUFSIZ,"_glue_id_%jd",(intmax_t)i);
//JK 
//JK       lua_rawgeti(L, -1, i);
//JK       luaL_checktype(L, -1, LUA_TTABLE);
//JK 
//JK       lua_pushstring(L, "type");
//JK       lua_gettable(L, -2);
//JK       luaL_checktype(L, -1, LUA_TSTRING);
//JK 
//JK       size_t match_len = BFAM_MIN(lua_strlen(L, -1),8);
//JK 
//JK       glue_info_type_t type = UNKNOWN;
//JK       if(0==strncmp(lua_tostring(L, -1), "friction", match_len))
//JK         type = FRICTION;
//JK       else if(0==strncmp(lua_tostring(L, -1), "boundary", match_len))
//JK         type = BOUNDARY;
//JK 
//JK       BFAM_ROOT_VERBOSE("  type: %d", type);
//JK       lua_pop(L, 1);
//JK 
//JK       lua_pushstring(L, "tag");
//JK       lua_gettable(L, -2);
//JK       if(lua_isstring(L, -1))
//JK       {
//JK         const char * tag = lua_tostring(L, -1);
//JK 
//JK         bfam_domain_add_tag((bfam_domain_t*)blade->domain, BFAM_DOMAIN_OR,
//JK             this_glue, tag);
//JK       }
//JK       else
//JK         BFAM_WARNING("No tag for glue_id %d", i);
//JK       lua_pop(L, 1);
//JK 
//JK       if(type==FRICTION)
//JK       {
//JK         for(int f = 0; sw_fields[f] != NULL; ++f)
//JK         {
//JK           bfam_real_t value = 0;
//JK 
//JK           lua_pushstring(L,sw_fields[f]);
//JK           lua_gettable(L,-2);
//JK           if(!lua_isnumber(L,-1))
//JK             BFAM_ROOT_WARNING(
//JK                 " glue %d does not contain `%s', using default %"BFAM_REAL_PRIe,
//JK                 i, sw_fields[f], value);
//JK           else
//JK             value = (bfam_real_t)lua_tonumber(L, -1);
//JK           lua_pop(L, 1);
//JK 
//JK           bfam_domain_add_field(domain, BFAM_DOMAIN_OR, this_glue, sw_fields[f]);
//JK           bfam_domain_init_field(domain, BFAM_DOMAIN_OR, this_glue, sw_fields[f],
//JK               0, field_set_const, &value);
//JK         }
//JK         bfam_domain_init_field(domain, BFAM_DOMAIN_OR, this_glue, "Tp1_0",
//JK             0, field_set_friction_init_stress, NULL);
//JK       }
//JK 
//JK 
//JK       lua_pop(L, 1);
//JK     }
//JK   }
//JK 
//JK   /*
//JK    * Set default boundary condition
//JK    */
//JK   const char *default_boundary_glue[] = {"_glue_boundary","_glue_id_-1", NULL};
//JK   bfam_domain_add_tag((bfam_domain_t*)blade->domain, BFAM_DOMAIN_AND,
//JK       default_boundary_glue, prefs->default_boundary_tag);
//JK 
//JK 
//JK   lua_pop(L,-1);
//JK   BFAM_ASSERT(top == lua_gettop(L));
}

/*
 * This function initializes the domain in bfam
 */
static void
init_domain(blade_t *blade, prefs_t *prefs)
{
  /* Set up the connectivity */
  if(prefs->brick_args != NULL)
    blade->conn = p4est_connectivity_new_brick(
        prefs->brick_args->nx,
        prefs->brick_args->ny,
#if DIM==3
        prefs->brick_args->nz,
#endif
        prefs->brick_args->periodic_x,
        prefs->brick_args->periodic_y
#if DIM==3
        ,prefs->brick_args->periodic_z
#endif
        );
  else BFAM_ABORT("no connectivity");

  /* let the user modify the connectivity vertices */
  for(int i = 0; i < blade->conn->num_vertices; i++)
  {
    bfam_real_t x = blade->conn->vertices[i*3+0];
    bfam_real_t y = blade->conn->vertices[i*3+1];
    bfam_real_t z = blade->conn->vertices[i*3+2];
    int result = lua_global_function_call(prefs->L,"connectivity_vertices",
        "rrr>rrr",x,y,z,&x,&y,&z);
    if(result != 0) break;
    blade->conn->vertices[i*3+0] = x;
    blade->conn->vertices[i*3+1] = y;
    blade->conn->vertices[i*3+2] = z;
  }
  /* create the domain */
  blade->domain = bfam_domain_pxest_new(blade->mpicomm, blade->conn);

  /* call user refinement function */
  lua_getglobal(prefs->L,"refinement_function");
  if(lua_isfunction(prefs->L,-1))
  {
    lua_pop(prefs->L, 1);
    void *current_user_pointer = blade->domain->pxest->user_pointer;
    blade->domain->pxest->user_pointer = prefs->L;
    p4est_refine(blade->domain->pxest, 1, refine_fn, NULL);

    blade->domain->pxest->user_pointer = current_user_pointer;
  }
  else BFAM_ROOT_WARNING("function `%s' not found in lua file",
      "refinement_function");

  p4est_balance(blade->domain->pxest, P4EST_CONNECT_CORNER, NULL);

  /*
   * This is to statically refine all cells of a balanced mesh, since it's
   * already balanced it will remain balanced
   */
  int stat_ref = lua_get_global_int(prefs->L, "static_refinement", 0);
  for(int i = 0; i < stat_ref; i++)
    p4est_refine(blade->domain->pxest, 0, static_refine_fn, NULL);

  p4est_partition(blade->domain->pxest, NULL);

  /* split the domain */
  split_domain(blade,prefs);

  /* add fields to the domnain */
  domain_add_fields(blade,prefs);
}

static void
stop_blade(blade_t *blade,prefs_t *prefs)
{
  //JK bfam_free(blade->comm_args);
  //JK if(prefs->lsrk_method != BFAM_TS_LSRK_NOOP)
  //JK   bfam_ts_lsrk_free((bfam_ts_lsrk_t*) blade->blade_ts);
  //JK if(prefs->adams_method != BFAM_TS_ADAMS_NOOP)
  //JK   bfam_ts_adams_free((bfam_ts_adams_t*) blade->blade_ts);
  //JK if(prefs->local_adams_method != BFAM_TS_LOCAL_ADAMS_NOOP)
  //JK   bfam_ts_local_adams_free((bfam_ts_local_adams_t*) blade->blade_ts);
  //JK bfam_free(blade->blade_ts);
  bfam_domain_pxest_free(blade->domain);
  bfam_free(blade->domain);
  p4est_connectivity_destroy(blade->conn);
}

/*
 * run the blade
 */
static void
run(MPI_Comm mpicomm, prefs_t *prefs)
{
  blade_t blade;
  blade.conn            = NULL;
  blade.domain          = NULL;
  blade.blade_ts        = NULL;
  blade.comm_args       = NULL;

  init_mpi(&blade, mpicomm);

  init_domain(&blade, prefs);

  //JK run_simulation(&blade, prefs);

  stop_blade(&blade,prefs);
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

#endif
