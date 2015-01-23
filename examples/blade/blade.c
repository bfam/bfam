#include <bfam.h>

#ifdef BLADE_DGX_DIMENSION

/* Handle the header files */
#if   BLADE_DGX_DIMENSION==2
#include <bfam_domain_pxest_2.h>
#include "blade_dgx_rhs_2.h"
#include <p4est_iterate.h>
#elif BLADE_DGX_DIMENSION==3
#include <bfam_domain_pxest_3.h>
#include "blade_dgx_rhs_3.h"
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
  bfam_domain_pxest_split_dgx_subdomains(domain, data.max_N, sub_ids, NULL, N,
      NULL, NULL, NULL);

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
#if DIM==2
  const char *fields[] = {"ux", "uy", "q", NULL};
#elif DIM==3
  const char *fields[] = {"ux", "uy", "uz", "q", NULL};
#else
#error "Bad Dimension"
#endif
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
//JK     int N_glueids = lua_objlen(L, -1);
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

  p4est_partition(blade->domain->pxest, 1, NULL);

  /* split the domain */
  split_domain(blade,prefs);

  /* add fields to the domnain */
  domain_add_fields(blade,prefs);
}

static void
stop_blade(blade_t *blade,prefs_t *prefs)
{
  bfam_free(blade->comm_args);
  if(prefs->lsrk_method != BFAM_TS_LSRK_NOOP)
    bfam_ts_lsrk_free((bfam_ts_lsrk_t*) blade->blade_ts);
  if(prefs->adams_method != BFAM_TS_ADAMS_NOOP)
    bfam_ts_adams_free((bfam_ts_adams_t*) blade->blade_ts);
  if(prefs->local_adams_method != BFAM_TS_LOCAL_ADAMS_NOOP)
    bfam_ts_local_adams_free((bfam_ts_local_adams_t*) blade->blade_ts);
  bfam_free(blade->blade_ts);
  bfam_domain_pxest_free(blade->domain);
  bfam_free(blade->domain);
  p4est_connectivity_destroy(blade->conn);
}

static void
compute_subdomain_dt(bfam_locidx_t npoints, const char *name, bfam_real_t time,
    bfam_real_t *restrict x, bfam_real_t *restrict y, bfam_real_t *restrict z,
    struct bfam_subdomain *s, void *arg, bfam_real_t *restrict JI)
{
  BFAM_ASSUME_ALIGNED(x, 32);
  BFAM_ASSUME_ALIGNED(y, 32);
  BFAM_ASSUME_ALIGNED(z, 32);
  BFAM_ASSUME_ALIGNED(JI, 32);
  bfam_real_t *dt = (bfam_real_t*)arg;

  bfam_real_t *restrict Jr0x0 =
    bfam_dictionary_get_value_ptr(&s->fields, "_grid_Jr0x0");
  BFAM_ASSUME_ALIGNED(Jr0x0, 32);

  bfam_real_t *restrict Jr0x1 =
    bfam_dictionary_get_value_ptr(&s->fields, "_grid_Jr0x1");
  BFAM_ASSUME_ALIGNED(Jr0x1, 32);

#if DIM==3
  bfam_real_t *restrict Jr0x2 =
    bfam_dictionary_get_value_ptr(&s->fields, "_grid_Jr0x2");
  BFAM_ASSUME_ALIGNED(Jr0x2, 32);
#endif

  bfam_real_t *restrict Jr1x0 =
    bfam_dictionary_get_value_ptr(&s->fields, "_grid_Jr1x0");
  BFAM_ASSUME_ALIGNED(Jr1x0, 32);

  bfam_real_t *restrict Jr1x1 =
    bfam_dictionary_get_value_ptr(&s->fields, "_grid_Jr1x1");
  BFAM_ASSUME_ALIGNED(Jr1x1, 32);

#if DIM==3
  bfam_real_t *restrict Jr1x2 =
    bfam_dictionary_get_value_ptr(&s->fields, "_grid_Jr1x2");
  BFAM_ASSUME_ALIGNED(Jr1x2, 32);
#endif

#if DIM==3
  bfam_real_t *restrict Jr2x0 =
    bfam_dictionary_get_value_ptr(&s->fields, "_grid_Jr2x0");
  BFAM_ASSUME_ALIGNED(Jr2x0, 32);

  bfam_real_t *restrict Jr2x1 =
    bfam_dictionary_get_value_ptr(&s->fields, "_grid_Jr2x1");
  BFAM_ASSUME_ALIGNED(Jr2x1, 32);

  bfam_real_t *restrict Jr2x2 =
    bfam_dictionary_get_value_ptr(&s->fields, "_grid_Jr2x2");
  BFAM_ASSUME_ALIGNED(Jr2x2, 32);
#endif

  bfam_real_t *restrict ux =
    bfam_dictionary_get_value_ptr(&s->fields, "ux");
  BFAM_ASSUME_ALIGNED(ux, 32);

  bfam_real_t *restrict uy =
    bfam_dictionary_get_value_ptr(&s->fields, "uy");
  BFAM_ASSUME_ALIGNED(uy, 32);

#if DIM==3
  bfam_real_t *restrict uz =
    bfam_dictionary_get_value_ptr(&s->fields, "uz");
  BFAM_ASSUME_ALIGNED(uz, 32);
#endif

  bfam_subdomain_dgx_t *sub = (bfam_subdomain_dgx_t *) s;
  bfam_real_t p2 = sub->N*sub->N;

  for(int n = 0; n < npoints; n++)
  {
#if   DIM==2
    bfam_real_t u = BFAM_REAL_SQRT(ux[n]*ux[n] + uy[n]*uy[n]);
    bfam_real_t hr1 =
      BFAM_REAL(1.0)/BFAM_REAL_SQRT(JI[n]*JI[n]*(Jr1x1[n]*Jr1x1[n]
                                                +Jr1x0[n]*Jr1x0[n]));
    bfam_real_t hr0 =
      BFAM_REAL(1.0)/BFAM_REAL_SQRT(JI[n]*JI[n]*(Jr0x1[n]*Jr0x1[n]
                                                +Jr0x0[n]*Jr0x0[n]));
    dt[0] = BFAM_MIN(dt[0], BFAM_MIN(hr1,hr0)/u/p2);
#elif DIM==3
    bfam_real_t u = BFAM_REAL_SQRT(ux[n]*ux[n] + uy[n]*uy[n] + uz[n]*uz[n]);
    bfam_real_t hr2 =
      BFAM_REAL(1.0)/BFAM_REAL_SQRT(JI[n]*JI[n]*(Jr2x2[n]*Jr2x2[n]
                                                +Jr2x1[n]*Jr2x1[n]
                                                +Jr2x0[n]*Jr2x0[n]));
    bfam_real_t hr1 =
      BFAM_REAL(1.0)/BFAM_REAL_SQRT(JI[n]*JI[n]*(Jr1x2[n]*Jr1x2[n]
                                                +Jr1x1[n]*Jr1x1[n]
                                                +Jr1x0[n]*Jr1x0[n]));
    bfam_real_t hr0 =
      BFAM_REAL(1.0)/BFAM_REAL_SQRT(JI[n]*JI[n]*(Jr0x2[n]*Jr0x2[n]
                                                +Jr0x1[n]*Jr0x1[n]
                                                +Jr0x0[n]*Jr0x0[n]));
    dt[0] = BFAM_MIN(dt[0], BFAM_MIN(BFAM_MIN(hr2,hr1),hr0)/u/p2);
#else
#error "bad dimension"
#endif
  }
}

static void
time_level_comm_info(bfam_subdomain_t* thisSubdomain, size_t *send_sz,
    size_t *recv_sz, void *comm_args)
{
  BFAM_ASSERT(comm_args);
  *send_sz += sizeof(int);
  *recv_sz += sizeof(int);
}

static void
time_level_put_send_buffer(bfam_subdomain_t* thisSubdomain, void *buffer,
    size_t send_sz, void *comm_args)
{
  /* just some sanity check */
  BFAM_ASSERT(send_sz == sizeof(int));
  BFAM_ASSERT(comm_args);
  BFAM_ASSERT(thisSubdomain->glue_m);
  BFAM_ASSERT(thisSubdomain->glue_m->sub_m);

  /* determine the level of the minus side and add tag to minus side */
  bfam_subdomain_comm_args_t *args =
    (bfam_subdomain_comm_args_t*) comm_args;
  BFAM_ASSERT(args->user_data);
  int max_lvl = *(int*)args->user_data;

  int lvl = 0;
  char tag[BFAM_BUFSIZ];
  for(; lvl <= max_lvl;lvl++)
  {
    bfam_ts_local_adams_fill_level_tag(tag,BFAM_BUFSIZ,lvl);
    if(bfam_subdomain_has_tag(thisSubdomain->glue_m->sub_m,tag)) break;
  }
  BFAM_ABORT_IF(lvl > max_lvl,"glue %s: "
      "max number of levels searched in minus side %s "
      "and no level tag found",
      thisSubdomain->name,thisSubdomain->glue_m->sub_m->name);

  bfam_subdomain_minus_add_tag(thisSubdomain,tag);

  /* store the level */
  *(int*)buffer = lvl;
}

static void
time_level_get_recv_buffer(bfam_subdomain_t *thisSubdomain, void *buffer,
    size_t recv_sz, void *comm_args)
{
  BFAM_ASSERT(comm_args);
  BFAM_ASSERT(recv_sz == sizeof(int));
  int lvl = *(int*)buffer;
  char tag[BFAM_BUFSIZ];
  bfam_ts_local_adams_fill_level_tag(tag,BFAM_BUFSIZ,lvl);
  bfam_subdomain_plus_add_tag(thisSubdomain,tag);
}

/* These are the functions required by the time steppers */
static void
field_zero(bfam_locidx_t npoints, const char *name, bfam_real_t time,
    bfam_real_t *restrict x, bfam_real_t *restrict y, bfam_real_t *restrict z,
    struct bfam_subdomain *s, void *arg, bfam_real_t *restrict field)
{
  BFAM_ASSUME_ALIGNED(x, 32);
  BFAM_ASSUME_ALIGNED(y, 32);
  BFAM_ASSUME_ALIGNED(z, 32);
  BFAM_ASSUME_ALIGNED(field, 32);

  for(bfam_locidx_t n=0; n < npoints; ++n) field[n] = 0;
}

void aux_rates (bfam_subdomain_t *thisSubdomain, const char *prefix)
{
  if(bfam_subdomain_has_tag(thisSubdomain,"_volume"))
  {
    const char *fields[] =
      {"q",NULL};

    char field[BFAM_BUFSIZ];
    for(int f = 0; fields[f]!=NULL; ++f)
    {
      snprintf(field,BFAM_BUFSIZ,"%s%s",prefix,fields[f]);
      thisSubdomain->field_add(thisSubdomain,field);
      bfam_subdomain_field_init(thisSubdomain, field, 0, field_zero, NULL);
    }
  }
}

void scale_rates_advection (bfam_subdomain_dgx_t *sub, const char *rate_prefix,
    const bfam_long_real_t a)
{
#if  DIM==2
#define X(order) \
  case order: blade_dgx_scale_rates_advection_2_##order(sub->N,sub, \
                  rate_prefix,a); break;
#elif  DIM==3
#define X(order) \
  case order: blade_dgx_scale_rates_advection_3_##order(sub->N,sub, \
                  rate_prefix,a); break;
#else
#error "Bad Dimension"
#endif

  switch(sub->N)
  {
    BFAM_LIST_OF_DGX_NORDERS
    default:
#if   DIM==2
      blade_dgx_scale_rates_advection_2_(sub->N,sub,rate_prefix,a);
#elif DIM==3
      blade_dgx_scale_rates_advection_3_(sub->N,sub,rate_prefix,a);
#else
#error "Bad Dimension"
#endif

      break;
  }
#undef X
}

void scale_rates (bfam_subdomain_t *thisSubdomain, const char *rate_prefix,
    const bfam_long_real_t a)
{
  BFAM_ASSERT(rate_prefix);
  BFAM_ASSERT(bfam_subdomain_has_tag(thisSubdomain,"_subdomain_dgx"));
  if(bfam_subdomain_has_tag(thisSubdomain,"_volume"))
  {
    bfam_subdomain_dgx_t *sub = (bfam_subdomain_dgx_t*)thisSubdomain;
    scale_rates_advection(sub,rate_prefix,a);
  }
  else if(bfam_subdomain_has_tag(thisSubdomain,"_glue_boundary")
      ||  bfam_subdomain_has_tag(thisSubdomain,"_glue_parallel")
      ||  bfam_subdomain_has_tag(thisSubdomain,"_glue_local"   ));
  else
    BFAM_ABORT("Unknown subdomain: %s",thisSubdomain->name);
}

static void
intra_rhs_advection(int N, bfam_subdomain_dgx_t *sub,
    const char *rate_prefix, const char *field_prefix, const bfam_long_real_t t)
{
#if  DIM==2
#define X(order) \
  case order: blade_dgx_intra_rhs_advection_2_##order(N,sub, \
                  rate_prefix,field_prefix,t); break;
#elif  DIM==3
#define X(order) \
  case order: blade_dgx_intra_rhs_advection_3_##order(N,sub, \
                  rate_prefix,field_prefix,t); break;
#else
#error "Bad Dimension"
#endif

  switch(N)
  {
    BFAM_LIST_OF_DGX_NORDERS
    default:
#if   DIM==2
      blade_dgx_intra_rhs_advection_2_(N,sub,rate_prefix, field_prefix,t);
#elif DIM==3
      blade_dgx_intra_rhs_advection_3_(N,sub,rate_prefix, field_prefix,t);
#else
#error "Bad Dimension"
#endif
      break;
  }
#undef X
}

void intra_rhs (bfam_subdomain_t *thisSubdomain, const char *rate_prefix,
    const char *minus_rate_prefix, const char *field_prefix,
    const bfam_long_real_t t)
{
  BFAM_ASSERT(bfam_subdomain_has_tag(thisSubdomain,"_subdomain_dgx"));

  bfam_subdomain_dgx_t *sub = (bfam_subdomain_dgx_t*) thisSubdomain;
  if(bfam_subdomain_has_tag(thisSubdomain,"_volume"))
  {
    intra_rhs_advection(sub->N,sub,rate_prefix,field_prefix,t);
  }
  else if(bfam_subdomain_has_tag(thisSubdomain,"_glue_boundary")
      ||  bfam_subdomain_has_tag(thisSubdomain,"_glue_parallel")
      ||  bfam_subdomain_has_tag(thisSubdomain,"_glue_local"   ));
  else
    BFAM_ABORT("Unknown subdomain: %s",thisSubdomain->name);
}

void inter_rhs (bfam_subdomain_t *thisSubdomain, const char *rate_prefix,
    const char *minus_rate_prefix, const char *field_prefix,
    const bfam_long_real_t t)
{
  BFAM_ASSERT(bfam_subdomain_has_tag(thisSubdomain,"_subdomain_dgx"));

  if(bfam_subdomain_has_tag(thisSubdomain,"_volume"));
  else if(bfam_subdomain_has_tag(thisSubdomain,"_glue_parallel")
      ||  bfam_subdomain_has_tag(thisSubdomain,"_glue_local"))
  {
    BFAM_ABORT("inter_rhs: not implemented");
    BFAM_ASSERT(minus_rate_prefix);
    //JK bfam_subdomain_dgx_t *sub =
    //JK   (bfam_subdomain_dgx_t*) thisSubdomain;
    //JK inter_rhs_interface(((bfam_subdomain_dgx_t*)sub->base.glue_m->sub_m)->N,
    //JK     sub,minus_rate_prefix,field_prefix,t);
  }
  else
    BFAM_ABORT("Unknown subdomain: %s",thisSubdomain->name);
}

void add_rates_advection (bfam_subdomain_dgx_t *sub,
    const char *field_prefix_lhs, const char *field_prefix_rhs,
    const char *rate_prefix, const bfam_long_real_t a)
{
#if   DIM==2
#define X(order) \
  case order: blade_dgx_add_rates_advection_2_##order(sub->N,sub, \
                  field_prefix_lhs,field_prefix_rhs,rate_prefix,a); break;
#elif DIM==3
#define X(order) \
  case order: blade_dgx_add_rates_advection_3_##order(sub->N,sub, \
                  field_prefix_lhs,field_prefix_rhs,rate_prefix,a); break;
#else
#error "bad dimension"
#endif

  switch(sub->N)
  {
    BFAM_LIST_OF_DGX_NORDERS
    default:
#if   DIM==2
      blade_dgx_add_rates_advection_2_(sub->N,sub,field_prefix_lhs,
          field_prefix_rhs,rate_prefix,a);
#elif DIM==3
      blade_dgx_add_rates_advection_3_(sub->N,sub,field_prefix_lhs,
          field_prefix_rhs,rate_prefix,a);
#else
#error "bad dimension"
#endif
      break;
  }
#undef X
}

void add_rates (bfam_subdomain_t *thisSubdomain, const char *field_prefix_lhs,
    const char *field_prefix_rhs, const char *rate_prefix,
    const bfam_long_real_t a)
{
  BFAM_ASSERT(bfam_subdomain_has_tag(thisSubdomain,"_subdomain_dgx"));
  if(bfam_subdomain_has_tag(thisSubdomain,"_volume"))
  {
    bfam_subdomain_dgx_t *sub = (bfam_subdomain_dgx_t*)thisSubdomain;
    add_rates_advection(sub,field_prefix_lhs,field_prefix_rhs,rate_prefix,a);
  }
  else if(bfam_subdomain_has_tag(thisSubdomain,"_glue_boundary")
      ||  bfam_subdomain_has_tag(thisSubdomain,"_glue_parallel")
      ||  bfam_subdomain_has_tag(thisSubdomain,"_glue_local"   ));
  else
    BFAM_ABORT("Unknown subdomain: %s",thisSubdomain->name);
}

void glue_rates (bfam_subdomain_t *thisSubdomain, const char *prefix)
{
  if(bfam_subdomain_has_tag(thisSubdomain,"_glue_parallel") ||
     bfam_subdomain_has_tag(thisSubdomain,"_glue_local"   ))
  {
    BFAM_LDEBUG("Adding glue rate to %s",thisSubdomain->name);
    for(int f = 0 ; comm_args_scalars[f] != NULL; f++)
    {
      char field[BFAM_BUFSIZ];
      snprintf(field,BFAM_BUFSIZ,"%s%s",prefix,comm_args_scalars[f]);
      thisSubdomain->field_minus_add(thisSubdomain,field);
      thisSubdomain->field_plus_add(thisSubdomain,field);
    }
    for(int f = 0 ; comm_args_vectors[f] != NULL; f++)
    {
      char field[BFAM_BUFSIZ];
      snprintf(field,BFAM_BUFSIZ, "%s%sn",prefix,comm_args_vectors[f]);
      thisSubdomain->field_minus_add(thisSubdomain,field);
      thisSubdomain->field_plus_add(thisSubdomain,field);

      snprintf(field,BFAM_BUFSIZ, "%s%sp1",prefix,comm_args_vectors[f]);
      thisSubdomain->field_minus_add(thisSubdomain,field);
      thisSubdomain->field_plus_add(thisSubdomain,field);

      snprintf(field,BFAM_BUFSIZ, "%s%sp2",prefix,comm_args_vectors[f]);
      thisSubdomain->field_minus_add(thisSubdomain,field);
      thisSubdomain->field_plus_add(thisSubdomain,field);

      snprintf(field,BFAM_BUFSIZ, "%s%sp3",prefix,comm_args_vectors[f]);
      thisSubdomain->field_minus_add(thisSubdomain,field);
      thisSubdomain->field_plus_add(thisSubdomain,field);
    }
    for(int f = 0 ; comm_args_tensors[f] != NULL; f++)
    {
      char field[BFAM_BUFSIZ];
      snprintf(field,BFAM_BUFSIZ, "%s%sn",prefix,comm_args_tensors[f]);
      thisSubdomain->field_minus_add(thisSubdomain,field);
      thisSubdomain->field_plus_add(thisSubdomain,field);

      snprintf(field,BFAM_BUFSIZ, "%s%sp1",prefix,comm_args_tensors[f]);
      thisSubdomain->field_minus_add(thisSubdomain,field);
      thisSubdomain->field_plus_add(thisSubdomain,field);

      snprintf(field,BFAM_BUFSIZ, "%s%sp2",prefix,comm_args_tensors[f]);
      thisSubdomain->field_minus_add(thisSubdomain,field);
      thisSubdomain->field_plus_add(thisSubdomain,field);

      snprintf(field,BFAM_BUFSIZ, "%s%sp3",prefix,comm_args_tensors[f]);
      thisSubdomain->field_minus_add(thisSubdomain,field);
      thisSubdomain->field_plus_add(thisSubdomain,field);
    }
  }
}

void add_rates_glue_p (bfam_subdomain_t *thisSubdomain,
    const char *field_prefix_lhs, const char *field_prefix_rhs,
    const char *rate_prefix, const bfam_long_real_t a)
{
  BFAM_ABORT("add_rates_glue_p: not implemented");
  BFAM_ASSERT(bfam_subdomain_has_tag(thisSubdomain,"_subdomain_dgx"));
  if(bfam_subdomain_has_tag(thisSubdomain,"_glue_parallel")
      ||  bfam_subdomain_has_tag(thisSubdomain,"_glue_local"   ))
  {
    //JK bfam_subdomain_dgx_t *sub = (bfam_subdomain_dgx_t*)thisSubdomain;
    //JK add_rates_advection_glue_p(sub,field_prefix_lhs,field_prefix_rhs,
    //JK                          rate_prefix,a);
  }
  else
    BFAM_ABORT("Unknown subdomain: %s",thisSubdomain->name);
}

static void
init_time_stepper(blade_t *blade, prefs_t *prefs, bfam_locidx_t num_lvl)
{
  blade->comm_args = bfam_malloc(sizeof(bfam_subdomain_comm_args_t));
  bfam_subdomain_comm_args_t *args = blade->comm_args;

  args->scalars_m           = comm_args_scalars;
  args->vectors_m           = comm_args_vectors;
  args->vector_components_m = comm_args_vector_components;
  args->tensors_m           = comm_args_tensors;
  args->tensor_components_m = comm_args_tensor_components;
  args->face_scalars_m      = comm_args_face_scalars;

  args->scalars_p           = comm_args_scalars;
  args->vectors_p           = comm_args_vectors;
  args->vector_components_p = comm_args_vector_components;
  args->tensors_p           = comm_args_tensors;
  args->tensor_components_p = comm_args_tensor_components;
  args->face_scalars_p      = comm_args_face_scalars;

  args->user_comm_info       = NULL;
  args->user_put_send_buffer = NULL;
  args->user_get_recv_buffer = NULL;
  args->user_data            = NULL;
  args->user_prefix_function = NULL;


  const char *timestep_tags[] = {"_volume","_glue_parallel","_glue_local",
    "_glue_boundary", NULL};
  const char *glue[]   = {"_glue_parallel", "_glue_local", NULL};

  if(prefs->lsrk_method != BFAM_TS_LSRK_NOOP)
    blade->blade_ts = (bfam_ts_t*)bfam_ts_lsrk_new(
        (bfam_domain_t*) blade->domain, prefs->lsrk_method,
        BFAM_DOMAIN_OR,timestep_tags, BFAM_DOMAIN_OR,glue, blade->mpicomm, 10,
        blade->comm_args, &aux_rates,&scale_rates,&intra_rhs,&inter_rhs,
        &add_rates);
  else if(prefs->adams_method != BFAM_TS_ADAMS_NOOP)
    blade->blade_ts = (bfam_ts_t*)bfam_ts_adams_new(
        (bfam_domain_t*) blade->domain, prefs->adams_method, BFAM_DOMAIN_OR,
        timestep_tags, BFAM_DOMAIN_OR,glue, blade->mpicomm, 10,
        blade->comm_args, &aux_rates,&scale_rates,&intra_rhs,&inter_rhs,
        &add_rates, lua_get_global_int(prefs->L, "RK_init", 1));
  else if(prefs->local_adams_method != BFAM_TS_LOCAL_ADAMS_NOOP)
    blade->blade_ts = (bfam_ts_t*)bfam_ts_local_adams_new(
        (bfam_domain_t*) blade->domain, prefs->local_adams_method, num_lvl,
        BFAM_DOMAIN_OR, timestep_tags, BFAM_DOMAIN_OR,glue, blade->mpicomm, 10,
        blade->comm_args, &aux_rates, &glue_rates, &scale_rates,&intra_rhs,
        &inter_rhs, &add_rates, &add_rates_glue_p,
        lua_get_global_int(prefs->L, "RK_init", 1));
}

static bfam_real_t
compute_domain_dt(blade_t *blade, prefs_t *prefs, const char *volume[],
    const char *glue[], int *nsteps_ptr, int *ndisp_ptr, int *noutput_ptr,
    int *nerr_ptr)
{

  /* first we get all the volume subdomains */
  bfam_subdomain_t **subdomains =
    bfam_malloc(((bfam_domain_t*)blade->domain)->numSubdomains
        *sizeof(bfam_subdomain_t**));

  bfam_locidx_t numSubdomains = 0;

  bfam_domain_get_subdomains(((bfam_domain_t*)blade->domain), BFAM_DOMAIN_OR,
      volume, ((bfam_domain_t*)blade->domain)->numSubdomains,
      subdomains, &numSubdomains);

  bfam_real_t ldt[numSubdomains];
  bfam_real_t min_ldt = INFINITY;

  for(bfam_locidx_t s = 0; s < numSubdomains; ++s)
  {
    ldt[s] = INFINITY;
    bfam_subdomain_field_init(subdomains[s], "_grid_JI", 0,
        compute_subdomain_dt, &ldt[s]);
    BFAM_INFO("local dt for %s is %"BFAM_REAL_FMTe,
        subdomains[s]->name, ldt[s]);
    min_ldt = BFAM_MIN(min_ldt,ldt[s]);
  }

  BFAM_INFO("min local dt is %"BFAM_REAL_FMTe, min_ldt);

  bfam_real_t min_global_dt = 0;
  BFAM_MPI_CHECK(MPI_Allreduce(&min_ldt,&min_global_dt,1,BFAM_REAL_MPI, MPI_MIN,
        blade->mpicomm));

  bfam_real_t dt = 0;
  int result = lua_global_function_call(prefs->L, "time_step_parameters",
      "r>riii", min_global_dt, &dt, nsteps_ptr, ndisp_ptr, noutput_ptr);
  BFAM_ABORT_IF_NOT(result == 0,
      "problem with lua call to 'time_step_parameters': "
      "should be a function that takes dt "
      "and returns dt, nsteps, ndisp, noutput");
  result = lua_global_function_call(prefs->L,"nerr","r>i",dt,nerr_ptr);
  if(*nerr_ptr > 0)
  {
      const char *err_flds[] = { "error_q",NULL};
      for(int f = 0; err_flds[f] != NULL; f++)
        bfam_domain_add_field ((bfam_domain_t*)blade->domain, BFAM_DOMAIN_OR,
            volume, err_flds[f]);
  }

  bfam_locidx_t num_time_lvl = 0;
  if(prefs->local_adams_method != BFAM_TS_LOCAL_ADAMS_NOOP)
  {
    bfam_real_t max_ldt = dt;

    bfam_long_real_t dt_fudge = (bfam_long_real_t)dt /
                                (bfam_long_real_t)min_global_dt;

    int max_time_level = lua_get_global_int(prefs->L, "max_time_level", 32);
    for(bfam_locidx_t s = 0; s < numSubdomains; ++s)
    {
      ldt[s] = dt_fudge*ldt[s];

      /* keep double time_step until the level is big enough */
      bfam_locidx_t time_level = 0;
      if(ldt[s] != INFINITY) for(;(1<<(time_level+1))*dt < ldt[s];time_level++);
      time_level = BFAM_MIN(time_level,max_time_level);
      BFAM_ASSERT(time_level >= 0);

      char tag[BFAM_BUFSIZ];
      bfam_ts_local_adams_fill_level_tag(tag,BFAM_BUFSIZ,time_level);
      bfam_subdomain_add_tag(subdomains[s],tag);

      /* set this guys real dt */
      ldt[s] = dt*(1<<time_level);

      BFAM_INFO("For %s is time level %d with local dt %"BFAM_REAL_FMTe,
          subdomains[s]->name, time_level, ldt[s]);

      max_ldt = BFAM_MAX(max_ldt,ldt[s]);

      num_time_lvl = BFAM_MAX(num_time_lvl, time_level+1);
    }
    BFAM_INFO("local number of time levels %"BFAM_LOCIDX_PRId,num_time_lvl);

    BFAM_INFO("max local dt is %"BFAM_REAL_FMTe, max_ldt);

    BFAM_MPI_CHECK(MPI_Allreduce(MPI_IN_PLACE,&num_time_lvl,1,BFAM_LOCIDX_MPI,
          MPI_MAX, blade->mpicomm));

    dt = (1<<(num_time_lvl-1))*dt;

    BFAM_INFO("number of time levels %"BFAM_LOCIDX_PRId
        " for global dt %"BFAM_REAL_FMTe, num_time_lvl, dt);

    bfam_communicator_t level_comm;

    bfam_subdomain_comm_args_t level_args;
    const char * level_NULL[]      = {NULL};
    level_args.scalars_m           = level_NULL;
    level_args.vectors_m           = level_NULL;
    level_args.vector_components_m = level_NULL;
    level_args.tensors_m           = level_NULL;
    level_args.tensor_components_m = level_NULL;
    level_args.face_scalars_m      = level_NULL;

    level_args.scalars_p           = level_NULL;
    level_args.vectors_p           = level_NULL;
    level_args.vector_components_p = level_NULL;
    level_args.tensors_p           = level_NULL;
    level_args.tensor_components_p = level_NULL;
    level_args.face_scalars_p      = level_NULL;

    level_args.user_comm_info       = time_level_comm_info;
    level_args.user_put_send_buffer = time_level_put_send_buffer;
    level_args.user_get_recv_buffer = time_level_get_recv_buffer;
    level_args.user_data            = &num_time_lvl;

    level_args.user_prefix_function = NULL;


    bfam_communicator_init(&level_comm,(bfam_domain_t*)blade->domain,
        BFAM_DOMAIN_OR,glue, blade->mpicomm,10,&level_args);
    bfam_communicator_start( &level_comm);
    bfam_communicator_finish(&level_comm);
    bfam_communicator_free(  &level_comm);
  }

  bfam_free(subdomains);

  init_time_stepper(blade, prefs, num_time_lvl);

  return dt;
}

static bfam_real_t
compute_energy(blade_t *blade, prefs_t *prefs, bfam_real_t t,
    const char *prefix)
{
  const char *tags[] = {"_volume",NULL};
  bfam_subdomain_t *subs[blade->domain->base.numSubdomains];
  bfam_locidx_t num_subs = 0;
  bfam_domain_get_subdomains((bfam_domain_t*) blade->domain,
      BFAM_DOMAIN_OR,tags,blade->domain->base.numSubdomains,
      subs,&num_subs);
  bfam_real_t energy = 0;
  bfam_real_t energy_local = 0;
  for(bfam_locidx_t s = 0; s<num_subs; s++)
  {
    bfam_subdomain_dgx_t *sub = (bfam_subdomain_dgx_t*) subs[s];
#if   DIM==2
#define X(order) \
    case order: blade_dgx_energy_2_##order(sub->N,&energy_local, \
                    sub,prefix); break;
#elif DIM==3
#define X(order) \
    case order: blade_dgx_energy_3_##order(sub->N,&energy_local, \
                    sub,prefix); break;
#else
#error "Bad Dimension"
#endif

    switch(sub->N)
    {
      BFAM_LIST_OF_DGX_NORDERS
      default:
#if   DIM==2
        blade_dgx_energy_2_(sub->N,&energy_local,sub,prefix);
#elif DIM==3
        blade_dgx_energy_3_(sub->N,&energy_local,sub,prefix);
#else
#error "Bad Dimension"
#endif
        break;
    }
#undef X
  }
  BFAM_MPI_CHECK(MPI_Reduce(&energy_local,&energy,1,BFAM_REAL_MPI,
         MPI_SUM,0,blade->mpicomm));
  if(blade->mpirank == 0)
    energy = BFAM_REAL_SQRT(energy);
  return energy;
}

typedef struct check_error_args
{
  char *field_prefix;
  lua_State *L;
} check_error_args_t;


static void
check_error(bfam_locidx_t npoints, const char *name, bfam_real_t t,
    bfam_real_t *restrict x, bfam_real_t *restrict y, bfam_real_t *restrict z,
    struct bfam_subdomain *s, void *arg, bfam_real_t *restrict err)
{
  BFAM_ASSUME_ALIGNED(x, 32);
  BFAM_ASSUME_ALIGNED(y, 32);
  BFAM_ASSUME_ALIGNED(z, 32);
  BFAM_ASSUME_ALIGNED(err, 32);

  check_error_args_t* err_args = (check_error_args_t*)arg;
  lua_State *L = err_args->L;
  lua_getglobal(L,name+6);
  BFAM_ABORT_IF_NOT(lua_isfunction(L,-1),
      "no callback function for initial condition and error: %s",name+6);
  lua_pop(L,1);
  char fname[BFAM_BUFSIZ];
  snprintf(fname,BFAM_BUFSIZ,"%s%s",err_args->field_prefix,name+6);
  bfam_real_t *restrict fld = bfam_dictionary_get_value_ptr(&s->fields, fname);
  BFAM_ABORT_IF(fld == NULL, "field '%s' not in fields for %s",fname,s->name);

#if DIM==2
  bfam_real_t tmpz = 0;
#endif
  for(bfam_locidx_t n=0; n < npoints; ++n)
  {
#if   DIM==2
    lua_global_function_call(L, name+6, "rrrr>r", x[n],y[n],tmpz,t,&err[n]);
#elif DIM==3
    lua_global_function_call(L, name+6, "rrrr>r", x[n],y[n],z[n],t,&err[n]);
#else
#error "Bad Dimension"
#endif
    err[n] -= fld[n];
  }
}

static void
run_simulation(blade_t *blade,prefs_t *prefs)
{
  const char *volume[] = {"_volume",NULL};
  const char *glue[]   = {"_glue_parallel", "_glue_local", NULL};

  int nsteps  = 0;
  int ndisp   = 0;
  int noutput = 0;
  int nerr = 0;

  bfam_real_t dt = compute_domain_dt(blade, prefs, volume, glue,
      &nsteps, &ndisp, &noutput, &nerr);

  BFAM_ROOT_INFO("dt        = %"BFAM_REAL_FMTe,dt);
  BFAM_ROOT_INFO("nsteps    = %d",nsteps);
  BFAM_ROOT_INFO("ndisp     = %d",ndisp);
  BFAM_ROOT_INFO("noutput   = %d",noutput);
  BFAM_ROOT_INFO("nerr      = %d",nerr);

  /* compute the initial energy */
  bfam_real_t initial_energy = compute_energy(blade,prefs,0,"");
  bfam_real_t energy = initial_energy;
  if(initial_energy < BFAM_REAL_EPS) initial_energy = -1;
  BFAM_ROOT_INFO("\x1B[%dm"
      "time: %10.5"BFAM_REAL_PRIe
      " energy: %10.5"BFAM_REAL_PRIe
      "\x1B[0m",
      34,
      0*dt,
      energy);

  if(noutput >= 0)
  {
   char output[BFAM_BUFSIZ];

   /* dump the pxest mesh */
   snprintf(output,BFAM_BUFSIZ,"%s/%s_pxest_mesh",
       prefs->data_directory,prefs->output_prefix);
   p4est_vtk_write_all(blade->domain->pxest, NULL,
                       1, 1, 1, 1, 0, 0, 0, output);

    const char *fields[] = {"q","ux","uy",
#if DIM==3
      "uz",
#endif
      NULL};
    snprintf(output,BFAM_BUFSIZ,"%s_%05d",prefs->output_prefix,0);
    bfam_vtk_write_file((bfam_domain_t*) blade->domain, BFAM_DOMAIN_OR,
        volume, prefs->data_directory, output, (0)*dt, fields, NULL, NULL,
        prefs->vtk_binary, prefs->vtk_compress, 0);
  }

  BFAM_ASSERT(blade->blade_ts);
  for(int s = 1; s <= nsteps; s++)
  {
    blade->blade_ts->step(blade->blade_ts,dt);
    if(s%ndisp == 0)
    {
      bfam_real_t new_energy = compute_energy(blade,prefs,s*dt,"");
      if(initial_energy < 0)
      {
        BFAM_ROOT_INFO("\x1B[%dm"
            "time: %10.5"BFAM_REAL_PRIe
            "\x1B[0m",
            34,
            s*dt);
      }
      else
      {
        int color = 32;
        if(new_energy > energy) color = 31;
        BFAM_ROOT_INFO("\x1B[%dm"
            "time: %"BFAM_REAL_FMTe"\n"
            " init energy: %"BFAM_REAL_FMTe
            " energy: %"BFAM_REAL_FMTe
            " norm energy: %"BFAM_REAL_FMTe"\n"
            " current delta energy: %+"BFAM_REAL_FMTe
            " initial delta energy: %+"BFAM_REAL_FMTe"\n"
            "\x1B[0m",
            color,
            s*dt,
            initial_energy,
            new_energy,
            new_energy/initial_energy,
            (new_energy-energy)/initial_energy,
            energy/initial_energy-1);
      }
      energy = new_energy;
    }
    if(noutput > 0 && s%noutput == 0)
    {
      const char *fields[] = {"q",NULL};
      char output[BFAM_BUFSIZ];
      snprintf(output,BFAM_BUFSIZ,"%s_%05d",prefs->output_prefix,s);
      bfam_vtk_write_file((bfam_domain_t*) blade->domain, BFAM_DOMAIN_OR,
          volume, prefs->data_directory, output, (s)*dt, fields, NULL, NULL,
          prefs->vtk_binary, prefs->vtk_compress, 0);
    }
    if(nerr > 0 && s%nerr == 0)
    {
      check_error_args_t err_args;
      err_args.L = prefs->L;
      char prefix[] = "";
      err_args.field_prefix = prefix;
      const char *err_flds[] = { "error_q", NULL};
      for(int f = 0; err_flds[f] != NULL; f++)
        bfam_domain_init_field((bfam_domain_t*)blade->domain, BFAM_DOMAIN_OR,
            volume, err_flds[f], s*dt, check_error, &err_args);
      bfam_real_t error = compute_energy(blade,prefs,s*dt,"error_");
      bfam_real_t new_energy = compute_energy(blade,prefs,s*dt,"");
      BFAM_ROOT_INFO(
          "time: %"BFAM_REAL_FMTe" error: %"BFAM_REAL_FMTe
          " d_energy: %"BFAM_REAL_FMTe,
          s*dt, error,(new_energy-energy)/initial_energy);
      if(noutput > 0)
      {
        char err_output[BFAM_BUFSIZ];
        snprintf(err_output,BFAM_BUFSIZ,"%s_error_%05d",prefs->output_prefix,s);
        bfam_vtk_write_file((bfam_domain_t*) blade->domain, BFAM_DOMAIN_OR,
            volume, prefs->data_directory, err_output, (s)*dt, err_flds, NULL,
            NULL, prefs->vtk_binary, prefs->vtk_compress, 0);
        energy = new_energy;
      }
    }
  }
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

  run_simulation(&blade, prefs);

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
