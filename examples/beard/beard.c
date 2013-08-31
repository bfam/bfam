#include <bfam.h>
#include "beard_dgx_rhs.h"
#include <p4est_iterate.h>

static int refine_level     = 0;
static int max_refine_level = 0;
static bfam_real_t energy_sq = 0;
static bfam_real_t init_energy_sq = 1;

/*
const char *comm_args_scalars[]           = {NULL};
const char *comm_args_vectors[]           = {"v",NULL};
const char *comm_args_vector_components[] = {"v1","v2","v3",NULL};
const char *comm_args_tensors[]           = {"T",NULL};
const char *comm_args_tensor_components[] = {"S11","S22","S33",
                                             "S12","S13","S23",NULL};
*/
const char *comm_args_scalars[]           = {"v1","v2","v3","S11","S22","S33",
                                             "S12","S13","S23",NULL};
const char *comm_args_vectors[]           = {NULL};
const char *comm_args_vector_components[] = {NULL};
const char *comm_args_tensors[]           = {NULL};
const char *comm_args_tensor_components[] = {NULL};

static void
beard_grid_glue(bfam_locidx_t npoints, const char *name, bfam_real_t time,
    bfam_real_t *restrict x, bfam_real_t *restrict y, bfam_real_t *restrict z,
    struct bfam_subdomain *s, void *arg, bfam_real_t *restrict field)
{
  bfam_subdomain_dgx_quad_glue_t* sub_g = (bfam_subdomain_dgx_quad_glue_t*) s;

  /* get the fields we will need */
  bfam_dictionary_t *fields_m    = &sub_g->base.fields_m;
  bfam_dictionary_t *fields_p    = &sub_g->base.fields_p;
  bfam_real_t *x_m = bfam_dictionary_get_value_ptr(fields_m, "_grid_x");
  bfam_real_t *x_p = bfam_dictionary_get_value_ptr(fields_p, "_grid_x");

  bfam_real_t *y_m = bfam_dictionary_get_value_ptr(fields_m, "_grid_y");
  bfam_real_t *y_p = bfam_dictionary_get_value_ptr(fields_p, "_grid_y");

  bfam_real_t *z_m = bfam_dictionary_get_value_ptr(fields_m, "_grid_z");
  bfam_real_t *z_p = bfam_dictionary_get_value_ptr(fields_p, "_grid_z");

  // BFAM_ASSERT((sub_g->N_m+1)*sub_g->K == npoints);
  // BFAM_ASSERT((sub_g->N_p+1)*sub_g->K == npoints);
  for(int n = 0;n < npoints;n++)
  {
    // BFAM_ASSERT(BFAM_REAL_APPROX_EQ(x_m[n], x_p[n], 10));
    // BFAM_ASSERT(BFAM_REAL_APPROX_EQ(y_m[n], y_p[n], 10));
    // BFAM_ASSERT(BFAM_REAL_APPROX_EQ(z_m[n], z_p[n], 10));
    x[n] = 0.5*(x_m[n]+x_p[n]);
    y[n] = 0.5*(y_m[n]+y_p[n]);
    z[n] = 0.5*(z_m[n]+z_p[n]);
  }

}

static int
get_global_int(lua_State *L, const char *name, int def, int warning)
{
  lua_getglobal(L, name);
  int result = def;
  if(!lua_isnumber(L, -1))
  { /* for some reason gcc combines the two if without these. Why? */
    if(warning) BFAM_ROOT_WARNING("`%s' not found, using default", name);
  }
  else
    result = (int)lua_tonumber(L, -1);
  lua_pop(L, 1);
  return result;
}

static bfam_real_t
get_global_real(lua_State *L, const char *name, bfam_real_t def, int warning)
{
  lua_getglobal(L, name);
  bfam_real_t result = def;
  if(!lua_isnumber(L, -1))
  { /* for some reason gcc combines the two if without these. Why?? */
    if(warning) BFAM_ROOT_WARNING("`%s' not found, using default", name);
  }
  else
    result = (bfam_real_t)lua_tonumber(L, -1);
  lua_pop(L, 1);
  return result;
}


/*
 * Uniform refinement function
 */
static int
refine_near_fault_fn(p4est_t * p4est, p4est_topidx_t which_tree,
                     p4est_quadrant_t * quadrant)
{
  /* if not at base level refine */
  if((int)quadrant->level >= max_refine_level) return 0;
  if((int)quadrant->level < refine_level) return 1;

  /* loop through corners to see if we are near the fault */
  for(int ix = 0; ix < 2; ix++)
    for(int iy = 0; iy < 2; iy++)
    {
      int ox = ix*(1 << (P4EST_MAXLEVEL-quadrant->level));
      int oy = iy*(1 << (P4EST_MAXLEVEL-quadrant->level));
      double vxyz[3];
      p4est_qcoord_to_vertex (p4est->connectivity, which_tree,
          quadrant->x+ox, quadrant->y+oy,vxyz);
      if(BFAM_REAL_ABS(vxyz[1])<1*(max_refine_level-quadrant->level)
          && BFAM_REAL_ABS(vxyz[0])<17+1*(max_refine_level-quadrant->level))
          return 1;
    }

  return 0;
}

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

typedef struct brick_args
{
  int mi;
  int ni;
  int periodic_a;
  int periodic_b;
  bfam_real_t Lx;
  bfam_real_t Ly;
  bfam_real_t rotate;
  int random;
} brick_args_t;

typedef struct prefs
{
  lua_State *L;

  int N;
  int N_fault;
  bfam_locidx_t num_subdomains;
  bfam_real_t dt_scale;

  p4est_connectivity_t * (*conn_fn) (void);

  bfam_real_t rho;
  bfam_real_t mu;
  bfam_real_t lam;

  bfam_ts_lsrk_method_t lsrk_method;

  brick_args_t *brick;
} prefs_t;

typedef struct beard
{
  MPI_Comm mpicomm;
  int      mpirank;
  int      mpisize;

  p4est_connectivity_t *conn;
  bfam_domain_p4est_t  *domain;
  bfam_ts_lsrk_t       *lsrk;
  bfam_subdomain_comm_args_t * comm_args;
} beard_t;

static void
init_mpi(beard_t *beard, MPI_Comm mpicomm)
{
  beard->mpicomm = mpicomm;
  BFAM_MPI_CHECK(MPI_Comm_rank(mpicomm, &beard->mpirank));
  BFAM_MPI_CHECK(MPI_Comm_size(mpicomm, &beard->mpisize));
}

typedef struct split_iter_data
{
  int loc;
  bfam_locidx_t *subdomain_id;
  bfam_real_t dt;
  int N_fault;
  int base_N;
} split_iter_data_t;

static void
split_domain_treeid_iter_volume(p4est_iter_volume_info_t *info, void *arg)
{
  bfam_locidx_t num_trees = info->p4est->connectivity->num_trees;
  split_iter_data_t *data = (split_iter_data_t*) arg;
  int N = BFAM_MIN(data->base_N,
      ceil(data->N_fault*sqrt(pow(2.0,max_refine_level-info->quad->level))));
  data->subdomain_id[data->loc] = info->treeid + num_trees*(N-1);
  data->loc++;
}

static void
split_domain_treeid(beard_t *beard, int base_N, int N_fault)
{
  BFAM_ROOT_INFO("Splitting p4est based on tree id");

  bfam_domain_p4est_t *domain = beard->domain;
  if(domain->p4est->local_num_quadrants == 0) BFAM_ABORT("No quadrants");
  bfam_locidx_t *subdomain_id =
    bfam_malloc(domain->p4est->local_num_quadrants*sizeof(bfam_locidx_t));
  bfam_locidx_t num_trees = domain->p4est->connectivity->num_trees;
  int Nmax =
    BFAM_MIN(base_N,ceil(N_fault*sqrt(pow(2.0,max_refine_level-refine_level))));
  bfam_real_t num_subdomains = num_trees*Nmax;
  bfam_locidx_t *N = bfam_malloc(num_subdomains*sizeof(int));
  for(int l = 0; l < Nmax;l++)
    for(int t = 0;t < num_trees;t++)
    {
      N[l*num_trees+t] = l+1;
    }

  split_iter_data_t data = {0,subdomain_id,1,N_fault,base_N};
  p4est_iterate (domain->p4est, NULL, &data,
      split_domain_treeid_iter_volume, NULL, NULL);

  bfam_domain_p4est_split_dgx_quad_subdomains(domain, num_subdomains,
      subdomain_id, N);

  bfam_free(subdomain_id);
  bfam_free(N);
}

static void
split_domain_arbitrary(beard_t *beard, int base_N, bfam_locidx_t num_subdomains)
{
  bfam_domain_p4est_t *domain = beard->domain;
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
    N[id] = base_N+id;/*+id;*/

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
    domain->p4est->global_first_quadrant[beard->mpirank];

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

typedef struct stress_free_box_params
{
  bfam_real_t A;
  bfam_real_t rho;
  bfam_real_t mu;
  int n_ap;
  int m_ap;
} stress_free_box_params_t;

static void
stress_free_box(bfam_locidx_t npoints, const char* name, bfam_real_t t,
    bfam_real_t *restrict x, bfam_real_t *restrict y, bfam_real_t *restrict z,
    struct bfam_subdomain *s, void *arg, bfam_real_t *restrict field)
{
  stress_free_box_params_t *params = (stress_free_box_params_t *) arg;
  bfam_real_t A = params->A;
  bfam_real_t rho = params->rho;
  bfam_real_t mu = params->mu;
  int n_ap = params->n_ap;
  int m_ap = params->m_ap;
  bfam_long_real_t pi = 4*atanl(1);

  /*
  for(bfam_locidx_t n=0; n < npoints; ++n)
    field[n] = 1;
  return;
  */
  if(strcmp(name,"v3")==0)
  {
    bfam_long_real_t kx = n_ap*pi;
    bfam_long_real_t ky = m_ap*pi;
    bfam_long_real_t w = -BFAM_REAL_SQRT((mu/rho)*(kx*kx+ky*ky));
    for(bfam_locidx_t n=0; n < npoints; ++n)
      field[n] = -(A*w/mu)*sin(kx*x[n])*sin(ky*y[n])*sin(w*t);
  }
  else if(strcmp(name,"S13")==0)
  {
    bfam_long_real_t kx = n_ap*pi;
    bfam_long_real_t ky = m_ap*pi;
    bfam_long_real_t w = -BFAM_REAL_SQRT((mu/rho)*(kx*kx+ky*ky));
    for(bfam_locidx_t n=0; n < npoints; ++n)
      field[n] = A*ky*sin(kx*x[n])*cos(ky*y[n])*cos(w*t);
  }
  else if(strcmp(name,"S23")==0)
  {
    bfam_long_real_t kx = n_ap*pi;
    bfam_long_real_t ky = m_ap*pi;
    bfam_long_real_t w = -BFAM_REAL_SQRT((mu/rho)*(kx*kx+ky*ky));
    for(bfam_locidx_t n=0; n < npoints; ++n)
      field[n] = A*kx*cos(kx*x[n])*sin(ky*y[n])*cos(w*t);
  }
  else
  {
    /* BFAM_ABORT("no stress_free_box for field %s",name);*/
    for(bfam_locidx_t n=0; n < npoints; ++n)
      field[n] = 100*exp(-(pow(x[n],2) + pow(y[n],2))/2.5);
  }
}

static void
print_field(bfam_locidx_t npoints, const char *name, bfam_real_t time,
    bfam_real_t *restrict x, bfam_real_t *restrict y, bfam_real_t *restrict z,
    struct bfam_subdomain *s, void *arg, bfam_real_t *restrict field)
{
  BFAM_ASSUME_ALIGNED(x, 32);
  BFAM_ASSUME_ALIGNED(y, 32);
  BFAM_ASSUME_ALIGNED(z, 32);
  BFAM_ASSUME_ALIGNED(field, 32);

  for(bfam_locidx_t n=0; n < npoints; ++n)
    BFAM_INFO("%s: %f",name,field[n]);
}

static void
field_set_val(bfam_locidx_t npoints, const char *name, bfam_real_t time,
    bfam_real_t *restrict x, bfam_real_t *restrict y, bfam_real_t *restrict z,
    struct bfam_subdomain *s, void *arg, bfam_real_t *restrict field)
{
  BFAM_ASSUME_ALIGNED(x, 32);
  BFAM_ASSUME_ALIGNED(y, 32);
  BFAM_ASSUME_ALIGNED(z, 32);
  BFAM_ASSUME_ALIGNED(field, 32);
  bfam_real_t val = 0;
  if(arg != NULL) val = *((bfam_real_t*)arg);

  for(bfam_locidx_t n=0; n < npoints; ++n)
    field[n] = val;
}

static void
compute_dt(bfam_locidx_t npoints, const char *name, bfam_real_t time,
    bfam_real_t *restrict x, bfam_real_t *restrict y, bfam_real_t *restrict z,
    struct bfam_subdomain *s, void *arg, bfam_real_t *restrict JI)
{
  BFAM_ASSUME_ALIGNED(x, 32);
  BFAM_ASSUME_ALIGNED(y, 32);
  BFAM_ASSUME_ALIGNED(z, 32);
  BFAM_ASSUME_ALIGNED(JI, 32);
  bfam_real_t *dt = (bfam_real_t*)arg;

  bfam_real_t *restrict Jrx =
    bfam_dictionary_get_value_ptr(&s->fields, "_grid_Jrx");
  BFAM_ASSUME_ALIGNED(Jrx, 32);

  bfam_real_t *restrict Jry =
    bfam_dictionary_get_value_ptr(&s->fields, "_grid_Jry");
  BFAM_ASSUME_ALIGNED(Jry, 32);

  bfam_real_t *restrict Jsx =
    bfam_dictionary_get_value_ptr(&s->fields, "_grid_Jsx");
  BFAM_ASSUME_ALIGNED(Jsx, 32);

  bfam_real_t *restrict Jsy =
    bfam_dictionary_get_value_ptr(&s->fields, "_grid_Jsy");
  BFAM_ASSUME_ALIGNED(Jsy, 32);

  bfam_real_t *restrict mu =
    bfam_dictionary_get_value_ptr(&s->fields, "mu");
  BFAM_ASSUME_ALIGNED(mu, 32);

  bfam_real_t *restrict lam =
    bfam_dictionary_get_value_ptr(&s->fields, "lam");
  BFAM_ASSUME_ALIGNED(lam, 32);

  bfam_real_t *restrict rho =
    bfam_dictionary_get_value_ptr(&s->fields, "rho");
  BFAM_ASSUME_ALIGNED(rho, 32);

  bfam_subdomain_dgx_quad_t *sub = (bfam_subdomain_dgx_quad_t *) s;
  bfam_real_t p2 = sub->N*sub->N;

  for(int n = 0; n < npoints; n++)
  {
    bfam_real_t cp = BFAM_REAL_SQRT((lam[n]+2*mu[n])/rho[n]);
    bfam_real_t hs =
      1.0/BFAM_REAL_SQRT(JI[n]*JI[n]*(Jsy[n]*Jsy[n]+Jsx[n]*Jsx[n]));
    bfam_real_t hr =
      1.0/BFAM_REAL_SQRT(JI[n]*JI[n]*(Jry[n]*Jry[n]+Jrx[n]*Jrx[n]));
    dt[0] = BFAM_MIN(dt[0], BFAM_MIN(hs,hr)/cp/p2);
  }
}


static void
compute_dt_grid(bfam_locidx_t npoints, const char *name, bfam_real_t time,
    bfam_real_t *restrict x, bfam_real_t *restrict y, bfam_real_t *restrict z,
    struct bfam_subdomain *s, void *arg, bfam_real_t *restrict JI)
{
  BFAM_ASSUME_ALIGNED(x, 32);
  BFAM_ASSUME_ALIGNED(y, 32);
  BFAM_ASSUME_ALIGNED(z, 32);
  bfam_real_t *dt = (bfam_real_t*)arg;

  bfam_real_t *restrict mu =
    bfam_dictionary_get_value_ptr(&s->fields, "mu");
  BFAM_ASSUME_ALIGNED(mu, 32);

  bfam_real_t *restrict lam =
    bfam_dictionary_get_value_ptr(&s->fields, "lam");
  BFAM_ASSUME_ALIGNED(lam, 32);

  bfam_real_t *restrict rho =
    bfam_dictionary_get_value_ptr(&s->fields, "rho");
  BFAM_ASSUME_ALIGNED(rho, 32);

  bfam_subdomain_dgx_quad_t *sub = (bfam_subdomain_dgx_quad_t *) s;

  bfam_locidx_t Nrp = sub->N+1;
  for(bfam_locidx_t k = 0; k < sub->K; k++)
  {
    bfam_locidx_t n = k*sub->Np;
    bfam_real_t cp = BFAM_REAL_SQRT((lam[n]+2*mu[n])/rho[n]);
    for(bfam_locidx_t j = 0; j < Nrp; j++)
      for(bfam_locidx_t i = 0; i < Nrp; i++)
      {
        bfam_real_t h =  INFINITY;
        if(i < sub->N)
        {
          bfam_real_t hx = x[n+j*Nrp+i+1]-x[n+j*Nrp+i];
          bfam_real_t hy = y[n+j*Nrp+i+1]-y[n+j*Nrp+i];
          h = BFAM_MIN(h,BFAM_REAL_SQRT(hx*hx+hy*hy));
        }
        if(j < sub->N)
        {
          bfam_real_t hx = x[n+(j+1)*Nrp+i]-x[n+j*Nrp+i];
          bfam_real_t hy = y[n+(j+1)*Nrp+i]-y[n+j*Nrp+i];
          h = BFAM_MIN(h,BFAM_REAL_SQRT(hx*hx+hy*hy));
        }
        dt[0] = BFAM_MIN(dt[0], h/cp);
      }
  }
}

typedef struct field_set_val_lua_args
{
  lua_State *L;
  bfam_real_t def;
  int warning;
  const char *lua_name;
} field_set_val_lua_args_t;

static void
lua_field_set_val_lua (lua_State *L, bfam_locidx_t npoints, const char *name,
    bfam_real_t time, const bfam_real_t *x, const bfam_real_t *y,
    const bfam_real_t *z, const bfam_real_t def, int warning,
    bfam_real_t *field)
{
  lua_getglobal(L, name);
  bfam_real_t val = def;
  int lua_func = 0;

  if(lua_isnumber(L,-1))
    val = (bfam_real_t)lua_tonumber(L, -1);
  else if(lua_isfunction(L,-1))
    lua_func = 1;
  else if(warning) BFAM_ROOT_WARNING("`%s' not found, using default", name);
  lua_pop(L,1);

  for(int n = 0; n < npoints; n++)
  {
    if(lua_func)
    {
      lua_getglobal(L, name);
      lua_pushnumber(L, time);
      lua_pushnumber(L, x[n]);
      lua_pushnumber(L, y[n]);
      lua_pushnumber(L, z[n]);
      if(lua_pcall(L,4,1,0) != 0)
        BFAM_ABORT("error running function %s: %s",name,lua_tostring(L,-1));
      if(!lua_isnumber(L,-1))
        BFAM_ABORT("function '%s' must return number",name);
      val = (bfam_real_t)lua_tonumber(L,-1);
      lua_pop(L,1);
    }
    field[n] = val;
  }
}

static void
field_set_val_lua (bfam_locidx_t npoints, const char *name, bfam_real_t time,
    bfam_real_t *restrict x, bfam_real_t *restrict y, bfam_real_t *restrict z,
    struct bfam_subdomain *s, void *args_in, bfam_real_t *restrict field)
{
  BFAM_ASSERT(args_in != NULL);
  field_set_val_lua_args_t *args = (field_set_val_lua_args_t *)args_in;
  lua_field_set_val_lua (args->L, npoints, args->lua_name, time, x, y, z,
      args->def, args->warning, field);
}

static void
field_set_val_normal_lua(bfam_locidx_t npoints, const char *name,
    bfam_real_t time, bfam_real_t *restrict x, bfam_real_t *restrict y,
    bfam_real_t *restrict z, struct bfam_subdomain *s, void *args_in,
    bfam_real_t *restrict field)
{
  BFAM_ASSERT(args_in != NULL);
  field_set_val_lua_args_t *args = (field_set_val_lua_args_t *)args_in;

  bfam_subdomain_dgx_quad_glue_t *sub_g = (bfam_subdomain_dgx_quad_glue_t*)s;

  bfam_dictionary_t *fields_face = &sub_g->sub_m->base.fields_face;
  bfam_real_t *n1 = bfam_dictionary_get_value_ptr(fields_face, "_grid_nx");
  bfam_real_t *n2 = bfam_dictionary_get_value_ptr(fields_face, "_grid_ny");

  // BFAM_ASSERT((sub_g->N_m+1)*sub_g->K == npoints);
  // BFAM_ASSERT((sub_g->N  +1)*sub_g->K == npoints);
  // BFAM_ASSERT((sub_g->N_p+1)*sub_g->K == npoints);

  lua_State *L = args->L;
  lua_getglobal(L, args->lua_name);

  if(!lua_isfunction(L,-1))
    BFAM_ABORT("`%s' must be a lua function", args->lua_name);
  lua_pop(L,1);

  for(int le = 0; le<sub_g->K; le++)
  {
    bfam_locidx_t e = sub_g->EToEm[le];
    int8_t face = sub_g->EToFm[le];
    /* WARNING ASSUMES STRAIGHT SIDED ELEMENTS!!!*/
    bfam_real_t nm[] = {n1[sub_g->sub_m->Nfp*(face+4*e)],
                        n2[sub_g->sub_m->Nfp*(face+4*e)],0};

    for(int pnt = 0; pnt < sub_g->N+1;pnt++)
    {
      bfam_locidx_t n = pnt+le*(sub_g->N+1);
      lua_getglobal (L, args->lua_name);
      lua_pushnumber(L, time);
      lua_pushnumber(L, x[n]);
      lua_pushnumber(L, y[n]);
      lua_pushnumber(L, z[n]);
      lua_pushnumber(L, nm[0]);
      lua_pushnumber(L, nm[1]);
      lua_pushnumber(L, nm[2]);
      if(lua_pcall(L,7,1,0) != 0)
        BFAM_ABORT("error running function %s: %s",
            args->lua_name,lua_tostring(L,-1));
      if(!lua_isnumber(L,-1))
        BFAM_ABORT("function '%s' must return number",args->lua_name);
      field[n] = (bfam_real_t)lua_tonumber(L,-1);
      lua_pop(L,1);
    }
  }
}


void aux_rates (bfam_subdomain_t *thisSubdomain, const char *prefix)
{
  if(bfam_subdomain_has_tag(thisSubdomain,"_volume"))
  {
    const char *fields[] =
      {"v1","v2","v3","S11","S22","S33","S12","S13","S23",NULL};

    char field[BFAM_BUFSIZ];
    for(int f = 0; fields[f]!=NULL; ++f)
    {
      snprintf(field,BFAM_BUFSIZ,"%s%s",prefix,fields[f]);
      thisSubdomain->field_add(thisSubdomain,field);
      bfam_subdomain_field_init(thisSubdomain, field, 0, field_set_val, NULL);
    }
  }
  else if(bfam_subdomain_has_tag(thisSubdomain,"slip weakening"))
  {
    const char *fields[] = {"Dp","Dn",NULL};

    char field[BFAM_BUFSIZ];
    for(int f = 0; fields[f]!=NULL; ++f)
    {
      snprintf(field,BFAM_BUFSIZ,"%s%s",prefix,fields[f]);
      thisSubdomain->field_add(thisSubdomain,field);
      bfam_subdomain_field_init(thisSubdomain, field, 0, field_set_val, NULL);
    }
  }
}

void scale_rates_elastic (bfam_subdomain_dgx_quad_t *sub,
    const char *rate_prefix, const bfam_long_real_t a)
{
#define X(order) \
  case order: beard_dgx_scale_rates_elastic_##order(sub->N,sub, \
                  rate_prefix,a); break;

  switch(sub->N)
  {
    BFAM_LIST_OF_DGX_QUAD_NORDERS
    default:
      beard_dgx_scale_rates_elastic_(sub->N,sub,rate_prefix,a);
      break;
  }
#undef X
}


void scale_rates_slip_weakening (bfam_subdomain_dgx_quad_glue_t *sub,
    const char *rate_prefix, const bfam_long_real_t a)
{
#define X(order) \
  case order: beard_dgx_scale_rates_slip_weakening_##order(sub->N,sub, \
                  rate_prefix,a); break;

  switch(sub->N)
  {
    BFAM_LIST_OF_DGX_QUAD_NORDERS
    default:
      beard_dgx_scale_rates_slip_weakening_(sub->N,sub,rate_prefix,a);
      break;
  }
#undef X
}

void scale_rates (bfam_subdomain_t *thisSubdomain, const char *rate_prefix,
    const bfam_long_real_t a)
{
  BFAM_ASSERT(bfam_subdomain_has_tag(thisSubdomain,"_subdomain_dgx_quad"));
  if(bfam_subdomain_has_tag(thisSubdomain,"_volume"))
  {
    bfam_subdomain_dgx_quad_t *sub = (bfam_subdomain_dgx_quad_t*)thisSubdomain;
    scale_rates_elastic(sub,rate_prefix,a);
  }
  else if(bfam_subdomain_has_tag(thisSubdomain,"slip weakening"))
  {
    bfam_subdomain_dgx_quad_glue_t *sub =
      (bfam_subdomain_dgx_quad_glue_t*)thisSubdomain;
    scale_rates_slip_weakening(sub,rate_prefix,a);
  }
  else if(bfam_subdomain_has_tag(thisSubdomain,"_glue_boundary"));
  else if(bfam_subdomain_has_tag(thisSubdomain,"_glue_parallel"));
  else if(bfam_subdomain_has_tag(thisSubdomain,"_glue_local"));
  else
    BFAM_ABORT("Uknown subdomain: %s",thisSubdomain->name);
}

static void
intra_rhs_elastic(int N, bfam_subdomain_dgx_quad_t *sub,
    const char *rate_prefix, const char *field_prefix, const bfam_long_real_t t)
{
#define X(order) \
  case order: beard_dgx_intra_rhs_elastic_##order(N,sub, \
                  rate_prefix,field_prefix,t); break;

  switch(N)
  {
    BFAM_LIST_OF_DGX_QUAD_NORDERS
    default:
      beard_dgx_intra_rhs_elastic_(N,sub,rate_prefix,
          field_prefix,t);
      break;
  }
#undef X
}

void intra_rhs (bfam_subdomain_t *thisSubdomain, const char *rate_prefix,
    const char *field_prefix, const bfam_long_real_t t)
{
  BFAM_ASSERT(bfam_subdomain_has_tag(thisSubdomain,"_subdomain_dgx_quad"));

  bfam_subdomain_dgx_quad_t *sub = (bfam_subdomain_dgx_quad_t*) thisSubdomain;
  if(bfam_subdomain_has_tag(thisSubdomain,"_volume"))
    intra_rhs_elastic(sub->N,sub,rate_prefix,field_prefix,t);
  else if(bfam_subdomain_has_tag(thisSubdomain,"_glue_boundary")
      ||  bfam_subdomain_has_tag(thisSubdomain,"_glue_parallel")
      ||  bfam_subdomain_has_tag(thisSubdomain,"_glue_local"   ));
  else
    BFAM_ABORT("Uknown subdomain: %s",thisSubdomain->name);
}

void inter_rhs_boundary(int N, bfam_subdomain_dgx_quad_glue_t *sub,
    const char *rate_prefix, const char *field_prefix, const bfam_long_real_t t,
    const bfam_real_t R)
{
#define X(order) \
  case order: beard_dgx_inter_rhs_boundary_##order(N,sub, \
                  rate_prefix,field_prefix,t, R); break;

  switch(N)
  {
    BFAM_LIST_OF_DGX_QUAD_NORDERS
    default:
      beard_dgx_inter_rhs_boundary_(N,sub,rate_prefix, field_prefix,t, R);
      break;
  }
#undef X
}

void inter_rhs_interface(int N, bfam_subdomain_dgx_quad_glue_t *sub,
    const char *rate_prefix, const char *field_prefix, const bfam_long_real_t t)
{
#define X(order) \
  case order: beard_dgx_inter_rhs_interface_##order(N,sub, \
                  rate_prefix,field_prefix,t); break;

  switch(N)
  {
    BFAM_LIST_OF_DGX_QUAD_NORDERS
    default:
      beard_dgx_inter_rhs_interface_(N,sub,rate_prefix,
          field_prefix,t);
      break;
  }
#undef X
}

void inter_rhs_slip_weakening_interface(int N, bfam_subdomain_dgx_quad_glue_t
    *sub, const char *rate_prefix, const char *field_prefix,
    const bfam_long_real_t t)
{
#define X(order) \
  case order: beard_dgx_inter_rhs_slip_weakening_interface_##order(N,sub, \
                  rate_prefix,field_prefix,t); break;

  switch(N)
  {
    BFAM_LIST_OF_DGX_QUAD_NORDERS
    default:
      beard_dgx_inter_rhs_slip_weakening_interface_(N,sub,rate_prefix,
          field_prefix,t);
      break;
  }
#undef X
}

void inter_rhs (bfam_subdomain_t *thisSubdomain, const char *rate_prefix,
    const char *field_prefix, const bfam_long_real_t t)
{
  BFAM_ASSERT(bfam_subdomain_has_tag(thisSubdomain,"_subdomain_dgx_quad"));

  bfam_subdomain_dgx_quad_glue_t *sub =
    (bfam_subdomain_dgx_quad_glue_t*) thisSubdomain;
  if(bfam_subdomain_has_tag(thisSubdomain,"_volume"));
  else if(bfam_subdomain_has_tag(thisSubdomain,"_glue_boundary"))
    inter_rhs_boundary(sub->sub_m->N,sub,rate_prefix,field_prefix,t,0);
  else if(bfam_subdomain_has_tag(thisSubdomain,"slip weakening"))
    inter_rhs_slip_weakening_interface(sub->sub_m->N,sub,rate_prefix,
        field_prefix,t);
  else if(bfam_subdomain_has_tag(thisSubdomain,"_glue_parallel")
      ||  bfam_subdomain_has_tag(thisSubdomain,"_glue_local"))
    inter_rhs_interface(sub->sub_m->N,sub,rate_prefix,field_prefix,t);
  else
    BFAM_ABORT("Uknown subdomain: %s",thisSubdomain->name);
}

void add_rates_slip_weakening (bfam_subdomain_dgx_quad_glue_t *sub,
    const char *field_prefix_lhs, const char *field_prefix_rhs,
    const char *rate_prefix, const bfam_long_real_t a)
{
#define X(order) \
  case order: beard_dgx_add_rates_slip_weakening_##order(sub->N,sub, \
                  field_prefix_lhs,field_prefix_rhs,rate_prefix,a); break;

  switch(sub->N)
  {
    BFAM_LIST_OF_DGX_QUAD_NORDERS
    default:
      beard_dgx_add_rates_slip_weakening_(sub->N,sub,field_prefix_lhs,
          field_prefix_rhs,rate_prefix,a);
      break;
  }
#undef X
}

void add_rates_elastic (bfam_subdomain_dgx_quad_t *sub,
    const char *field_prefix_lhs, const char *field_prefix_rhs,
    const char *rate_prefix, const bfam_long_real_t a)
{
#define X(order) \
  case order: beard_dgx_add_rates_elastic_##order(sub->N,sub, \
                  field_prefix_lhs,field_prefix_rhs,rate_prefix,a); break;

  switch(sub->N)
  {
    BFAM_LIST_OF_DGX_QUAD_NORDERS
    default:
      beard_dgx_add_rates_elastic_(sub->N,sub,field_prefix_lhs,
          field_prefix_rhs,rate_prefix,a);
      break;
  }
#undef X
}

void add_rates (bfam_subdomain_t *thisSubdomain, const char *field_prefix_lhs,
    const char *field_prefix_rhs, const char *rate_prefix,
    const bfam_long_real_t a)
{
  BFAM_ASSERT(bfam_subdomain_has_tag(thisSubdomain,"_subdomain_dgx_quad"));
  if(bfam_subdomain_has_tag(thisSubdomain,"_volume"))
  {
    bfam_subdomain_dgx_quad_t *sub = (bfam_subdomain_dgx_quad_t*)thisSubdomain;
    add_rates_elastic(sub,field_prefix_lhs,field_prefix_rhs,rate_prefix,a);
  }
  else if(bfam_subdomain_has_tag(thisSubdomain,"slip weakening"))
  {
    bfam_subdomain_dgx_quad_glue_t *sub =
      (bfam_subdomain_dgx_quad_glue_t*)thisSubdomain;
    add_rates_slip_weakening(sub,field_prefix_lhs,field_prefix_rhs,rate_prefix,a);
  }
  else if(bfam_subdomain_has_tag(thisSubdomain,"_glue_boundary")
      ||  bfam_subdomain_has_tag(thisSubdomain,"_glue_parallel")
      ||  bfam_subdomain_has_tag(thisSubdomain,"_glue_local"   ));
  else
    BFAM_ABORT("Uknown subdomain: %s",thisSubdomain->name);
}

static void
init_lsrk(beard_t *beard, prefs_t *prefs)
{
  beard->comm_args = bfam_malloc(sizeof(bfam_subdomain_comm_args_t));
  bfam_subdomain_comm_args_t *args = beard->comm_args;

  args->scalars_m           = comm_args_scalars;
  args->vectors_m           = comm_args_vectors;
  args->vector_components_m = comm_args_vector_components;
  args->tensors_m           = comm_args_tensors;
  args->tensor_components_m = comm_args_tensor_components;

  args->scalars_p           = comm_args_scalars;
  args->vectors_p           = comm_args_vectors;
  args->vector_components_p = comm_args_vector_components;
  args->tensors_p           = comm_args_tensors;
  args->tensor_components_p = comm_args_tensor_components;


  const char *timestep_tags[] = {"_volume","_glue_parallel","_glue_local",
    "_glue_boundary", NULL};
  const char *glue[]   = {"_glue_parallel", "_glue_local", NULL};

  beard->lsrk = bfam_ts_lsrk_new((bfam_domain_t*) beard->domain,
      prefs->lsrk_method, BFAM_DOMAIN_OR,timestep_tags, BFAM_DOMAIN_OR,glue,
      beard->mpicomm, 10, beard->comm_args,
      &aux_rates,&scale_rates,&intra_rhs,&inter_rhs, &add_rates);
}

static void
compute_energy(beard_t *beard, prefs_t *prefs, bfam_real_t t, int init_energy)
{
  const char *tags[] = {"_volume",NULL};
  bfam_subdomain_t *subs[beard->domain->base.numSubdomains];
  bfam_locidx_t num_subs = 0;
  bfam_domain_get_subdomains((bfam_domain_t*) beard->domain,
      BFAM_DOMAIN_OR,tags,beard->domain->base.numSubdomains,
      subs,&num_subs);
  bfam_real_t energy_sq_old = energy_sq;
  bfam_real_t energy_sq_local = 0;
  for(bfam_locidx_t s = 0; s<num_subs; s++)
  {
    bfam_subdomain_dgx_quad_t *sub = (bfam_subdomain_dgx_quad_t*) subs[s];
#define X(order) \
    case order: beard_dgx_energy_##order(sub->N,&energy_sq_local, \
                    sub,""); break;

    switch(sub->N)
    {
      BFAM_LIST_OF_DGX_QUAD_NORDERS
      default:
        beard_dgx_energy_(sub->N,&energy_sq_local,sub,"");
        break;
    }
#undef X
  }
  BFAM_MPI_CHECK(MPI_Reduce(&energy_sq_local,&energy_sq,1,BFAM_REAL_MPI,
         MPI_SUM,0,beard->mpicomm));

  if(!init_energy)
    if(energy_sq > energy_sq_old)
      BFAM_ROOT_INFO("\x1B[31m"
          "time: %f energy_sq: %e D_energy_sq: %+e"
          " NDI_energy: %+e"
          "\x1B[0m",
          t,
          BFAM_REAL_SQRT(energy_sq/init_energy_sq),
          BFAM_REAL_SQRT(energy_sq/init_energy_sq)
          - BFAM_REAL_SQRT(energy_sq_old/init_energy_sq),
          BFAM_REAL_SQRT(energy_sq/init_energy_sq)-1);
    else
      BFAM_ROOT_INFO("\x1B[32m"
          "time: %f energy_sq: %e D_energy_sq: %+e"
          " NDI_energy: %+e"
          "\x1B[0m",
          t,
          BFAM_REAL_SQRT(energy_sq/init_energy_sq),
          BFAM_REAL_SQRT(energy_sq/init_energy_sq)
          - BFAM_REAL_SQRT(energy_sq_old/init_energy_sq),
          BFAM_REAL_SQRT(energy_sq/init_energy_sq)-1);
  else
    init_energy_sq = energy_sq;
}

static void
init_domain(beard_t *beard, prefs_t *prefs)
{
  if(prefs->brick != NULL)
  {
    /* just so we always get the same random grid for convergence tests */
    srandom(8);
    beard->conn = p4est_connectivity_new_brick(
        prefs->brick->mi,prefs->brick->ni,
        prefs->brick->periodic_a,prefs->brick->periodic_b);
    for(int i = 0; i < beard->conn->num_vertices; i++)
    {
      int x = beard->conn->vertices[i*3+0];
      int y = beard->conn->vertices[i*3+1];
      if(x > 0 && x < prefs->brick->mi && y > 0 && y < prefs->brick->ni
          && prefs->brick->random)
      {
        bfam_real_t r = 2*(random() / (bfam_real_t) RAND_MAX - 0.5);
        beard->conn->vertices[i*3+0] = x + 0.7*0.5*r;
      }
      beard->conn->vertices[i*3+0] -= 0.5*prefs->brick->mi;
      beard->conn->vertices[i*3+0] *= 2*prefs->brick->Lx / prefs->brick->mi;

      if(x > 0 && x < prefs->brick->mi && y > 0 && y < prefs->brick->ni
          && prefs->brick->random)
      {
        bfam_real_t r = 2*(random() / (bfam_real_t) RAND_MAX - 0.5);
        beard->conn->vertices[i*3+1] = y + 0.7*0.5*r;
      }
      beard->conn->vertices[i*3+1] -= 0.5*prefs->brick->ni;
      beard->conn->vertices[i*3+1] *= 2*prefs->brick->Ly / prefs->brick->ni;

      if(prefs->brick->rotate > 0)
      {
        double q = prefs->brick->rotate;
        bfam_real_t x = beard->conn->vertices[i*3+0];
        bfam_real_t y = beard->conn->vertices[i*3+1];
        beard->conn->vertices[i*3+0] = x*cos(q) + y*sin(q);
        beard->conn->vertices[i*3+1] =-x*sin(q) + y*cos(q);
      }

      // BFAM_INFO("%e %e %e",
      //     beard->conn->vertices[i*3+0],
      //     beard->conn->vertices[i*3+1],
      //     beard->conn->vertices[i*3+2]);
    }
  }
  else
    beard->conn = prefs->conn_fn();

  beard->domain = bfam_domain_p4est_new(beard->mpicomm, beard->conn);

  // p4est_refine(beard->domain->p4est, 2, refine_fn, NULL);
  p4est_refine(beard->domain->p4est, 2, refine_near_fault_fn, NULL);
  p4est_balance(beard->domain->p4est, P4EST_CONNECT_CORNER, NULL);
  p4est_partition(beard->domain->p4est, NULL);

  p4est_vtk_write_file(beard->domain->p4est, NULL, "p4est_mesh");

  // split_domain_arbitrary(beard, prefs->N, prefs->num_subdomains);
  split_domain_treeid(beard,prefs->N,prefs->N_fault);

  const char *volume[] = {"_volume", NULL};
  const char *glue[]   = {"_glue_parallel", "_glue_local", NULL};

  bfam_domain_t *domain = (bfam_domain_t*)beard->domain;

  bfam_domain_add_field(domain, BFAM_DOMAIN_OR, volume, "rho");
  bfam_domain_add_field(domain, BFAM_DOMAIN_OR, volume, "rho_inv");
  bfam_domain_add_field(domain, BFAM_DOMAIN_OR, volume, "lam");
  bfam_domain_add_field(domain, BFAM_DOMAIN_OR, volume, "mu");
  bfam_domain_add_field(domain, BFAM_DOMAIN_OR, volume, "Zs");
  bfam_domain_add_field(domain, BFAM_DOMAIN_OR, volume, "Zp");

  /* set material properties */
  bfam_domain_init_field(domain, BFAM_DOMAIN_OR, volume, "rho", 0,
      field_set_val, &prefs->rho);

  bfam_real_t rho_inv = 1.0/prefs->rho;
  bfam_domain_init_field(domain, BFAM_DOMAIN_OR, volume, "rho_inv", 0,
      field_set_val, &rho_inv);

  bfam_domain_init_field(domain, BFAM_DOMAIN_OR, volume, "mu", 0,
      field_set_val, &prefs->mu);

  bfam_domain_init_field(domain, BFAM_DOMAIN_OR, volume, "lam", 0,
      field_set_val, &prefs->lam);

  bfam_real_t Zs = sqrt(prefs->rho*prefs->mu);
  bfam_domain_init_field(domain, BFAM_DOMAIN_OR, volume, "Zs", 0,
      field_set_val, &Zs);

  bfam_real_t Zp = sqrt(prefs->rho*(prefs->lam+2*prefs->mu));
  bfam_domain_init_field(domain, BFAM_DOMAIN_OR, volume, "Zp", 0,
      field_set_val, &Zp);


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

  /* zero out fields */
  bfam_domain_init_field(domain, BFAM_DOMAIN_OR, volume, "v1", 0,
      field_set_val, NULL);
  bfam_domain_init_field(domain, BFAM_DOMAIN_OR, volume, "v2", 0,
      field_set_val, NULL);
  bfam_domain_init_field(domain, BFAM_DOMAIN_OR, volume, "v3", 0,
      field_set_val, NULL);
  bfam_domain_init_field(domain, BFAM_DOMAIN_OR, volume, "S11", 0,
      field_set_val, NULL);
  bfam_domain_init_field(domain, BFAM_DOMAIN_OR, volume, "S22", 0,
      field_set_val, NULL);
  bfam_domain_init_field(domain, BFAM_DOMAIN_OR, volume, "S33", 0,
      field_set_val, NULL);
  bfam_domain_init_field(domain, BFAM_DOMAIN_OR, volume, "S12", 0,
      field_set_val, NULL);
  bfam_domain_init_field(domain, BFAM_DOMAIN_OR, volume, "S13", 0,
      field_set_val, NULL);
  bfam_domain_init_field(domain, BFAM_DOMAIN_OR, volume, "S23", 0,
      field_set_val, NULL);


  /* exchange material properties to glue */
  const char *glue_mat[] = {"Zs","Zp","_grid_x","_grid_y","_grid_z",NULL};
  for(bfam_locidx_t g = 0; glue_mat[g] != NULL; g++)
  {
    bfam_domain_add_minus_field(domain, BFAM_DOMAIN_OR, glue, glue_mat[g]);
    bfam_domain_add_plus_field( domain, BFAM_DOMAIN_OR, glue, glue_mat[g]);
  }
  bfam_domain_add_field( domain, BFAM_DOMAIN_OR, glue, "_grid_x");
  bfam_domain_add_field( domain, BFAM_DOMAIN_OR, glue, "_grid_y");
  bfam_domain_add_field( domain, BFAM_DOMAIN_OR, glue, "_grid_z");
  bfam_communicator_t material_comm;

  bfam_subdomain_comm_args_t mat_args;
  const char * mat_NULL[]      = {NULL};
  mat_args.scalars_m           = glue_mat;
  mat_args.vectors_m           = mat_NULL;
  mat_args.vector_components_m = mat_NULL;
  mat_args.tensors_m           = mat_NULL;
  mat_args.tensor_components_m = mat_NULL;

  mat_args.scalars_p           = glue_mat;
  mat_args.vectors_p           = mat_NULL;
  mat_args.vector_components_p = mat_NULL;
  mat_args.tensors_p           = mat_NULL;
  mat_args.tensor_components_p = mat_NULL;

  bfam_communicator_init(&material_comm,domain,BFAM_DOMAIN_OR,glue,
      beard->mpicomm,10,&mat_args);
  bfam_communicator_start( &material_comm);
  bfam_communicator_finish(&material_comm);
  bfam_communicator_free(  &material_comm);

  /*
   * we can trick init fields into handling locations
   * to glue
   */
  bfam_domain_init_field(domain, BFAM_DOMAIN_OR, glue, "_grid_x", 0,
      beard_grid_glue, NULL);

  /* Set up some problem specific problems */
  lua_getglobal(prefs->L, "problem");
  if(lua_isstring(prefs->L, -1))
  {
    const char *prob_name = lua_tostring(prefs->L, -1);

    if(strcmp(prob_name,"stress free box") == 0)
    {

      stress_free_box_params_t field_params;
      field_params.A    = 1;
      field_params.rho  = prefs->rho;
      field_params.mu   = prefs->mu;
      field_params.n_ap = 2;
      field_params.m_ap = 2;

      bfam_domain_init_field(domain, BFAM_DOMAIN_OR, volume, "v3", 0,
          stress_free_box, &field_params);
      bfam_domain_init_field(domain, BFAM_DOMAIN_OR, volume, "S13", 0,
          stress_free_box, &field_params);
      bfam_domain_init_field(domain, BFAM_DOMAIN_OR, volume, "S23", 0,
          stress_free_box, &field_params);
      /*
      bfam_domain_init_field(domain, BFAM_DOMAIN_OR, volume, "v1", 0,
          stress_free_box, &field_params);
      bfam_domain_init_field(domain, BFAM_DOMAIN_OR, volume, "v2", 0,
          stress_free_box, &field_params);
      */


    }
    else if(strcmp(prob_name,"slip weakening") == 0)
    {
      bfam_subdomain_t * fault_sub[domain->numSubdomains];
      bfam_locidx_t num_subs;
      char friction_name[2][BFAM_BUFSIZ];
      const char *friction_tag[] = {friction_name[0],friction_name[1],NULL};
      const int num_trees = beard->conn->num_trees;
      snprintf(friction_name[0],BFAM_BUFSIZ,"_glue_%d_%d",
          num_trees*(prefs->N_fault-1),2+num_trees*(prefs->N_fault-1));
      snprintf(friction_name[1],BFAM_BUFSIZ,"_glue_%d_%d",
          1+num_trees*(prefs->N_fault-1),3+num_trees*(prefs->N_fault-1));
      bfam_domain_get_subdomains(domain, BFAM_DOMAIN_OR, friction_tag,
          domain->numSubdomains, fault_sub, &num_subs);
      for(int s = 0; s < num_subs; s++)
        bfam_subdomain_add_tag(fault_sub[s],"slip weakening");

      const char *fric[]   = {"slip weakening", NULL};

      const char *fault_fields[] = {"Dp", "Dn", "V", "Tp1_0", "Tp2_0", "Tp3_0",
        "Tn_0", "Tp1", "Tp2", "Tp3", "Tn",  "fs",  "fd", "Dc",NULL};
      for(int fld = 0; fault_fields[fld] != NULL; fld++)
      {
        bfam_domain_add_field( domain, BFAM_DOMAIN_OR, fric, fault_fields[fld]);
        bfam_domain_init_field(domain, BFAM_DOMAIN_OR, fric, fault_fields[fld],
            0, field_set_val, NULL);
      }

      const char * fault_lua_field[] = {"fs", "fd", "Dc", NULL};
      const char * fault_lua_names[] = {"fault_fs", "fault_fd",
        "fault_Dc", NULL};
      for(int fld = 0; fault_lua_field[fld] != NULL; fld++)
      {
        BFAM_ASSERT(fault_lua_names[fld] != NULL);
        field_set_val_lua_args_t args;
        args.L = prefs->L;
        args.def = 0;
        args.warning = 1;
        args.lua_name = fault_lua_names[fld];
        bfam_domain_init_field(domain, BFAM_DOMAIN_OR, fric,
            fault_lua_field[fld], 0, field_set_val_lua, &args);
      }

      const char * fault_lua_field_normal[] = {"Tp1_0", "Tp2_0", "Tp3_0",
        "Tn_0", NULL};
      const char * fault_lua_names_normal[] = {"fault_Tp1_0", "fault_Tp2_0",
        "fault_Tp3_0", "fault_Tn_0", NULL};
      for(int fld = 0; fault_lua_field_normal[fld] != NULL; fld++)
      {
        BFAM_ASSERT(fault_lua_names_normal[fld] != NULL);
        field_set_val_lua_args_t args;
        args.L = prefs->L;
        args.def = 0;
        args.warning = 1;
        args.lua_name = fault_lua_names_normal[fld];
        bfam_domain_init_field(domain, BFAM_DOMAIN_OR, fric,
            fault_lua_field_normal[fld], 0, field_set_val_normal_lua, &args);
      }
    }
  }
  lua_pop(prefs->L, 1);
}

static void
shave_beard(beard_t *beard)
{
  bfam_free(beard->comm_args);
  bfam_ts_lsrk_free(beard->lsrk);
  bfam_free(beard->lsrk);
  bfam_domain_p4est_free(beard->domain);
  bfam_free(beard->domain);
  p4est_connectivity_destroy(beard->conn);
}

static void
run(MPI_Comm mpicomm, prefs_t *prefs)
{
  const char *fields[] = {"v1", "v2", "v3",
    "S11", "S22", "S33", "S12", "S13", "S23", NULL};
  // const char *fields[] = {"v1", NULL};
  const char *volume[] = {"_volume", NULL};
  const char *fault_fields[] = {"Dp", "Dn", "V",NULL};
  const char *fric_tags[]   = {"slip weakening", NULL};

  beard_t beard;

  init_mpi(&beard, mpicomm);

  init_domain(&beard, prefs);

  init_lsrk(&beard, prefs);

  const char directory[] = "output";
  char output[BFAM_BUFSIZ];
  snprintf(output,BFAM_BUFSIZ,"solution_%05d",0);
  bfam_vtk_write_file((bfam_domain_t*) beard.domain, BFAM_DOMAIN_OR, volume,
      directory,output,0, fields, NULL, NULL, 1, 1);

  snprintf(output,BFAM_BUFSIZ,"solution_friction_%05d",0);
  bfam_vtk_write_file((bfam_domain_t*) beard.domain, BFAM_DOMAIN_OR, fric_tags,
      directory,output,0, fault_fields, NULL, NULL, 0, 0);

  bfam_real_t ldt = INFINITY;
  bfam_domain_init_field((bfam_domain_t*) beard.domain, BFAM_DOMAIN_OR, volume,
      "_grid_JI", 0, compute_dt, &ldt);
  /*
  bfam_domain_init_field((bfam_domain_t*) beard.domain, BFAM_DOMAIN_OR, volume,
      "_grid_JI", 0, compute_dt_grid, &ldt);
  */
  ldt  *= prefs->dt_scale;
  bfam_real_t dt = 0;
  BFAM_MPI_CHECK(MPI_Allreduce(&ldt,&dt,1,BFAM_REAL_MPI, MPI_MIN,beard.mpicomm));

  BFAM_INFO("global dt = %e (local dt = %e)", dt, ldt);

  int nsteps = get_global_int(prefs->L,"nsteps",1000,0);
  int ndisp  = get_global_int(prefs->L,"ndisp" ,  10,0);
  int noutput_body  = get_global_int(prefs->L,"noutput_body",   100,0);
  int noutput_fault  = get_global_int(prefs->L,"noutput_fault",  10,0);

  lua_getglobal(prefs->L, "dt_modify");
  if(lua_isfunction(prefs->L,-1))
  {
    lua_pushnumber(prefs->L, dt);
    if(lua_pcall(prefs->L,1,5,0) != 0)
      BFAM_ABORT("error running function %s: %s",
          "dt_modify",lua_tostring(prefs->L,-1));

    dt = (bfam_real_t)lua_tonumber(prefs->L,-1);
    lua_pop(prefs->L,1);

    nsteps = (int)lua_tonumber(prefs->L,-1);
    lua_pop(prefs->L,1);

    ndisp = (int)lua_tonumber(prefs->L,-1);
    lua_pop(prefs->L,1);

    noutput_body = (int)lua_tonumber(prefs->L,-1);
    lua_pop(prefs->L,1);

    noutput_fault = (int)lua_tonumber(prefs->L,-1);
    lua_pop(prefs->L,1);
  }
  else lua_pop(prefs->L,1);

  BFAM_ROOT_INFO("dt:     %"BFAM_REAL_FMTe,dt);
  BFAM_ROOT_INFO("nsteps: %04d",nsteps);
  BFAM_ROOT_INFO("ndisp:  %04d",ndisp);
  BFAM_ROOT_INFO("nout_b: %04d",noutput_body);
  BFAM_ROOT_INFO("nout_f: %04d",noutput_fault);

  compute_energy(&beard,prefs,0,0);

  for(int s = 1; s <= nsteps; s++)
  {
    beard.lsrk->base.step((bfam_ts_t*) beard.lsrk,dt);
    if(s%noutput_body == 0)
    {
      snprintf(output,BFAM_BUFSIZ,"solution_%05d",s);
      bfam_vtk_write_file((bfam_domain_t*) beard.domain, BFAM_DOMAIN_OR, volume,
          directory,output,(s)*dt, fields, NULL, NULL, 1, 1);
    }
    if(s%noutput_fault == 0)
    {
      snprintf(output,BFAM_BUFSIZ,"solution_friction_%05d",s);
      bfam_vtk_write_file((bfam_domain_t*) beard.domain, BFAM_DOMAIN_OR, fric_tags,
          directory,output,(s)*dt, fault_fields, NULL, NULL, 0, 0);
    }
    if(s%ndisp == 0)
    {
      compute_energy(&beard,prefs,(s)*dt,0);
    }
  }

  shave_beard(&beard);
}

static void
print_order(int N)
{
#define X(order) \
  case order: beard_dgx_print_order_##order(N); break;

  switch(N)
  {
    BFAM_LIST_OF_DGX_QUAD_NORDERS
    default:
      beard_dgx_print_order_(N);
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

  prefs->N = get_global_int(L, "N", 5, 1);
  prefs->N_fault = get_global_int(L, "N_fault", 1, 1);
  prefs->dt_scale = get_global_real(L, "dt_scale", 1, 1);
  prefs->num_subdomains = get_global_int(L, "num_subdomains", 1, 1);
  refine_level = get_global_int(L, "refine_level", 0, 1);
  max_refine_level = get_global_int(L, "max_refine_level", refine_level, 1);
  prefs->rho = get_global_real(L, "rho", 1, 1);
  prefs->mu  = get_global_real(L, "mu" , 1, 1);
  prefs->lam = get_global_real(L, "lam", 1, 1);
  BFAM_ASSERT(lua_gettop(L)==0);

  lua_getglobal(L, "connectivity");
  prefs->conn_fn = conn_table[0].conn_fn;
  prefs->brick = NULL;
  if(lua_isstring(L, -1))
  {
    int i;
    const char *conn_name = lua_tostring(L, -1);

    if(strcmp(conn_name,"brick") == 0)
    {
      prefs->conn_fn = NULL;

      prefs->brick = bfam_malloc(sizeof(brick_args_t));
      prefs->brick->ni = get_global_int(L,"brick_n",1,1);
      prefs->brick->mi = get_global_int(L,"brick_m",1,1);
      prefs->brick->periodic_a = get_global_int(L,"brick_a",0,1);
      prefs->brick->periodic_b = get_global_int(L,"brick_b",0,1);
      prefs->brick->random = get_global_int(L,"brick_random",0,1);
      prefs->brick->rotate = get_global_int(L,"brick_rotate",0,1);
      prefs->brick->Lx = get_global_real(L,"brick_Lx",1,1);
      prefs->brick->Ly = get_global_real(L,"brick_Ly",1,0);
      BFAM_ROOT_INFO("%e",(double)prefs->brick->Lx);
      BFAM_ROOT_INFO("%e",(double)prefs->brick->Ly);
    }
    else
    {
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
  }
  else
  {
    BFAM_ROOT_WARNING("`connectivity' not found, using default");
  }
  BFAM_ASSERT(prefs->conn_fn != NULL || prefs->brick != NULL);
  lua_pop(L, 1);

  lua_getglobal(L, "lsrk_method");
  prefs->lsrk_method = lsrk_table[0].lsrk_method;
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
      BFAM_LERROR("invalid lsrk method name: `%s'; using default", lsrk_name);
    else
    {
      prefs->lsrk_method = lsrk_table[i].lsrk_method;
    }
  }
  else
  {
    BFAM_ROOT_WARNING("`lsrk method' not found, using default");
  }
  BFAM_ASSERT(prefs->lsrk_method != BFAM_TS_LSRK_NOOP);
  lua_pop(L, 1);


  BFAM_ASSERT(lua_gettop(L)==0);

  return prefs;
}

static void
print_prefs(prefs_t *prefs)
{
  BFAM_ROOT_INFO("----------Preferences----------");
  BFAM_ROOT_INFO("N=%d", prefs->N);
  BFAM_ROOT_INFO("N_fault=%d", prefs->N_fault);
  BFAM_ROOT_INFO("dt_scale=%f", prefs->dt_scale);
  BFAM_ROOT_INFO("num_subdomains=%"BFAM_LOCIDX_PRId, prefs->num_subdomains);
  for(int i=0; conn_table[i].name!=NULL; ++i)
    if(conn_table[i].conn_fn == prefs->conn_fn)
      BFAM_ROOT_INFO("conn_fn=`%s'", conn_table[i].name);
  for(int i=0; lsrk_table[i].name!=NULL; ++i)
    if(lsrk_table[i].lsrk_method == prefs->lsrk_method)
      BFAM_ROOT_INFO("lsrk_method=`%s'", lsrk_table[i].name);
  BFAM_ROOT_INFO("rho=%f", prefs->rho);
  BFAM_ROOT_INFO("lam=%f", prefs->lam);
  BFAM_ROOT_INFO("mu =%f", prefs->mu );
  BFAM_ROOT_INFO("-------------------------------");
}

static void
free_prefs(prefs_t *prefs)
{
  if(prefs->brick != NULL) bfam_free(prefs->brick);
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
