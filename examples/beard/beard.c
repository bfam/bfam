#include <bfam.h>
#include "beard_dgx_rhs.h"

static int refine_level = 0;
static bfam_real_t energy_sq = 0;

/*
const char *comm_args_scalars[]           = {NULL};
const char *comm_args_vectors[]           = {"v",NULL};
const char *comm_args_vector_components[] = {"v1","v2","v3",NULL};
const char *comm_args_tensors[]           = {"T",NULL};
const char *comm_args_tensor_components[] = {"S11","S22","S33",
                                             "S12","S13","S23",NULL};
*/
const char *comm_args_scalars[]           = {"v1","v2","v3","S11","S22","S33",
                                             "S12","S13","S23","Zs","Zp",NULL};
const char *comm_args_vectors[]           = {NULL};
const char *comm_args_vector_components[] = {NULL};
const char *comm_args_tensors[]           = {NULL};
const char *comm_args_tensor_components[] = {NULL};

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

typedef struct brick_args
{
  int mi;
  int ni;
  int periodic_a;
  int periodic_b;
} brick_args_t;

typedef struct prefs
{
  lua_State *L;

  int N;
  bfam_locidx_t num_subdomains;

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
    N[id] = base_N;/*+id;*/

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
    bfam_long_real_t w = -sqrt((mu/rho)*(kx*kx+ky*ky));
    for(bfam_locidx_t n=0; n < npoints; ++n)
      field[n] = -(A*w/mu)*sin(kx*x[n])*sin(ky*y[n])*sin(w*t);
  }
  else if(strcmp(name,"S13")==0)
  {
    bfam_long_real_t kx = n_ap*pi;
    bfam_long_real_t ky = m_ap*pi;
    bfam_long_real_t w = -sqrt((mu/rho)*(kx*kx+ky*ky));
    for(bfam_locidx_t n=0; n < npoints; ++n)
      field[n] = A*ky*sin(kx*x[n])*cos(ky*y[n])*cos(w*t);
  }
  else if(strcmp(name,"S23")==0)
  {
    bfam_long_real_t kx = n_ap*pi;
    bfam_long_real_t ky = m_ap*pi;
    bfam_long_real_t w = -sqrt((mu/rho)*(kx*kx+ky*ky));
    for(bfam_locidx_t n=0; n < npoints; ++n)
      field[n] = A*kx*cos(kx*x[n])*sin(ky*y[n])*cos(w*t);
  }
  else
  {
    /* BFAM_ABORT("no stress_free_box for field %s",name);*/
    for(bfam_locidx_t n=0; n < npoints; ++n)
      field[n] = exp(-(pow(x[n],2) + pow(y[n],2))*50);
  }
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

void scale_rates (bfam_subdomain_t *thisSubdomain, const char *rate_prefix,
    const bfam_long_real_t a)
{
  BFAM_ASSERT(bfam_subdomain_has_tag(thisSubdomain,"_subdomain_dgx_quad"));
  bfam_subdomain_dgx_quad_t *sub = (bfam_subdomain_dgx_quad_t*)thisSubdomain;
  if(bfam_subdomain_has_tag(thisSubdomain,"_volume"))
    scale_rates_elastic(sub,rate_prefix,a);
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
    const char *rate_prefix, const char *field_prefix, const bfam_long_real_t t)
{
#define X(order) \
  case order: beard_dgx_inter_rhs_boundary_##order(N,sub, \
                  rate_prefix,field_prefix,t, 1); break;

  switch(N)
  {
    BFAM_LIST_OF_DGX_QUAD_NORDERS
    default:
      beard_dgx_inter_rhs_boundary_(N,sub,rate_prefix,
          field_prefix,t, 1);
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

void inter_rhs (bfam_subdomain_t *thisSubdomain, const char *rate_prefix,
    const char *field_prefix, const bfam_long_real_t t)
{
  BFAM_ASSERT(bfam_subdomain_has_tag(thisSubdomain,"_subdomain_dgx_quad"));

  bfam_subdomain_dgx_quad_glue_t *sub =
    (bfam_subdomain_dgx_quad_glue_t*) thisSubdomain;
  if(bfam_subdomain_has_tag(thisSubdomain,"_volume"));
  else if(bfam_subdomain_has_tag(thisSubdomain,"_glue_boundary"))
    inter_rhs_boundary(sub->N,sub,rate_prefix,field_prefix,t);
  else if(bfam_subdomain_has_tag(thisSubdomain,"_glue_parallel")
      ||  bfam_subdomain_has_tag(thisSubdomain,"_glue_local"))
    inter_rhs_interface(sub->N,sub,rate_prefix,field_prefix,t);
  else
    BFAM_ABORT("Uknown subdomain: %s",thisSubdomain->name);
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
  bfam_subdomain_dgx_quad_t *sub = (bfam_subdomain_dgx_quad_t*)thisSubdomain;
  if(bfam_subdomain_has_tag(thisSubdomain,"_volume"))
    add_rates_elastic(sub,field_prefix_lhs,field_prefix_rhs,rate_prefix,a);
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
compute_energy(beard_t *beard, prefs_t *prefs)
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

  BFAM_ROOT_INFO("energy: %e delta energy: %e", energy_sq, energy_sq-energy_sq_old);
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
      if(x > 0 && x < prefs->brick->mi)
      {
        bfam_real_t r = random() / (bfam_real_t) RAND_MAX;
        beard->conn->vertices[i*3+0] = x + (r-0.5)/2;
      }
      beard->conn->vertices[i*3+0] -= 0.5*prefs->brick->mi;
      beard->conn->vertices[i*3+0] /= 0.5*prefs->brick->mi;

      int y = beard->conn->vertices[i*3+1];
      if(y > 0 && y < prefs->brick->ni)
      {
        bfam_real_t r = random() / (bfam_real_t) RAND_MAX;
        beard->conn->vertices[i*3+1] = y + (r-0.5)/2;
      }
      beard->conn->vertices[i*3+1] -= 0.5*prefs->brick->ni;
      beard->conn->vertices[i*3+1] /= 0.5*prefs->brick->ni;

      // BFAM_INFO("%e %e %e",
      //     beard->conn->vertices[i*3+0],
      //     beard->conn->vertices[i*3+1],
      //     beard->conn->vertices[i*3+2]);
    }
  }
  else
    beard->conn = prefs->conn_fn();

  beard->domain = bfam_domain_p4est_new(beard->mpicomm, beard->conn);

  p4est_refine(beard->domain->p4est, 2, refine_fn, NULL);
  p4est_balance(beard->domain->p4est, P4EST_CONNECT_CORNER, NULL);
  p4est_partition(beard->domain->p4est, NULL);

  p4est_vtk_write_file(beard->domain->p4est, NULL, "p4est_mesh");

  split_domain_arbitrary(beard, prefs->N, prefs->num_subdomains);

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
  bfam_domain_init_field(domain, BFAM_DOMAIN_OR, volume, "v1", 0,
      stress_free_box, &field_params);
  bfam_domain_init_field(domain, BFAM_DOMAIN_OR, volume, "v2", 0,
      stress_free_box, &field_params);

  /* exchange material properties to glue */
  const char *glue_mat[] = {"Zs","Zp",NULL};
  for(bfam_locidx_t g = 0; glue_mat[g] != NULL; g++)
  {
    bfam_domain_add_minus_field(domain, BFAM_DOMAIN_OR, glue, glue_mat[g]);
    bfam_domain_add_plus_field( domain, BFAM_DOMAIN_OR, glue, glue_mat[g]);
  }

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

  beard_t beard;

  init_mpi(&beard, mpicomm);

  init_domain(&beard, prefs);

  init_lsrk(&beard, prefs);

  char output[BFAM_BUFSIZ];

  snprintf(output,BFAM_BUFSIZ,"fields_%05d",0);
  bfam_vtk_write_file((bfam_domain_t*) beard.domain, BFAM_DOMAIN_OR, volume,
      output, fields, NULL, NULL, 0, 0);
  bfam_real_t dt = 0.01/pow(2,refine_level);
  if(prefs->brick != NULL)
    dt /= (bfam_real_t) BFAM_MAX(prefs->brick->mi,prefs->brick->ni);
  int nsteps = 10/dt;
  int ndisp  = 0.1 / dt;
  // nsteps = 1/dt;
  // ndisp  = 1;
  for(int s = 0; s < nsteps; s++)
  {
    if(s%ndisp == 0)
    {
      snprintf(output,BFAM_BUFSIZ,"fields_%05d",s+1);
      bfam_vtk_write_file((bfam_domain_t*) beard.domain, BFAM_DOMAIN_OR, volume,
          output, fields, NULL, NULL, 0, 0);
      compute_energy(&beard,prefs);
    }
    beard.lsrk->base.step((bfam_ts_t*) beard.lsrk,dt);
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

static bfam_real_t
get_global_real(lua_State *L, const char *name, bfam_real_t def)
{
  lua_getglobal(L, name);
  bfam_real_t result = def;
  if(!lua_isnumber(L, -1))
    BFAM_ROOT_WARNING("`%s' not found, using default", name);
  else
    result = (bfam_real_t)lua_tonumber(L, -1);
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
  prefs->rho = get_global_real(L, "rho", 1);
  prefs->mu  = get_global_real(L, "mu" , 1);
  prefs->lam = get_global_real(L, "lam", 1);
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
      prefs->brick->ni = get_global_int(L,"brick_n",1);
      prefs->brick->mi = get_global_int(L,"brick_m",1);
      prefs->brick->periodic_a = get_global_int(L,"brick_a",0);
      prefs->brick->periodic_b = get_global_int(L,"brick_b",0);
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
