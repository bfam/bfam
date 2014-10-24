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
  p8est_partition(domain->pxest, 0, NULL);

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

  commargs.user_comm_info = NULL;
  commargs.user_get_recv_buffer = NULL;
  commargs.user_put_send_buffer = NULL;
  commargs.user_data = NULL;
  commargs.user_prefix_function = NULL;

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
                      NULL,"volume",0, NULL, NULL, NULL, 0, 0, 0);

  const char *glue_id_0[]   = {"_glue_id_0", NULL};

  bfam_vtk_write_file(d, BFAM_DOMAIN_OR, glue_id_0,
                      NULL,"glue_id",0, NULL, NULL, NULL, 0, 0, 0);

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

static void
poly1_field(bfam_locidx_t npoints, const char* name,
    bfam_real_t time, bfam_real_t *restrict x, bfam_real_t *restrict y,
    bfam_real_t *restrict z, struct bfam_subdomain *s, void *arg,
    bfam_real_t *restrict field)
{
  BFAM_ASSUME_ALIGNED(x, 32);
  BFAM_ASSUME_ALIGNED(y, 32);
  BFAM_ASSUME_ALIGNED(z, 32);
  BFAM_ASSUME_ALIGNED(field, 32);

  for(bfam_locidx_t n=0; n < npoints; ++n)
    field[n] = x[n];
}

static void
poly2_field(bfam_locidx_t npoints, const char* name,
    bfam_real_t time, bfam_real_t *restrict x, bfam_real_t *restrict y,
    bfam_real_t *restrict z, struct bfam_subdomain *s, void *arg,
    bfam_real_t *restrict field)
{
  BFAM_ASSUME_ALIGNED(x, 32);
  BFAM_ASSUME_ALIGNED(y, 32);
  BFAM_ASSUME_ALIGNED(z, 32);
  BFAM_ASSUME_ALIGNED(field, 32);

  for(bfam_locidx_t n=0; n < npoints; ++n)
    field[n] = y[n];
}

static void
poly3_field(bfam_locidx_t npoints, const char* name,
    bfam_real_t time, bfam_real_t *restrict x, bfam_real_t *restrict y,
    bfam_real_t *restrict z, struct bfam_subdomain *s, void *arg,
    bfam_real_t *restrict field)
{
  BFAM_ASSUME_ALIGNED(x, 32);
  BFAM_ASSUME_ALIGNED(y, 32);
  BFAM_ASSUME_ALIGNED(z, 32);
  BFAM_ASSUME_ALIGNED(field, 32);

  for(bfam_locidx_t n=0; n < npoints; ++n)
    field[n] = z[n];
}

static int
check_pm(bfam_subdomain_dgx_t *sub, const char *name, bfam_real_t fac)
{
  int failures = 0;
  bfam_real_t *f_m = bfam_dictionary_get_value_ptr(&sub->base.glue_m->fields,
      name);
  bfam_real_t *f_p = bfam_dictionary_get_value_ptr(&sub->base.glue_p->fields,
      name);

  BFAM_ASSERT(f_m != NULL);
  BFAM_ASSERT(f_p != NULL);

  BFAM_LDEBUG("Testing subdomain (%2jd, %2jd) -- (%2jd, %2jd) field %s",
      (intmax_t)sub->base.glue_m->rank, (intmax_t)sub->base.glue_m->id_s,
      (intmax_t)sub->base.glue_p->rank, (intmax_t)sub->base.glue_p->id_s,
      name);

  bfam_subdomain_dgx_glue_data_t* glue_p =
    (bfam_subdomain_dgx_glue_data_t*) sub->base.glue_p;

  for(bfam_locidx_t i = 0; i < sub->K; ++i)
  {
    BFAM_LDEBUG("Testing element %2jd face %d h %d o %d",
        (intmax_t)glue_p->EToEm[i], glue_p->EToFm[i], glue_p->EToHm[i],
        glue_p->EToOm[i]);
    // For(bfam_locidx_t j = 0; j < sub->Np; ++j)
    //   BFAM_LDEBUG("fm[%2d][%2d] = %20"BFAM_REAL_PRIe
    //         "    fp[%2d][%2d] = %20"BFAM_REAL_PRIe,
    //       i, j, f_m[i*sub->Np + j], i, j, fac*f_p[i*sub->Np + j]);

    for(int j=0; j<sub->Np; ++j)
    {
      size_t idx = i*sub->Np + j;
      int fail = !REAL_APPROX_EQ(f_m[idx], fac*f_p[idx], 100);

      if(fail)
      BFAM_LDEBUG("Fail Match: fm[%2d][%2d] = %20"BFAM_REAL_PRIe
            "    fp[%2d][%2d] = %20"BFAM_REAL_PRIe
            "    fp-fm = %20"BFAM_REAL_PRIe,
            i, j, f_m[i*sub->Np + j], i, j, fac*f_p[i*sub->Np + j],
            f_m[i*sub->Np + j]-fac*f_p[i*sub->Np + j]);

      failures += fail;
    }
  }

  if(failures > 0)
    BFAM_WARNING("FAIL! %s",name);
  return failures;
}

static int
check_vmaps(bfam_subdomain_dgx_t *sub, const char *name)
{
  int failures = 0;
  bfam_real_t *f = bfam_dictionary_get_value_ptr(&sub->base.fields, name);

  BFAM_ABORT_IF(f == NULL, "subdomain %s: doesn't have a field %s", sub->base.name,
      name);

  for(bfam_locidx_t i = 0; i < sub->K * sub->Ngp[0] * sub->Ng[0]; ++i)
  {
    int fail = !REAL_APPROX_EQ(f[sub->vmapM[i]], f[sub->vmapP[i]], 1000);

    if(fail)
      BFAM_LDEBUG("Fail Match fm[%2jd] = %20"BFAM_REAL_PRIe
                         "    fp[%2jd] = %20"BFAM_REAL_PRIe,
                         (intmax_t)i, f[sub->vmapM[i]],
                         (intmax_t)i, f[sub->vmapP[i]]);

    failures += fail;
  }

  if(failures > 0)
    BFAM_WARNING("FAIL! %s",name);
  return failures;
}

static int
test_conn(MPI_Comm mpicomm, state_t *state)
{
  int failures = 0;
  bfam_domain_pxest_t_3 *domain = state->domain;

  const char *volume[] = {"_volume", NULL};
  const char *glue[]   = {"_glue_parallel", "_glue_local", NULL};

  bfam_domain_add_field((bfam_domain_t*)domain, BFAM_DOMAIN_OR, volume, "p1");
  bfam_domain_add_field((bfam_domain_t*)domain, BFAM_DOMAIN_OR, volume, "p2");
  bfam_domain_add_field((bfam_domain_t*)domain, BFAM_DOMAIN_OR, volume, "p3");
  bfam_domain_add_field((bfam_domain_t*)domain, BFAM_DOMAIN_OR, volume, "p4");
  bfam_domain_add_field((bfam_domain_t*)domain, BFAM_DOMAIN_OR, volume, "p5");
  bfam_domain_add_field((bfam_domain_t*)domain, BFAM_DOMAIN_OR, volume, "p6");

  bfam_domain_add_field((bfam_domain_t*)domain, BFAM_DOMAIN_OR, glue, "p1");
  bfam_domain_add_field((bfam_domain_t*)domain, BFAM_DOMAIN_OR, glue, "p2");
  bfam_domain_add_field((bfam_domain_t*)domain, BFAM_DOMAIN_OR, glue, "p3");
  bfam_domain_add_field((bfam_domain_t*)domain, BFAM_DOMAIN_OR, glue, "p4");
  bfam_domain_add_field((bfam_domain_t*)domain, BFAM_DOMAIN_OR, glue, "p5");
  bfam_domain_add_field((bfam_domain_t*)domain, BFAM_DOMAIN_OR, glue, "p6");

  bfam_domain_init_field((bfam_domain_t*)domain, BFAM_DOMAIN_OR, volume, "p1",
      0, poly1_field, NULL);
  bfam_domain_init_field((bfam_domain_t*)domain, BFAM_DOMAIN_OR, volume, "p2",
      0, poly2_field, NULL);
  bfam_domain_init_field((bfam_domain_t*)domain, BFAM_DOMAIN_OR, volume, "p3",
      0, poly3_field, NULL);
  bfam_domain_init_field((bfam_domain_t*)domain, BFAM_DOMAIN_OR, volume, "p4",
      0, poly1_field, NULL);
  bfam_domain_init_field((bfam_domain_t*)domain, BFAM_DOMAIN_OR, volume, "p5",
      0, poly2_field, NULL);
  bfam_domain_init_field((bfam_domain_t*)domain, BFAM_DOMAIN_OR, volume, "p6",
      0, poly3_field, NULL);


  bfam_subdomain_comm_args_t commargs;

  const char *comm_args_face_scalars[]      = {"_grid_nx0", "_grid_nx1",
                                               "_grid_nx2", NULL};
  const char *comm_args_scalars[]           = {"p1","p2","p3","p4","p5","p6", NULL};
  const char *comm_args_vectors[]           = {"v","u",NULL};
  const char *comm_args_vector_components[] = {"p3","p1","p2",
                                               "p4","p5","p6",NULL};
  const char *comm_args_tensors[]           = {"T","S",NULL};
  const char *comm_args_tensor_components[] = {"p3", "p2", "p1",
                                               "p6", "p4", "p5",
                                               "p2", "p3", "p1",
                                               "p5", "p4", "p6", NULL};
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

  commargs.user_comm_info = NULL;
  commargs.user_get_recv_buffer = NULL;
  commargs.user_put_send_buffer = NULL;
  commargs.user_data = NULL;
  commargs.user_prefix_function = NULL;

  /* add glue fields */
  for(int f = 0 ; comm_args_face_scalars[f] != NULL; f++)
  {
    bfam_domain_add_minus_field((bfam_domain_t*) domain, BFAM_DOMAIN_OR, glue,
        comm_args_face_scalars[f]);
    bfam_domain_add_plus_field( (bfam_domain_t*) domain, BFAM_DOMAIN_OR, glue,
        comm_args_face_scalars[f]);
  }
  for(int f = 0 ; comm_args_scalars[f] != NULL; f++)
  {
    bfam_domain_add_minus_field((bfam_domain_t*) domain, BFAM_DOMAIN_OR, glue,
        comm_args_scalars[f]);
    bfam_domain_add_plus_field( (bfam_domain_t*) domain, BFAM_DOMAIN_OR, glue,
        comm_args_scalars[f]);
  }
  for(int f = 0 ; comm_args_vectors[f] != NULL; f++)
  {
    char name[BFAM_BUFSIZ];
    snprintf(name,BFAM_BUFSIZ, "%sn",comm_args_vectors[f]);
    bfam_domain_add_minus_field((bfam_domain_t*) domain, BFAM_DOMAIN_OR, glue,
        name);
    bfam_domain_add_plus_field( (bfam_domain_t*) domain, BFAM_DOMAIN_OR, glue,
        name);
    snprintf(name,BFAM_BUFSIZ, "%sp1",comm_args_vectors[f]);
    bfam_domain_add_minus_field((bfam_domain_t*) domain, BFAM_DOMAIN_OR, glue,
        name);
    bfam_domain_add_plus_field( (bfam_domain_t*) domain, BFAM_DOMAIN_OR, glue,
        name);
    snprintf(name,BFAM_BUFSIZ, "%sp2",comm_args_vectors[f]);
    bfam_domain_add_minus_field((bfam_domain_t*) domain, BFAM_DOMAIN_OR, glue,
        name);
    bfam_domain_add_plus_field( (bfam_domain_t*) domain, BFAM_DOMAIN_OR, glue,
        name);
    snprintf(name,BFAM_BUFSIZ, "%sp3",comm_args_vectors[f]);
    bfam_domain_add_minus_field((bfam_domain_t*) domain, BFAM_DOMAIN_OR, glue,
        name);
    bfam_domain_add_plus_field( (bfam_domain_t*) domain, BFAM_DOMAIN_OR, glue,
        name);
  }
  for(int f = 0 ; comm_args_tensors[f] != NULL; f++)
  {
    char name[BFAM_BUFSIZ];
    snprintf(name,BFAM_BUFSIZ, "%sn",comm_args_tensors[f]);
    bfam_domain_add_minus_field((bfam_domain_t*) domain, BFAM_DOMAIN_OR, glue,
        name);
    bfam_domain_add_plus_field( (bfam_domain_t*) domain, BFAM_DOMAIN_OR, glue,
        name);
    snprintf(name,BFAM_BUFSIZ, "%sp1",comm_args_tensors[f]);
    bfam_domain_add_minus_field((bfam_domain_t*) domain, BFAM_DOMAIN_OR, glue,
        name);
    bfam_domain_add_plus_field( (bfam_domain_t*) domain, BFAM_DOMAIN_OR, glue,
        name);
    snprintf(name,BFAM_BUFSIZ, "%sp2",comm_args_tensors[f]);
    bfam_domain_add_minus_field((bfam_domain_t*) domain, BFAM_DOMAIN_OR, glue,
        name);
    bfam_domain_add_plus_field( (bfam_domain_t*) domain, BFAM_DOMAIN_OR, glue,
        name);
    snprintf(name,BFAM_BUFSIZ, "%sp3",comm_args_tensors[f]);
    bfam_domain_add_minus_field((bfam_domain_t*) domain, BFAM_DOMAIN_OR, glue,
        name);
    bfam_domain_add_plus_field( (bfam_domain_t*) domain, BFAM_DOMAIN_OR, glue,
        name);
  }

  bfam_communicator_t* communicator =
    bfam_communicator_new((bfam_domain_t*)domain, BFAM_DOMAIN_OR, glue,
        mpicomm, 10, &commargs);

  /* start recv_send */
  bfam_communicator_start(communicator);

  /* finish recv */
  bfam_communicator_finish(communicator);

{
    bfam_locidx_t numElements = domain->base.numSubdomains;
    bfam_subdomain_t **subdomains =
      bfam_malloc(numElements*sizeof(bfam_subdomain_t*));
    bfam_locidx_t numSubdomains;

    bfam_domain_get_subdomains((bfam_domain_t*)domain, BFAM_DOMAIN_OR, glue,
        numElements, subdomains, &numSubdomains);

    for(bfam_locidx_t i = 0; i < numSubdomains; ++i)
    {
      bfam_subdomain_t *s = subdomains[i];
      bfam_subdomain_dgx_t *sub = (bfam_subdomain_dgx_t*) s;

      bfam_real_t *restrict p1 = bfam_dictionary_get_value_ptr(&s->fields, "p1");
      bfam_real_t *restrict p2 = bfam_dictionary_get_value_ptr(&s->fields, "p2");
      bfam_real_t *restrict p3 = bfam_dictionary_get_value_ptr(&s->fields, "p3");
      bfam_real_t *restrict p4 = bfam_dictionary_get_value_ptr(&s->fields, "p4");
      bfam_real_t *restrict p5 = bfam_dictionary_get_value_ptr(&s->fields, "p5");
      bfam_real_t *restrict p6 = bfam_dictionary_get_value_ptr(&s->fields, "p6");

      bfam_real_t *restrict p1_m = bfam_dictionary_get_value_ptr(&s->glue_m->fields, "p1");
      bfam_real_t *restrict p2_m = bfam_dictionary_get_value_ptr(&s->glue_m->fields, "p2");
      bfam_real_t *restrict p3_m = bfam_dictionary_get_value_ptr(&s->glue_m->fields, "p3");
      bfam_real_t *restrict p4_m = bfam_dictionary_get_value_ptr(&s->glue_m->fields, "p4");
      bfam_real_t *restrict p5_m = bfam_dictionary_get_value_ptr(&s->glue_m->fields, "p5");
      bfam_real_t *restrict p6_m = bfam_dictionary_get_value_ptr(&s->glue_m->fields, "p6");

      bfam_real_t *restrict p1_p = bfam_dictionary_get_value_ptr(&s->glue_p->fields, "p1");
      bfam_real_t *restrict p2_p = bfam_dictionary_get_value_ptr(&s->glue_p->fields, "p2");
      bfam_real_t *restrict p3_p = bfam_dictionary_get_value_ptr(&s->glue_p->fields, "p3");
      bfam_real_t *restrict p4_p = bfam_dictionary_get_value_ptr(&s->glue_p->fields, "p4");
      bfam_real_t *restrict p5_p = bfam_dictionary_get_value_ptr(&s->glue_p->fields, "p5");
      bfam_real_t *restrict p6_p = bfam_dictionary_get_value_ptr(&s->glue_p->fields, "p6");

      for(bfam_locidx_t j = 0; j < sub->Np * sub->K; ++j)
      {
        p1[j] = (p1_m[j] + p1_p[j])/2;
        p2[j] = (p2_m[j] + p2_p[j])/2;
        p3[j] = (p3_m[j] + p3_p[j])/2;
        p4[j] = (p4_m[j] + p4_p[j])/2;
        p5[j] = (p5_m[j] + p5_p[j])/2;
        p6[j] = (p6_m[j] + p6_p[j])/2;
      }
    }

    bfam_free(subdomains);
  }

  const char *glue_id_0[]   = {"_glue_id_0", NULL};

  const char *ps[] = {"p1", "p2", "p3", "p4", "p5", "p6", NULL};

  bfam_vtk_write_file((bfam_domain_t*)domain, BFAM_DOMAIN_OR, glue_id_0,
                      NULL,"glue_id",0, ps, NULL, NULL, 0, 0, 0);

  /*
   * Check local subdomain vmaps
   */
  {
    bfam_subdomain_t **subdomains =
      bfam_malloc(domain->base.numSubdomains*sizeof(bfam_subdomain_t**));

    bfam_locidx_t numSubdomains = 0;

    bfam_domain_get_subdomains((bfam_domain_t*)domain, BFAM_DOMAIN_OR,
        volume, domain->base.numSubdomains, subdomains, &numSubdomains);

    BFAM_LDEBUG("Number of volume subdomains %jd", (intmax_t) numSubdomains);

    for(bfam_locidx_t s = 0; s < numSubdomains; ++s)
    {
      failures +=
        check_vmaps((bfam_subdomain_dgx_t*)subdomains[s], "p1");
      failures +=
        check_vmaps((bfam_subdomain_dgx_t*)subdomains[s], "p2");
      failures +=
        check_vmaps((bfam_subdomain_dgx_t*)subdomains[s], "p3");
      failures +=
        check_vmaps((bfam_subdomain_dgx_t*)subdomains[s], "p4");
      failures +=
        check_vmaps((bfam_subdomain_dgx_t*)subdomains[s], "p5");
      failures +=
        check_vmaps((bfam_subdomain_dgx_t*)subdomains[s], "p6");
    }

    bfam_free(subdomains);
  }

  {
    bfam_subdomain_t **subdomains =
      bfam_malloc(domain->base.numSubdomains*sizeof(bfam_subdomain_t**));

    bfam_locidx_t numSubdomains = 0;

    bfam_domain_get_subdomains((bfam_domain_t*)domain, BFAM_DOMAIN_OR,
        glue, domain->base.numSubdomains, subdomains, &numSubdomains);

    BFAM_LDEBUG("Number of local and parallel glue grids %jd",
        (intmax_t) numSubdomains);

    for(bfam_locidx_t s = 0; s < numSubdomains; ++s)
    {
      failures +=
        check_pm((bfam_subdomain_dgx_t*)subdomains[s], "p1", 1);
      failures +=
        check_pm((bfam_subdomain_dgx_t*)subdomains[s], "p2", 1);
      failures +=
        check_pm((bfam_subdomain_dgx_t*)subdomains[s], "p3", 1);
      failures +=
        check_pm((bfam_subdomain_dgx_t*)subdomains[s], "p4", 1);
      failures +=
        check_pm((bfam_subdomain_dgx_t*)subdomains[s], "p5", 1);
      failures +=
        check_pm((bfam_subdomain_dgx_t*)subdomains[s], "p6", 1);

      failures +=
        check_pm((bfam_subdomain_dgx_t*)subdomains[s], "Tn",   1);
      failures +=
        check_pm((bfam_subdomain_dgx_t*)subdomains[s], "Tp1", -1);
      failures +=
        check_pm((bfam_subdomain_dgx_t*)subdomains[s], "Tp2", -1);
      failures +=
        check_pm((bfam_subdomain_dgx_t*)subdomains[s], "Tp3", -1);

      failures +=
        check_pm((bfam_subdomain_dgx_t*)subdomains[s], "Sn",   1);
      failures +=
        check_pm((bfam_subdomain_dgx_t*)subdomains[s], "Sp1", -1);
      failures +=
        check_pm((bfam_subdomain_dgx_t*)subdomains[s], "Sp2", -1);
      failures +=
        check_pm((bfam_subdomain_dgx_t*)subdomains[s], "Sp3", -1);

      failures +=
        check_pm((bfam_subdomain_dgx_t*)subdomains[s], "vn", -1);
      failures +=
        check_pm((bfam_subdomain_dgx_t*)subdomains[s], "vp1", 1);
      failures +=
        check_pm((bfam_subdomain_dgx_t*)subdomains[s], "vp2", 1);
      failures +=
        check_pm((bfam_subdomain_dgx_t*)subdomains[s], "vp3", 1);

      failures +=
        check_pm((bfam_subdomain_dgx_t*)subdomains[s], "un", -1);
      failures +=
        check_pm((bfam_subdomain_dgx_t*)subdomains[s], "up1", 1);
      failures +=
        check_pm((bfam_subdomain_dgx_t*)subdomains[s], "up2", 1);
      failures +=
        check_pm((bfam_subdomain_dgx_t*)subdomains[s], "up3", 1);

      failures +=
        check_pm((bfam_subdomain_dgx_t*)subdomains[s], "_grid_nx0", -1);
      failures +=
        check_pm((bfam_subdomain_dgx_t*)subdomains[s], "_grid_nx1", -1);
      failures +=
        check_pm((bfam_subdomain_dgx_t*)subdomains[s], "_grid_nx2", -1);
    }
    bfam_free(subdomains);
  }


  bfam_communicator_free(communicator);
  bfam_free(communicator);
  return failures;
}

static int
run_tests(MPI_Comm mpicomm)
{
  int failures = 0;

  state_t astate; state_t *state = &astate;

  failures += build_state(mpicomm, state);
  failures += test_conn(mpicomm, state);

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

