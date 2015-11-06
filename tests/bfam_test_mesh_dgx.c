#include <bfam.h>
#include <bfam_domain_pxest_2.h>

#define REAL_APPROX_EQ(x, y, K)                                                \
  BFAM_APPROX_EQ((x), (y), (K), BFAM_REAL_ABS, BFAM_REAL_EPS, (K)*BFAM_REAL_EPS)

static int refine_level = 0;

#define DIM 2

#define bfam_domain_pxest_t BFAM_APPEND_EXPAND(bfam_domain_pxest_t_, DIM)
#define bfam_domain_pxest_init_callback                                        \
  BFAM_APPEND_EXPAND(bfam_domain_pxest_init_callback_, DIM)
#define bfam_domain_pxest_new BFAM_APPEND_EXPAND(bfam_domain_pxest_new_, DIM)
#define bfam_domain_pxest_init BFAM_APPEND_EXPAND(bfam_domain_pxest_init_, DIM)
#define bfam_domain_pxest_free BFAM_APPEND_EXPAND(bfam_domain_pxest_free_, DIM)
#define bfam_domain_pxest_split_dgx_subdomains                                 \
  BFAM_APPEND_EXPAND(bfam_domain_pxest_split_dgx_subdomains_, DIM)
#define bfam_domain_pxest_create_mesh                                          \
  BFAM_APPEND_EXPAND(bfam_domain_pxest_create_mesh_, DIM)

static int refine_fn(p4est_t *pxest, p4est_topidx_t which_tree,
                     p4est_quadrant_t *quadrant)
{
  if ((int)quadrant->level >= (refine_level - (int)(which_tree % 3)))
  {
    return 0;
  }
  if (quadrant->level == 1 && p4est_quadrant_child_id(quadrant) == 3)
  {
    return 1;
  }
  if (quadrant->x == P4EST_LAST_OFFSET(2) &&
      quadrant->y == P4EST_LAST_OFFSET(2))
  {
    return 1;
  }
  if (quadrant->x >= P4EST_QUADRANT_LEN(2))
  {
    return 0;
  }

  return 1;
}

static void poly0_field(bfam_locidx_t npoints, const char *name,
                        bfam_real_t time, bfam_real_t *restrict x,
                        bfam_real_t *restrict y, bfam_real_t *restrict z,
                        struct bfam_subdomain *s, void *arg,
                        bfam_real_t *restrict field)
{
  BFAM_ASSUME_ALIGNED(field, 32);

  for (bfam_locidx_t n = 0; n < npoints; ++n)
    field[n] = 4;
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

  for (bfam_locidx_t n = 0; n < npoints; ++n)
    field[n] = -x[n] - y[n] * x[n] * y[n];
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

  for (bfam_locidx_t n = 0; n < npoints; ++n)
    field[n] = x[n] + y[n];
}

static void poly3_field(bfam_locidx_t npoints, const char *name,
                        bfam_real_t time, bfam_real_t *restrict x,
                        bfam_real_t *restrict y, bfam_real_t *restrict z,
                        struct bfam_subdomain *s, void *arg,
                        bfam_real_t *restrict field)
{
  BFAM_ASSUME_ALIGNED(x, 32);
  BFAM_ASSUME_ALIGNED(y, 32);
  BFAM_ASSUME_ALIGNED(z, 32);
  BFAM_ASSUME_ALIGNED(field, 32);

  for (bfam_locidx_t n = 0; n < npoints; ++n)
    field[n] = x[n] + y[n] * x[n] + x[n] * x[n] * y[n];
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

  for (bfam_locidx_t n = 0; n < npoints; ++n)
    field[n] = y[n] + y[n] * x[n] + x[n] * x[n] * y[n];
}

static void poly5_field(bfam_locidx_t npoints, const char *name,
                        bfam_real_t time, bfam_real_t *restrict x,
                        bfam_real_t *restrict y, bfam_real_t *restrict z,
                        struct bfam_subdomain *s, void *arg,
                        bfam_real_t *restrict field)
{
  BFAM_ASSUME_ALIGNED(x, 32);
  BFAM_ASSUME_ALIGNED(y, 32);
  BFAM_ASSUME_ALIGNED(z, 32);
  BFAM_ASSUME_ALIGNED(field, 32);

  for (bfam_locidx_t n = 0; n < npoints; ++n)
    field[n] = x[n] * y[n] + y[n] * x[n] + x[n] * x[n] * y[n];
}

static void poly6_field(bfam_locidx_t npoints, const char *name,
                        bfam_real_t time, bfam_real_t *restrict x,
                        bfam_real_t *restrict y, bfam_real_t *restrict z,
                        struct bfam_subdomain *s, void *arg,
                        bfam_real_t *restrict field)
{
  BFAM_ASSUME_ALIGNED(x, 32);
  BFAM_ASSUME_ALIGNED(y, 32);
  BFAM_ASSUME_ALIGNED(z, 32);
  BFAM_ASSUME_ALIGNED(field, 32);

  for (bfam_locidx_t n = 0; n < npoints; ++n)
    field[n] = x[n] * x[n] + y[n] * y[n] + x[n] + y[n];
}

static int check_pm(bfam_subdomain_dgx_t *sub, const char *name,
                    bfam_real_t fac)
{
  int failures = 0;
  bfam_real_t *f_m =
      bfam_dictionary_get_value_ptr(&sub->base.glue_m->fields, name);
  bfam_real_t *f_p =
      bfam_dictionary_get_value_ptr(&sub->base.glue_p->fields, name);

  BFAM_ASSERT(f_m != NULL);
  BFAM_ASSERT(f_p != NULL);

  BFAM_LDEBUG(
      "Testing subdomain (%2jd, %2jd) -- (%2jd, %2jd) field %s",
      (intmax_t)sub->base.glue_m->rank, (intmax_t)sub->base.glue_m->id_s,
      (intmax_t)sub->base.glue_p->rank, (intmax_t)sub->base.glue_p->id_s, name);

  bfam_subdomain_dgx_glue_data_t *glue_p =
      (bfam_subdomain_dgx_glue_data_t *)sub->base.glue_p;
  bfam_subdomain_dgx_glue_data_t *glue_m =
      (bfam_subdomain_dgx_glue_data_t *)sub->base.glue_m;

  for (bfam_locidx_t i = 0; i < sub->K; ++i)
  {
    BFAM_LDEBUG("Testing element %2jd face %d h %d o %d",
                (intmax_t)glue_m->EToEm[i], glue_m->EToFm[i], glue_m->EToHm[i],
                glue_p->EToOp[i]);

    for (int j = 0; j < sub->Np; ++j)
    {
      size_t idx = i * sub->Np + j;
      int fail = !REAL_APPROX_EQ(f_m[idx], fac * f_p[idx], 1000);

      if (fail)
        BFAM_LDEBUG("Fail Match fm[%2d][%2d] = %20" BFAM_REAL_PRIe
                    "    fp[%2d][%2d] = %20" BFAM_REAL_PRIe,
                    i, j, f_m[i * sub->Np + j], i, j,
                    fac * f_p[i * sub->Np + j]);

      failures += fail;
    }
  }

  if (failures > 0)
    BFAM_WARNING("FAIL! %s", name);
  return failures;
}

static int check_vmaps(bfam_subdomain_dgx_t *sub, const char *name)
{
  int failures = 0;
  bfam_real_t *f = bfam_dictionary_get_value_ptr(&sub->base.fields, name);

  BFAM_ASSERT(f != NULL);

  for (bfam_locidx_t i = 0; i < sub->K * sub->Ngp[0] * sub->Ng[0]; ++i)
  {
    int fail = !REAL_APPROX_EQ(f[sub->vmapM[i]], f[sub->vmapP[i]], 1000);

    if (fail)
      BFAM_LDEBUG("Fail Match fm[%2jd] = %20" BFAM_REAL_PRIe
                  "    fp[%2jd] = %20" BFAM_REAL_PRIe,
                  (intmax_t)i, f[sub->vmapM[i]], (intmax_t)i, f[sub->vmapP[i]]);

    failures += fail;
  }

  if (failures > 0)
    BFAM_WARNING("FAIL! %s", name);
  return failures;
}

#if 0
static int
check_back(bfam_subdomain_dgx_quad_glue_t *sub, const char *name)
{
  int failures = 0;
  bfam_real_t *f_m = bfam_dictionary_get_value_ptr(&sub->base.fields_m, name);
  bfam_real_t *f_p = bfam_dictionary_get_value_ptr(&sub->base.fields_p, name);
  bfam_real_t *sub_m_f = bfam_dictionary_get_value_ptr(&sub->sub_m->base.fields,
                                                       name);

  BFAM_ASSERT(f_m != NULL);
  BFAM_ASSERT(f_p != NULL);
  BFAM_ASSERT(sub_m_f != NULL);

  BFAM_LDEBUG("Testing subdomain (%2jd, %2jd) -- (%2jd, %2jd)",
      (intmax_t)sub->rank_m, (intmax_t)sub->s_m,
      (intmax_t)sub->rank_p, (intmax_t)sub->s_p);

  /* check the projection back */
  /*
   * only check when orders match b/c exact mass is only exact for plus side
   * when the orders are the same 
   */
  bfam_subdomain_dgx_quad_t * sub_m = sub->sub_m;
  if(sub->N == sub->N_m)
  {
    bfam_locidx_t Nfp = sub_m->Nfp;
    bfam_real_t field_m[Nfp];
    bfam_real_t Mf[Nfp];
    bfam_real_t MPf[Nfp];
    bfam_real_t* mass = sub->exact_mass;
    bfam_real_t** MP   = sub->massprojection;
    for(bfam_locidx_t le = 0; le < sub->K; ++le)
    {
      bfam_locidx_t e = sub->EToEm[le];
      int8_t face = sub->EToFm[le];

      /* make sure that we have both the hanging guys on this glue grid */
      if(sub->EToHm[le] == 1 && sub->EToHm[le+1] == 2 &&  e == sub->EToEm[le+1])
      {
        for(bfam_locidx_t n = 0;n < Nfp;n++)
        {
          bfam_locidx_t f = n + Nfp*(face + 4*e);
          bfam_locidx_t iM = sub_m->vmapM[f];
          field_m[n] = sub_m_f[iM];
          Mf[n] = 0;
        }

        for(bfam_locidx_t j = 0; j < Nfp; ++j)
          for(bfam_locidx_t i = 0; i < Nfp; ++i)
            Mf[i] += mass[j * Nfp + i] * field_m[j];

        for(bfam_locidx_t i = 0; i < Nfp; i++) MPf[i] = 0;
        for(bfam_locidx_t i = 0; i < Nfp; i++)
          for(bfam_locidx_t j = 0; j < Nfp; j++)
          {
            MPf[i] += MP[1][i+j*Nfp]*f_m[le*sub->Np + j];
            MPf[i] += MP[2][i+j*Nfp]*f_m[(le+1)*sub->Np + j];
          }

        for(bfam_locidx_t n = 0; n < Nfp; n++)
        {
          int fail = !REAL_APPROX_EQ(Mf[n], MPf[n], 10);
          if(fail)
            BFAM_LDEBUG("mass projection fail on element %d face %d node %d:"
                "M*f = %e and MP[1]*I[1]*f + MP[2]*I[2]*f = %e",
                (int)e, (int)face,(int)n,(double)Mf[n],(double)MPf[n]);
          failures += fail;
        }
      }
      else if(sub->EToHm[le] == 0)
      {
        for(bfam_locidx_t n = 0;n < Nfp;n++)
        {
          bfam_locidx_t f = n + Nfp*(face + 4*e);
          bfam_locidx_t iM = sub_m->vmapM[f];
          field_m[n] = sub_m_f[iM];
          Mf[n] = 0;
        }

        for(bfam_locidx_t j = 0; j < Nfp; ++j)
          for(bfam_locidx_t i = 0; i < Nfp; ++i)
            Mf[i] += mass[j * Nfp + i] * field_m[j];

        for(bfam_locidx_t i = 0; i < Nfp; i++) MPf[i] = 0;
        for(bfam_locidx_t i = 0; i < Nfp; i++)
          for(bfam_locidx_t j = 0; j < Nfp; j++)
          {
            MPf[i] += MP[0][i+j*Nfp]*f_m[le*sub->Np + j];
          }

        for(bfam_locidx_t n = 0; n < Nfp; n++)
        {
          int fail = !REAL_APPROX_EQ(Mf[n], MPf[n], 10);
          if(fail)
            BFAM_LDEBUG("mass projection fail on element %d face %d node %d:"
                "M*f = %e and MP[0]*I[0]*f = %e",
                (int)e, (int)face,(int)n,(double)Mf[n],(double)MPf[n]);
          failures += fail;
        }
      }
    }
  }

  return failures;
}
#endif

static int build_mesh(MPI_Comm mpicomm)
{
  int failures = 0;
  int rank;
  BFAM_MPI_CHECK(MPI_Comm_rank(mpicomm, &rank));

  p4est_connectivity_t *conn = p4est_connectivity_new_disk();

  bfam_domain_pxest_t *domain = bfam_domain_pxest_new(mpicomm, conn);

  refine_level = 4;
  p4est_refine(domain->pxest, 2, refine_fn, bfam_domain_pxest_init_callback);
  p4est_balance(domain->pxest, P4EST_CONNECT_CORNER,
                bfam_domain_pxest_init_callback);
  p4est_partition(domain->pxest, 1, NULL);

  /*
  p4est_vtk_write_file(domain->pxest, NULL, "p4est_mesh");
  */

  bfam_locidx_t num_subdomains = 5;
  bfam_locidx_t *subdomainID =
      bfam_malloc(domain->pxest->local_num_quadrants * sizeof(bfam_locidx_t));
  bfam_locidx_t *N = bfam_malloc(num_subdomains * sizeof(int));

  /*
   * Create an arbitrary splitting of the domain to test things.
   *
   * When use a subdomain id independent of MPI partition.  In practice
   * the subdomain id will be selected based on physics, element type, element
   * order, etc.
   *
   * For no particular reason increase element order with id
   */
  BFAM_ROOT_INFO("Splitting pxest into %jd DG Quad subdomains",
                 (intmax_t)num_subdomains);
  for (bfam_locidx_t id = 0; id < num_subdomains; ++id)
  {
    N[id] = 3 + id;

    p4est_gloidx_t first = p4est_partition_cut_gloidx(
        domain->pxest->global_num_quadrants, id, num_subdomains);

    p4est_gloidx_t last =
        p4est_partition_cut_gloidx(domain->pxest->global_num_quadrants, id + 1,
                                   num_subdomains) -
        1;

    BFAM_ROOT_INFO("  id:%jd N:%d GIDs:%jd--%jd", (intmax_t)id, N[id],
                   (intmax_t)first, (intmax_t)last);
  }

  p4est_gloidx_t gkOffset = domain->pxest->global_first_quadrant[rank];

  bfam_locidx_t idStart = 0;
  while (gkOffset >
         p4est_partition_cut_gloidx(domain->pxest->global_num_quadrants,
                                    idStart + 1, num_subdomains) -
             1)
    ++idStart;

  for (p4est_locidx_t lk = 0, id = idStart;
       lk < domain->pxest->local_num_quadrants; ++lk)
  {
    p4est_gloidx_t gk = gkOffset + lk;

    if (gk > p4est_partition_cut_gloidx(domain->pxest->global_num_quadrants,
                                        id + 1, num_subdomains) -
                 1)
      ++id;

    BFAM_ASSERT(
        (gk >= p4est_partition_cut_gloidx(domain->pxest->global_num_quadrants,
                                          id, num_subdomains)) &&
        (gk < p4est_partition_cut_gloidx(domain->pxest->global_num_quadrants,
                                         id + 1, num_subdomains)));

    subdomainID[lk] = id;
  }

  bfam_domain_pxest_split_dgx_subdomains(domain, num_subdomains, subdomainID,
                                         NULL, N, NULL, NULL, NULL);
  bfam_domain_pxest_create_mesh(domain, NULL, NULL);

  const char *volume[] = {"_volume", NULL};
  const char *glue[] = {"_glue_parallel", "_glue_local", NULL};

  bfam_domain_add_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, volume, "p0");
  bfam_domain_add_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, volume, "p1");
  bfam_domain_add_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, volume, "p2");
  bfam_domain_add_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, volume, "p3");
  bfam_domain_add_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, volume, "p4");
  bfam_domain_add_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, volume, "p5");
  bfam_domain_add_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, volume, "p6");

  bfam_domain_add_minus_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, glue,
                              "p1");
  bfam_domain_add_minus_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, glue,
                              "p2");
  bfam_domain_add_minus_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, glue,
                              "p3");
  bfam_domain_add_minus_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, glue,
                              "p4");
  bfam_domain_add_minus_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, glue,
                              "p5");
  bfam_domain_add_minus_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, glue,
                              "p6");

  bfam_domain_add_plus_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, glue,
                             "p1");
  bfam_domain_add_plus_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, glue,
                             "p2");
  bfam_domain_add_plus_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, glue,
                             "p3");
  bfam_domain_add_plus_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, glue,
                             "p4");
  bfam_domain_add_plus_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, glue,
                             "p5");
  bfam_domain_add_plus_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, glue,
                             "p6");

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
  bfam_domain_init_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, volume, "p5",
                         0, poly5_field, NULL);
  bfam_domain_init_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, volume, "p6",
                         0, poly6_field, NULL);

  bfam_subdomain_comm_args_t commargs;

  const char *comm_args_face_scalars[] = {"_grid_nx0", "_grid_nx1", NULL};
  const char *comm_args_scalars[] = {"p1", "p2", "p3", "p4", "p5", "p6", NULL};
  const char *comm_args_vectors[] = {"v", "u", NULL};
  const char *comm_args_vector_components[] = {"p1", "p2", "p3", "p4",
                                               "p5", "p6", NULL};
  const char *comm_args_tensors[] = {"T", "S", NULL};
  const char *comm_args_tensor_components[] = {"p1", "p2", "p3", "p4", "p5",
                                               "p6", "p1", "p3", "p5", "p2",
                                               "p4", "p6", NULL};
  commargs.scalars_m = comm_args_scalars;
  commargs.scalars_p = comm_args_scalars;

  commargs.vectors_m = comm_args_vectors;
  commargs.vectors_p = comm_args_vectors;
  commargs.vector_components_m = comm_args_vector_components;
  commargs.vector_components_p = comm_args_vector_components;

  commargs.tensors_m = comm_args_tensors;
  commargs.tensors_p = comm_args_tensors;
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
  for (int f = 0; comm_args_face_scalars[f] != NULL; f++)
  {
    bfam_domain_add_minus_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, glue,
                                comm_args_face_scalars[f]);
    bfam_domain_add_plus_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, glue,
                               comm_args_face_scalars[f]);
  }
  for (int f = 0; comm_args_scalars[f] != NULL; f++)
  {
    bfam_domain_add_minus_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, glue,
                                comm_args_scalars[f]);
    bfam_domain_add_plus_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, glue,
                               comm_args_scalars[f]);
  }
  for (int f = 0; comm_args_vectors[f] != NULL; f++)
  {
    char name[BFAM_BUFSIZ];
    snprintf(name, BFAM_BUFSIZ, "%sn", comm_args_vectors[f]);
    bfam_domain_add_minus_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, glue,
                                name);
    bfam_domain_add_plus_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, glue,
                               name);
    snprintf(name, BFAM_BUFSIZ, "%sp1", comm_args_vectors[f]);
    bfam_domain_add_minus_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, glue,
                                name);
    bfam_domain_add_plus_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, glue,
                               name);
    snprintf(name, BFAM_BUFSIZ, "%sp2", comm_args_vectors[f]);
    bfam_domain_add_minus_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, glue,
                                name);
    bfam_domain_add_plus_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, glue,
                               name);
    snprintf(name, BFAM_BUFSIZ, "%sp3", comm_args_vectors[f]);
    bfam_domain_add_minus_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, glue,
                                name);
    bfam_domain_add_plus_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, glue,
                               name);
  }
  for (int f = 0; comm_args_tensors[f] != NULL; f++)
  {
    char name[BFAM_BUFSIZ];
    snprintf(name, BFAM_BUFSIZ, "%sn", comm_args_tensors[f]);
    bfam_domain_add_minus_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, glue,
                                name);
    bfam_domain_add_plus_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, glue,
                               name);
    snprintf(name, BFAM_BUFSIZ, "%sp1", comm_args_tensors[f]);
    bfam_domain_add_minus_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, glue,
                                name);
    bfam_domain_add_plus_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, glue,
                               name);
    snprintf(name, BFAM_BUFSIZ, "%sp2", comm_args_tensors[f]);
    bfam_domain_add_minus_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, glue,
                                name);
    bfam_domain_add_plus_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, glue,
                               name);
    snprintf(name, BFAM_BUFSIZ, "%sp3", comm_args_tensors[f]);
    bfam_domain_add_minus_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, glue,
                                name);
    bfam_domain_add_plus_field((bfam_domain_t *)domain, BFAM_DOMAIN_OR, glue,
                               name);
  }

  bfam_communicator_t *communicator = bfam_communicator_new(
      (bfam_domain_t *)domain, BFAM_DOMAIN_OR, glue, mpicomm, 10, &commargs);

  /* start recv_send */
  bfam_communicator_start(communicator);

  /* finish recv */
  bfam_communicator_finish(communicator);

  /*
  const char *ps[] = {"p1", "p2", "p3", "p4", "p5", "p6", NULL};

  bfam_vtk_write_file((bfam_domain_t*)domain, BFAM_DOMAIN_OR, volume,
                       NULL,"ps",0, ps, NULL, NULL, 0, 0, 0);
  */

  /*
   * Check local subdomain vmaps
   */
  {
    bfam_subdomain_t **subdomains =
        bfam_malloc(domain->base.num_subdomains * sizeof(bfam_subdomain_t **));

    bfam_locidx_t num_subdomains = 0;

    bfam_domain_get_subdomains((bfam_domain_t *)domain, BFAM_DOMAIN_OR, volume,
                               domain->base.num_subdomains, subdomains,
                               &num_subdomains);

    BFAM_LDEBUG("Number of volume subdomains %jd", (intmax_t)num_subdomains);

    for (bfam_locidx_t s = 0; s < num_subdomains; ++s)
    {
      failures += check_vmaps((bfam_subdomain_dgx_t *)subdomains[s], "p1");
      failures += check_vmaps((bfam_subdomain_dgx_t *)subdomains[s], "p2");
      failures += check_vmaps((bfam_subdomain_dgx_t *)subdomains[s], "p3");
      failures += check_vmaps((bfam_subdomain_dgx_t *)subdomains[s], "p4");
      failures += check_vmaps((bfam_subdomain_dgx_t *)subdomains[s], "p5");
      failures += check_vmaps((bfam_subdomain_dgx_t *)subdomains[s], "p6");
    }
  }

  /*
   * Check to see if neighboring values got communicated
   */
  {
    bfam_subdomain_t **subdomains =
        bfam_malloc(domain->base.num_subdomains * sizeof(bfam_subdomain_t **));

    bfam_locidx_t num_subdomains = 0;

    bfam_domain_get_subdomains((bfam_domain_t *)domain, BFAM_DOMAIN_OR, glue,
                               domain->base.num_subdomains, subdomains,
                               &num_subdomains);

    BFAM_LDEBUG("Number of local and parallel glue grids %jd",
                (intmax_t)num_subdomains);

    for (bfam_locidx_t s = 0; s < num_subdomains; ++s)
    {
      // failures +=
      //   check_back((bfam_subdomain_dgx_t*)subdomains[s], "p1");
      // failures +=
      //   check_back((bfam_subdomain_dgx_t*)subdomains[s], "p2");
      // failures +=
      //   check_back((bfam_subdomain_dgx_t*)subdomains[s], "p3");
      // failures +=
      //   check_back((bfam_subdomain_dgx_t*)subdomains[s], "p4");
      // failures +=
      //   check_back((bfam_subdomain_dgx_t*)subdomains[s], "p5");
      // failures +=
      //   check_back((bfam_subdomain_dgx_t*)subdomains[s], "p6");

      /* last argument lets us change the sign if necessary */
      failures += check_pm((bfam_subdomain_dgx_t *)subdomains[s], "p1", 1);
      failures += check_pm((bfam_subdomain_dgx_t *)subdomains[s], "p2", 1);
      failures += check_pm((bfam_subdomain_dgx_t *)subdomains[s], "p3", 1);
      failures += check_pm((bfam_subdomain_dgx_t *)subdomains[s], "p4", 1);
      failures += check_pm((bfam_subdomain_dgx_t *)subdomains[s], "p5", 1);
      failures += check_pm((bfam_subdomain_dgx_t *)subdomains[s], "p6", 1);

      failures += check_pm((bfam_subdomain_dgx_t *)subdomains[s], "Tn", 1);
      failures += check_pm((bfam_subdomain_dgx_t *)subdomains[s], "Tp1", -1);
      failures += check_pm((bfam_subdomain_dgx_t *)subdomains[s], "Tp2", -1);
      failures += check_pm((bfam_subdomain_dgx_t *)subdomains[s], "Tp3", -1);

      failures += check_pm((bfam_subdomain_dgx_t *)subdomains[s], "Sn", 1);
      failures += check_pm((bfam_subdomain_dgx_t *)subdomains[s], "Sp1", -1);
      failures += check_pm((bfam_subdomain_dgx_t *)subdomains[s], "Sp2", -1);
      failures += check_pm((bfam_subdomain_dgx_t *)subdomains[s], "Sp3", -1);

      failures += check_pm((bfam_subdomain_dgx_t *)subdomains[s], "vn", -1);
      failures += check_pm((bfam_subdomain_dgx_t *)subdomains[s], "vp1", 1);
      failures += check_pm((bfam_subdomain_dgx_t *)subdomains[s], "vp2", 1);
      failures += check_pm((bfam_subdomain_dgx_t *)subdomains[s], "vp3", 1);

      failures += check_pm((bfam_subdomain_dgx_t *)subdomains[s], "un", -1);
      failures += check_pm((bfam_subdomain_dgx_t *)subdomains[s], "up1", 1);
      failures += check_pm((bfam_subdomain_dgx_t *)subdomains[s], "up2", 1);
      failures += check_pm((bfam_subdomain_dgx_t *)subdomains[s], "up3", 1);

      failures +=
          check_pm((bfam_subdomain_dgx_t *)subdomains[s], "_grid_nx0", -1);
      failures +=
          check_pm((bfam_subdomain_dgx_t *)subdomains[s], "_grid_nx1", -1);
    }

    bfam_free(subdomains);
  }

  /* clean up */
  bfam_communicator_free(communicator);
  bfam_free(communicator);

  bfam_free(subdomainID);
  bfam_free(N);

  bfam_domain_pxest_free(domain);
  bfam_free(domain);
  p4est_connectivity_destroy(conn);

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
      "  there are four possible options to this program, some of which have\n"
      "  multiple names:\n"
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

  failures += build_mesh(comm);

  sc_finalize();
  BFAM_MPI_CHECK(MPI_Finalize());

  bfam_gopt_free(options);

  if (failures > 0)
    BFAM_WARNING("FAIL! with failures %d", failures);
  return failures;
}
