#include <bfam.h>

static int          refine_level = 0;

static int
refine_fn(p4est_t * p4est, p4est_topidx_t which_tree,
          p4est_quadrant_t * quadrant)
{
  if ((int) quadrant->level >= (refine_level - (int) (which_tree % 3)))
  {
    return 0;
  }
  if (quadrant->level == 1 && p4est_quadrant_child_id(quadrant) == 3)
  {
    return 1;
  }
  if (quadrant->x == P4EST_LAST_OFFSET (2) &&
      quadrant->y == P4EST_LAST_OFFSET (2))
  {
    return 1;
  }
  if (quadrant->x >= P4EST_QUADRANT_LEN (2))
  {
    return 0;
  }

  return 1;
}

int
main (int argc, char *argv[])
{
  bfam_mpicomm_t comm = { .comm = MPI_COMM_WORLD, .isMember = 1 };
  p4est_connectivity_t *conn;
  int rank;

  BFAM_MPI_CHECK(MPI_Init(&argc,&argv));
  BFAM_MPI_CHECK(MPI_Comm_rank(comm.comm, &rank));

  bfam_log_init(rank, stdout, BFAM_LL_DEFAULT);
  bfam_signal_handler_set();

  sc_init(comm.comm, 0, 0, NULL, SC_LP_DEFAULT);
  p4est_init(NULL, SC_LP_DEFAULT);

  conn = p4est_connectivity_new_corner();

  bfam_domain_p4est_t* domain = bfam_domain_p4est_new(&comm, conn);

  refine_level = 1;
  p4est_refine(domain->p4est, 1, refine_fn, NULL);
  p4est_balance(domain->p4est, P4EST_CONNECT_FACE, NULL);
  p4est_partition(domain->p4est, NULL);

  p4est_vtk_write_file(domain->p4est, NULL, "p4est_mesh");

  bfam_locidx_t numSubdomains = 4;
  bfam_locidx_t *subdomainID =
    bfam_malloc(domain->p4est->local_num_quadrants*sizeof(bfam_locidx_t));
  bfam_locidx_t *N = bfam_malloc(numSubdomains*sizeof(int));

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
      (intmax_t) numSubdomains);
  for(bfam_locidx_t id = 0; id < numSubdomains; ++id)
  {
    N[id] = 3+id;

    p4est_gloidx_t first =
      p4est_partition_cut_gloidx(domain->p4est->global_num_quadrants,
          id, numSubdomains);

    p4est_gloidx_t last =
      p4est_partition_cut_gloidx(domain->p4est->global_num_quadrants,
          id + 1, numSubdomains) - 1;

    BFAM_ROOT_INFO("  id:%jd N:%d GIDs:%jd--%jd", (intmax_t) id, N[id],
        (intmax_t) first, (intmax_t) last);
  }

  p4est_gloidx_t gkOffset = domain->p4est->global_first_quadrant[rank];

  bfam_locidx_t idStart = 0;
  while(gkOffset >
      p4est_partition_cut_gloidx(domain->p4est->global_num_quadrants,
        idStart + 1, numSubdomains) - 1) ++idStart;

  for(p4est_locidx_t lk = 0, id = idStart;
      lk < domain->p4est->local_num_quadrants;
      ++lk)
  {
    p4est_gloidx_t gk = gkOffset + lk;

    if(gk > p4est_partition_cut_gloidx(domain->p4est->global_num_quadrants,
                                       id + 1, numSubdomains) - 1)
      ++id;

    BFAM_ASSERT(
      (gk >= p4est_partition_cut_gloidx(domain->p4est->global_num_quadrants,
                                   id, numSubdomains)) &&
      (gk < p4est_partition_cut_gloidx(domain->p4est->global_num_quadrants,
                                   id + 1, numSubdomains)));

    subdomainID[lk] = id;
  }

  bfam_domain_p4est_split_dgx_quad_subdomains(domain, numSubdomains,
      subdomainID, N);

  const char *volume[] = {"_volume", NULL};
  bfam_vtk_write_file((bfam_domain_t*)domain, BFAM_DOMAIN_AND, volume,
                       "bfam_mesh", NULL, NULL, NULL);

  bfam_free(subdomainID);
  bfam_free(N);

  bfam_domain_p4est_free(domain);
  bfam_free(domain);
  p4est_connectivity_destroy(conn);

  sc_finalize();
  BFAM_MPI_CHECK(MPI_Finalize());

  return EXIT_SUCCESS;
}

