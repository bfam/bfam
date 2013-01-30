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
  MPI_Comm comm = MPI_COMM_WORLD;
  p4est_t *p4est;
  p4est_connectivity_t *conn;
  p4est_ghost_t *ghost_layer;
  p4est_lnodes_t *lnodes;
  int rank;
  const int degree = 1;

  BFAM_MPI_CHECK(MPI_Init(&argc,&argv));
  BFAM_MPI_CHECK(MPI_Comm_rank(comm, &rank));

  bfam_log_init(rank, stdout, BFAM_LL_DEFAULT);
  bfam_signal_handler_set();

  sc_init(comm, 0, 0, NULL, SC_LP_DEFAULT);
  p4est_init(NULL, SC_LP_DEFAULT);

  conn = p4est_connectivity_new_corner();
  p4est = p4est_new_ext(comm, conn, 0, 0, 0, 0, NULL, NULL);

  refine_level = 1;
  p4est_refine(p4est, 1, refine_fn, NULL);
  p4est_balance(p4est, P4EST_CONNECT_FACE, NULL);
  p4est_partition(p4est, NULL);

  ghost_layer = p4est_ghost_new(p4est, P4EST_CONNECT_FULL);
  lnodes = p4est_lnodes_new(p4est, ghost_layer, degree);


  /*
   * Output the mesh.  It can be read using something like following command:
   *
   * mpirun -np 3 ./bfam_exam_p4est | grep MESH | sort -n -k 2 | sort -n -k 5 | gvim -
   */
  fflush(stdout);
  BFAM_MPI_CHECK(MPI_Barrier(comm));
  BFAM_ROOT_INFO("MESH 0 ------------ Mesh Begin ------------");
  BFAM_ROOT_INFO("MESH 1 degree  = %d", lnodes->degree);
  BFAM_ROOT_INFO("MESH 2 vnodes = %d", lnodes->vnodes);
  BFAM_INFO("MESH 3 num_local_elements  = %jd", (intmax_t)lnodes->num_local_elements);
  BFAM_INFO("MESH 4 num_local_nodes = %jd", (intmax_t)lnodes->num_local_nodes);
  BFAM_INFO("MESH 5 owned_count = %jd", (intmax_t)lnodes->owned_count);
  BFAM_INFO("MESH 6 global_offset = %jd", (intmax_t)lnodes->global_offset);


  sc_array_t *global_nodes = sc_array_new(sizeof (p4est_gloidx_t));
  sc_array_resize(global_nodes, lnodes->num_local_nodes);
  for(size_t zz = 0; zz < global_nodes->elem_count; ++zz)
  {
    *((p4est_gloidx_t *) sc_array_index(global_nodes, zz)) =
      p4est_lnodes_global_index(lnodes, zz);
  }

  p4est_lnodes_share_owned(global_nodes, lnodes);

  for(size_t zz = 0; zz < global_nodes->elem_count; ++zz)
  {
    const p4est_gloidx_t gn =
      *((p4est_gloidx_t *)sc_array_index(global_nodes, zz));
    SC_CHECK_ABORT (gn == p4est_lnodes_global_index(lnodes, zz),
        "Lnodes: bad global index across procesors");
    BFAM_INFO("MESH 7 global_nodes[%zu] = %jd", zz, (intmax_t)gn);
  }

  sc_array_destroy(global_nodes);

  BFAM_ROOT_INFO("MESH 8 ------------ Mesh End ------------");



  p4est_vtk_write_file(p4est, NULL, "mesh");

  p4est_lnodes_destroy(lnodes);
  p4est_ghost_destroy(ghost_layer);
  p4est_destroy(p4est);
  p4est_connectivity_destroy(conn);

  sc_finalize();
  BFAM_MPI_CHECK(MPI_Finalize());

  return EXIT_SUCCESS;
}

