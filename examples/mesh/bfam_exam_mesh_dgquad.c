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
  int rank;

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

  p4est_vtk_write_file(p4est, NULL, "p4est_mesh");

  p4est_destroy(p4est);
  p4est_connectivity_destroy(conn);

  sc_finalize();
  BFAM_MPI_CHECK(MPI_Finalize());

  return EXIT_SUCCESS;
}

