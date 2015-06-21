#include <bfam.h>

typedef struct bfam_subdomain_comm_test
{
  bfam_subdomain_t base;
  bfam_locidx_t np;   /*neighbor ID*/
  bfam_locidx_t ns;   /*neighbor sub ID*/
  bfam_locidx_t ms;   /*my sub ID*/
  bfam_locidx_t rank; /*my rank*/
  bfam_locidx_t send_num;
  bfam_locidx_t recv_num;
  bfam_real_t send_data;
  bfam_real_t recv_data;
} bfam_subdomain_comm_test_t;

void bfam_subdomain_glue_put(bfam_subdomain_t *thisSubdomain, void *send_buf,
                             size_t send_sz, void *args)
{
  int num_reals = send_sz / sizeof(bfam_real_t);
  bfam_real_t *buffer = (bfam_real_t *)send_buf;
  bfam_subdomain_comm_test_t *sub = (bfam_subdomain_comm_test_t *)thisSubdomain;
  BFAM_ABORT_IF_NOT(num_reals == sub->send_num, "Wrong send size");
  for (int i = 0; i < num_reals; i++)
    buffer[i] = sub->send_data;
}

void bfam_subdomain_glue_get(bfam_subdomain_t *thisSubdomain, void *recv_buf,
                             size_t recv_sz, void *args)
{
  int num_reals = recv_sz / sizeof(bfam_real_t);
  bfam_real_t *buffer = (bfam_real_t *)recv_buf;
  bfam_subdomain_comm_test_t *sub = (bfam_subdomain_comm_test_t *)thisSubdomain;
  BFAM_ABORT_IF_NOT(num_reals == sub->recv_num, "Wrong recv size");
  for (int i = 0; i < num_reals; i++)
    BFAM_ABORT_IF_NOT(buffer[i] == sub->recv_data,
                      "Recv wrong data: got: %f wanted: %f from %d", buffer[i],
                      sub->recv_data, sub->np);
}

void bfam_subdomain_comm_test_info(bfam_subdomain_t *thisSubdomain, int *rank,
                                   bfam_locidx_t *sort, int num_sort,
                                   size_t *send_sz, size_t *recv_sz, void *args)
{
  BFAM_ASSERT(num_sort > 1);
  bfam_subdomain_comm_test_t *sub = (bfam_subdomain_comm_test_t *)thisSubdomain;
  *rank = sub->np;
  sort[0] = sub->ns; /* neighbor ID */
  sort[1] = sub->ms; /* my ID */

  *send_sz = sub->send_num * sizeof(bfam_real_t);
  *recv_sz = sub->recv_num * sizeof(bfam_real_t);
}

int main(int argc, char *argv[])
{
  int rank, size;

  BFAM_MPI_CHECK(MPI_Init(&argc, &argv));
  BFAM_MPI_CHECK(MPI_Comm_size(MPI_COMM_WORLD, &size));
  BFAM_MPI_CHECK(MPI_Comm_rank(MPI_COMM_WORLD, &rank));
  bfam_log_init(rank, stdout, BFAM_LL_VERBOSE);

  MPI_Comm comm = MPI_COMM_WORLD;

  bfam_domain_t domain;
  bfam_domain_init(&domain, comm);

  /* Set up some fake glue grids */
  if (rank == 0)
  {
    for (int i = 0; i < size; i++)
    {
      char tmp[BFAM_BUFSIZ];
      sprintf(tmp, "%d_%d_%d_%d", rank, i, i, 10 * i);
      bfam_subdomain_comm_test_t *newSub =
          bfam_malloc(sizeof(bfam_subdomain_comm_test_t));
      bfam_subdomain_init((bfam_subdomain_t *)newSub, i, -1, tmp);
      newSub->ns = 10 * i;
      newSub->ms = i;
      newSub->np = i;
      newSub->rank = rank;
      newSub->send_num = i + 1;
      newSub->recv_num = 10 * (i + 1);
      if (i == 0)
        newSub->recv_num = (i + 1);
      newSub->send_data = i + 1;
      newSub->recv_data = 10 * (i + 1);
      if (i == 0)
        newSub->recv_data = (i + 1);
      newSub->base.glue_comm_info = &bfam_subdomain_comm_test_info;
      newSub->base.glue_put_send_buffer = &bfam_subdomain_glue_put;
      newSub->base.glue_get_recv_buffer = &bfam_subdomain_glue_get;
      bfam_subdomain_add_tag((bfam_subdomain_t *)newSub, "_glue");
      bfam_domain_add_subdomain(&domain, (bfam_subdomain_t *)newSub);
    }

    {
      char tmp[BFAM_BUFSIZ];
      sprintf(tmp, "%d_%d_%d_%d", rank, rank, size, size + 1);
      bfam_subdomain_comm_test_t *newSub =
          bfam_malloc(sizeof(bfam_subdomain_comm_test_t));
      bfam_subdomain_init((bfam_subdomain_t *)newSub, size, -1, tmp);
      newSub->ns = size + 1;
      newSub->ms = size;
      newSub->np = 0;
      newSub->rank = rank;
      newSub->send_num = size;
      newSub->recv_num = size + 1;
      newSub->send_data = size;
      newSub->recv_data = size + 1;
      newSub->base.glue_comm_info = &bfam_subdomain_comm_test_info;
      newSub->base.glue_put_send_buffer = &bfam_subdomain_glue_put;
      newSub->base.glue_get_recv_buffer = &bfam_subdomain_glue_get;
      bfam_subdomain_add_tag((bfam_subdomain_t *)newSub, "_glue");
      bfam_domain_add_subdomain(&domain, (bfam_subdomain_t *)newSub);
    }

    {
      char tmp[BFAM_BUFSIZ];
      sprintf(tmp, "%d_%d_%d_%d", rank, rank, size + 1, size);
      bfam_subdomain_comm_test_t *newSub =
          bfam_malloc(sizeof(bfam_subdomain_comm_test_t));
      bfam_subdomain_init((bfam_subdomain_t *)newSub, size + 1, -1, tmp);
      newSub->ns = size;
      newSub->ms = size + 1;
      newSub->np = 0;
      newSub->rank = rank;
      newSub->send_num = size + 1;
      newSub->recv_num = size;
      newSub->send_data = size + 1;
      newSub->recv_data = size;
      newSub->base.glue_comm_info = &bfam_subdomain_comm_test_info;
      newSub->base.glue_put_send_buffer = &bfam_subdomain_glue_put;
      newSub->base.glue_get_recv_buffer = &bfam_subdomain_glue_get;
      bfam_subdomain_add_tag((bfam_subdomain_t *)newSub, "_glue");
      bfam_domain_add_subdomain(&domain, (bfam_subdomain_t *)newSub);
    }
  }
  else
  {
    {
      char tmp[BFAM_BUFSIZ];
      sprintf(tmp, "%d_%d_%d_%d", rank, 0, 10 * rank, rank);
      bfam_subdomain_comm_test_t *newSub =
          bfam_malloc(sizeof(bfam_subdomain_comm_test_t));
      bfam_subdomain_init((bfam_subdomain_t *)newSub, 10 * rank, -1, tmp);
      newSub->ns = rank;
      newSub->ms = 10 * rank;
      newSub->np = 0;
      newSub->rank = rank;
      newSub->send_num = 10 * (rank + 1);
      newSub->recv_num = (rank + 1);
      newSub->send_data = 10 * (rank + 1);
      newSub->recv_data = (rank + 1);
      newSub->base.glue_comm_info = &bfam_subdomain_comm_test_info;
      newSub->base.glue_put_send_buffer = &bfam_subdomain_glue_put;
      newSub->base.glue_get_recv_buffer = &bfam_subdomain_glue_get;
      bfam_subdomain_add_tag((bfam_subdomain_t *)newSub, "_glue");
      bfam_domain_add_subdomain(&domain, (bfam_subdomain_t *)newSub);
    }
    {
      char tmp[BFAM_BUFSIZ];
      sprintf(tmp, "%d_%d_%d_%d", rank, rank, rank, 111 * rank);
      bfam_subdomain_comm_test_t *newSub =
          bfam_malloc(sizeof(bfam_subdomain_comm_test_t));
      bfam_subdomain_init((bfam_subdomain_t *)newSub, 111 * rank, -1, tmp);
      newSub->ns = rank;
      newSub->ms = -rank;
      newSub->np = rank;
      newSub->rank = rank;
      newSub->send_num = 111;
      newSub->recv_num = 222;
      newSub->send_data = 111;
      newSub->recv_data = 222;
      newSub->base.glue_comm_info = &bfam_subdomain_comm_test_info;
      newSub->base.glue_put_send_buffer = &bfam_subdomain_glue_put;
      newSub->base.glue_get_recv_buffer = &bfam_subdomain_glue_get;
      bfam_subdomain_add_tag((bfam_subdomain_t *)newSub, "_glue");
      bfam_domain_add_subdomain(&domain, (bfam_subdomain_t *)newSub);
    }
    {
      char tmp[BFAM_BUFSIZ];
      sprintf(tmp, "%d_%d_%d_%d", rank, rank, rank, 222 * rank);
      bfam_subdomain_comm_test_t *newSub =
          bfam_malloc(sizeof(bfam_subdomain_comm_test_t));
      bfam_subdomain_init((bfam_subdomain_t *)newSub, 222 * rank, -1, tmp);
      newSub->ns = -rank;
      newSub->ms = rank;
      newSub->np = rank;
      newSub->rank = rank;
      newSub->send_num = 222;
      newSub->recv_num = 111;
      newSub->send_data = 222;
      newSub->recv_data = 111;
      newSub->base.glue_comm_info = &bfam_subdomain_comm_test_info;
      newSub->base.glue_put_send_buffer = &bfam_subdomain_glue_put;
      newSub->base.glue_get_recv_buffer = &bfam_subdomain_glue_get;
      bfam_subdomain_add_tag((bfam_subdomain_t *)newSub, "_glue");
      bfam_domain_add_subdomain(&domain, (bfam_subdomain_t *)newSub);
    }
  }

  /* set up communicator */
  {
    const char *tags[] = {"_glue", NULL};
    bfam_communicator_t *communicator = bfam_communicator_new(
        &domain, BFAM_DOMAIN_AND, tags, MPI_COMM_WORLD, 10, 0, NULL);

    /* start recv_send */
    bfam_communicator_start(communicator);

    /* finish recv */
    bfam_communicator_finish(communicator);

    /* clean up */
    bfam_communicator_free(communicator);
    bfam_free(communicator);
  }

  {
    /* set up communicator */
    const char *tags[] = {"_glue", NULL};
    bfam_communicator_t *communicator = bfam_communicator_new(
        &domain, BFAM_DOMAIN_AND, tags, MPI_COMM_WORLD, 10, 1, NULL);

    /* start recv_send */
    bfam_communicator_start(communicator);

    /* finish recv */
    bfam_communicator_finish(communicator);

    /* clean up */
    bfam_communicator_free(communicator);
    bfam_free(communicator);
  }

  bfam_domain_free(&domain);
  BFAM_MPI_CHECK(MPI_Finalize());
}
