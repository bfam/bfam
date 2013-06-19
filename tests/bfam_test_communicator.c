#include <bfam.h>

typedef struct bfam_subdomain_comm_test
{
  bfam_subdomain_t base;
  bfam_locidx_t    ns; /*neighbor sub ID*/
  bfam_locidx_t    np; /*neighbor ID*/
  bfam_locidx_t  rank; /*my rank*/
} bfam_subdomain_comm_test_t;

void
bfam_subdomain_glue_put(bfam_subdomain_t *thisSubdomain,
    void *send_buf, size_t send_sz)
{
  int num_reals = send_sz/sizeof(bfam_real_t);
  bfam_real_t* buffer = (bfam_real_t*) send_buf;
  bfam_subdomain_comm_test_t* sub = (bfam_subdomain_comm_test_t*) thisSubdomain;
  BFAM_ABORT_IF_NOT(num_reals == (1+sub->base.id),
      "Wrong send size");
  for(int i = 0;i < num_reals;i++)
    buffer[i] = sub->rank*sub->base.id;
}

void
bfam_subdomain_glue_get(bfam_subdomain_t *thisSubdomain,
    void *recv_buf, size_t recv_sz)
{
  int num_reals = recv_sz/sizeof(bfam_real_t);
  bfam_real_t* buffer = (bfam_real_t*) recv_buf;
  bfam_subdomain_comm_test_t* sub = (bfam_subdomain_comm_test_t*) thisSubdomain;
  BFAM_ABORT_IF_NOT(num_reals == (1+sub->ns),
      "Wrong recv size");
  for(int i = 0;i < num_reals;i++)
    BFAM_ABORT_IF_NOT(buffer[i] == sub->np*sub->ns,"Death");
}

void
bfam_subdomain_comm_test_info(bfam_subdomain_t *thisSubdomain,
    int *rank, bfam_locidx_t *my_id, bfam_locidx_t *neigh_id,
    size_t *send_sz, size_t *recv_sz)
{
  bfam_subdomain_comm_test_t *sub = (bfam_subdomain_comm_test_t*) thisSubdomain;
  *rank = sub->np;
  *neigh_id = sub->ns;
  *my_id = sub->base.id;

  *send_sz = (1+sub->base.id)*sizeof(bfam_real_t);
  *recv_sz = (1+sub->ns)*sizeof(bfam_real_t);
}

int
main (int argc, char *argv[])
{
  int rank, size;

  BFAM_MPI_CHECK(MPI_Init(&argc,&argv));
  BFAM_MPI_CHECK(MPI_Comm_size(MPI_COMM_WORLD, &size));
  BFAM_MPI_CHECK(MPI_Comm_rank(MPI_COMM_WORLD, &rank));
  bfam_log_init(rank,stdout,BFAM_LL_VERBOSE);

  bfam_mpicomm_t comm;
  comm.comm = MPI_COMM_WORLD;
  comm.isMember = (0==0);

  bfam_domain_t domain;
  bfam_domain_init(&domain,&comm);


  /* Set up some fake glue grids */
  if(rank == 0)
  {
    for(int i = 0; i < size; i++)
    {
      char tmp[BFAM_BUFSIZ];
      sprintf(tmp,"%d_%d_%d_%d",rank,i,i,10*i);
      bfam_subdomain_comm_test_t* newSub =
        bfam_malloc(sizeof(bfam_subdomain_comm_test_t));
      bfam_subdomain_init((bfam_subdomain_t*) newSub,i,tmp);
      newSub->ns = 10*i;
      newSub->np =    i;
      newSub->rank =  rank;
      newSub->base.glue_comm_info = &bfam_subdomain_comm_test_info;
      newSub->base.glue_put_send_buffer = &bfam_subdomain_glue_put;
      newSub->base.glue_get_recv_buffer = &bfam_subdomain_glue_get;
      bfam_subdomain_add_tag((bfam_subdomain_t*) newSub,"_glue");
      bfam_domain_add_subdomain(&domain,(bfam_subdomain_t*) newSub);
    }

    {
      char tmp[BFAM_BUFSIZ];
      sprintf(tmp,"%d_%d_%d_%d",rank,rank,size,size+1);
      bfam_subdomain_comm_test_t* newSub =
        bfam_malloc(sizeof(bfam_subdomain_comm_test_t));
      bfam_subdomain_init((bfam_subdomain_t*) newSub,size,tmp);
      newSub->ns = size+1;
      newSub->np = 0;
      newSub->rank =  rank;
      newSub->base.glue_comm_info = &bfam_subdomain_comm_test_info;
      newSub->base.glue_put_send_buffer = &bfam_subdomain_glue_put;
      newSub->base.glue_get_recv_buffer = &bfam_subdomain_glue_get;
      bfam_subdomain_add_tag((bfam_subdomain_t*) newSub,"_glue");
      bfam_domain_add_subdomain(&domain,(bfam_subdomain_t*) newSub);
    }

    {
      char tmp[BFAM_BUFSIZ];
      sprintf(tmp,"%d_%d_%d_%d",rank,rank,size+1,size);
      bfam_subdomain_comm_test_t* newSub =
        bfam_malloc(sizeof(bfam_subdomain_comm_test_t));
      bfam_subdomain_init((bfam_subdomain_t*) newSub,size+1,tmp);
      newSub->ns = size;
      newSub->np = 0;
      newSub->rank =  rank;
      newSub->base.glue_comm_info = &bfam_subdomain_comm_test_info;
      newSub->base.glue_put_send_buffer = &bfam_subdomain_glue_put;
      newSub->base.glue_get_recv_buffer = &bfam_subdomain_glue_get;
      bfam_subdomain_add_tag((bfam_subdomain_t*) newSub,"_glue");
      bfam_domain_add_subdomain(&domain,(bfam_subdomain_t*) newSub);
    }
  }
  else
  {
    char tmp[BFAM_BUFSIZ];
    sprintf(tmp,"%d_%d_%d_%d",rank,0,10*rank,rank);
    bfam_subdomain_comm_test_t* newSub =
      bfam_malloc(sizeof(bfam_subdomain_comm_test_t));
    bfam_subdomain_init((bfam_subdomain_t*) newSub,10*rank,tmp);
    newSub->ns = rank;
    newSub->np =    0;
    newSub->rank =  rank;
    newSub->base.glue_comm_info = &bfam_subdomain_comm_test_info;
      newSub->base.glue_put_send_buffer = &bfam_subdomain_glue_put;
      newSub->base.glue_get_recv_buffer = &bfam_subdomain_glue_get;
    bfam_subdomain_add_tag((bfam_subdomain_t*) newSub,"_glue");
    bfam_domain_add_subdomain(&domain,(bfam_subdomain_t*) newSub);
  }

  /* set up communicator */
  const char* tags[] = {"_glue",NULL};
  bfam_communicator_t* communicator = bfam_communicator_new(&domain,
      BFAM_DOMAIN_AND,tags,MPI_COMM_WORLD,10);

  /* start recv_send */
  bfam_communicator_start(communicator);

  /* finish recv */
  bfam_communicator_finish(communicator);

  /* clean up */
  bfam_communicator_free(communicator);
  bfam_free(communicator);
  bfam_domain_free(&domain);
  BFAM_MPI_CHECK(MPI_Finalize());
}
