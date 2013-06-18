#include <bfam.h>

typedef struct bfam_subdomain_comm_test
{
  bfam_subdomain_t base;
  bfam_locidx_t    ns;
  bfam_locidx_t    np;
} bfam_subdomain_comm_test_t;

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
      sprintf(tmp,"%d_%d",rank,i);
      bfam_subdomain_comm_test_t* newSub =
        bfam_malloc(sizeof(bfam_subdomain_comm_test_t));
      bfam_subdomain_init((bfam_subdomain_t*) newSub,i,tmp);
      newSub->ns = 10*i;
      newSub->np =    i;
      bfam_subdomain_add_tag((bfam_subdomain_t*) newSub,"_glue");
      bfam_domain_add_subdomain(&domain,(bfam_subdomain_t*) newSub);
    }

    {
      char tmp[BFAM_BUFSIZ];
      sprintf(tmp,"%d_%d_%d",rank,rank,size);
      bfam_subdomain_comm_test_t* newSub =
        bfam_malloc(sizeof(bfam_subdomain_comm_test_t));
      bfam_subdomain_init((bfam_subdomain_t*) newSub,size,tmp);
      newSub->ns = size+1;
      newSub->np = 0;
      bfam_subdomain_add_tag((bfam_subdomain_t*) newSub,"_glue");
      bfam_domain_add_subdomain(&domain,(bfam_subdomain_t*) newSub);
    }

    {
      char tmp[BFAM_BUFSIZ];
      sprintf(tmp,"%d_%d_%d",rank,rank,size+1);
      bfam_subdomain_comm_test_t* newSub =
        bfam_malloc(sizeof(bfam_subdomain_comm_test_t));
      bfam_subdomain_init((bfam_subdomain_t*) newSub,size+1,tmp);
      newSub->ns = size;
      newSub->np = 0;
      bfam_subdomain_add_tag((bfam_subdomain_t*) newSub,"_glue");
      bfam_domain_add_subdomain(&domain,(bfam_subdomain_t*) newSub);
    }
  }
  else
  {
    char tmp[BFAM_BUFSIZ];
    sprintf(tmp,"%d_%d",rank,0);
    bfam_subdomain_comm_test_t* newSub =
      bfam_malloc(sizeof(bfam_subdomain_comm_test_t));
    bfam_subdomain_init((bfam_subdomain_t*) newSub,10*rank,tmp);
    newSub->ns = rank;
    newSub->np =    0;
    bfam_subdomain_add_tag((bfam_subdomain_t*) newSub,"_glue");
    bfam_domain_add_subdomain(&domain,(bfam_subdomain_t*) newSub);
  }

  /* set up communicator */
  const char* tags[] = {"_glue",NULL};
  bfam_communicator_t* communicator = bfam_communicator_new(&domain,
      BFAM_DOMAIN_AND,tags,MPI_COMM_WORLD);

  bfam_domain_free(&domain);
  BFAM_MPI_CHECK(MPI_Finalize());
}
