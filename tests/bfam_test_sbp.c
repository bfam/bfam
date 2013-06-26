#include <bfam.h>
#include <bfam_subdomain_sbp.h>

void
simple_load_balance_1d(bfam_locidx_t *Nl, bfam_gloidx_t *gx,
    bfam_gloidx_t N, bfam_locidx_t size, bfam_locidx_t rank)
{
  Nl[0] = (N+1)/size;
  gx[0] = Nl[0]*rank;

  bfam_locidx_t rem = (N+1) % size;
  if(rank >= rem)
  {
    Nl[0] -= 1;
    gx[0] += rem;
  }
  else
  {
    gx[0] += rank;
  }
}

void
simple_load_balance(bfam_locidx_t *Nl, bfam_gloidx_t *gx,
    const bfam_gloidx_t *N, const bfam_locidx_t *pd, bfam_locidx_t rank)
{
  simple_load_balance_1d(Nl,gx,N[0],pd[0],rank%pd[0]);
  simple_load_balance_1d(&Nl[1],&gx[1],N[1],pd[1],rank/pd[0]);
}

int
main (int argc, char *argv[])
{
  int rank, size;

  /* startup MPI */
  BFAM_MPI_CHECK(MPI_Init(&argc,&argv));
  BFAM_MPI_CHECK(MPI_Comm_size(MPI_COMM_WORLD, &size));
  BFAM_MPI_CHECK(MPI_Comm_rank(MPI_COMM_WORLD, &rank));
  bfam_log_init(rank,stdout,BFAM_LL_VERBOSE);

  /* setup the domain*/
  bfam_domain_t domain;
  bfam_domain_init(&domain,MPI_COMM_WORLD);

  /* set up the block information */
  bfam_long_real_t x0[4] = {0,0.25,0,0.25};
  bfam_long_real_t y0[4] = {0,0,0.5,0.25};
  bfam_gloidx_t    N0[2] = {100,200};

  bfam_long_real_t x1[4] = {0.25,1,0.25,0.5};
  bfam_long_real_t y1[4] = {0,0,0.25,0.5};
  bfam_gloidx_t    N1[2] = {300,200};

  bfam_long_real_t x2[4] = {0.25,0.5,0,0};
  bfam_long_real_t y2[4] = {0.25,0.5,0.5,1};
  bfam_gloidx_t    N2[2] = {300,100};

  /* figure out which part of the domain I handle */
  if(size <= 6)
  {
    bfam_locidx_t Nl[2] = {-1,-1};
    bfam_gloidx_t gx[2] = {-1,-1};
    bfam_locidx_t Nb[4] = {0,0,0,0};

    // Subdomain 0
    bfam_locidx_t pd[2] = {1,size};
    simple_load_balance(Nl,gx,N0,pd,rank);
    BFAM_VERBOSE("s0: Nl: %d %d gx: %d %d",Nl[0],Nl[1],gx[0],gx[1]);
    bfam_subdomain_sbp_t *sub0 =
      bfam_subdomain_sbp_new(0,rank,size,"sub0",2,N0,Nl,Nb,gx,x0,y0,NULL);
    bfam_domain_add_subdomain(&domain,(bfam_subdomain_t*)sub0);

    // Subdomain 1
    pd[0] = size; pd[1] = 1;
    simple_load_balance(Nl,gx,N1,pd,rank);
    BFAM_VERBOSE("s1: Nl: %d %d gx: %d %d",Nl[0],Nl[1],gx[0],gx[1]);
    bfam_subdomain_sbp_t *sub1 =
      bfam_subdomain_sbp_new(1,rank,size,"sub1",2,N1,Nl,Nb,gx,x1,y1,NULL);
    bfam_domain_add_subdomain(&domain,(bfam_subdomain_t*)sub1);

    // Subdomain 2
    pd[0] = size; pd[1] = 1;
    simple_load_balance(Nl,gx,N2,pd,rank);
    BFAM_VERBOSE("s2: Nl: %d %d gx: %d %d",Nl[0],Nl[1],gx[0],gx[1]);
    bfam_subdomain_sbp_t *sub2 =
      bfam_subdomain_sbp_new(2,rank,size,"sub2",2,N2,Nl,Nb,gx,x2,y2,NULL);
    bfam_domain_add_subdomain(&domain,(bfam_subdomain_t*)sub2);

  }
  else
  {
  }

  /* clean up */
  bfam_domain_free(&domain);
  BFAM_MPI_CHECK(MPI_Finalize());
}
