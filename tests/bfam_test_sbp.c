#include <bfam.h>
#include <bfam_subdomain_sbp.h>

void
simple_load_balance_1d(bfam_locidx_t *Nl, bfam_gloidx_t *gx, bfam_locidx_t *Nb,
    bfam_gloidx_t N, bfam_locidx_t size, bfam_locidx_t rank,
    bfam_locidx_t bufsz)
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
  if(rank == 0) Nb[0] = 0;
  else Nb[0] = bufsz;
  if(rank == size-1) Nb[1] = 0;
  else Nb[1] = bufsz;
}

void
simple_load_balance(bfam_locidx_t *Nl, bfam_gloidx_t *gx, bfam_locidx_t *Nb,
    const bfam_gloidx_t *N, const bfam_locidx_t *pd, bfam_locidx_t rank,
    bfam_locidx_t bufsz)
{
  simple_load_balance_1d(Nl,gx,Nb,N[0],pd[0],rank%pd[0],bufsz);
  simple_load_balance_1d(&Nl[1],&gx[1],&Nb[2],N[1],pd[1],rank/pd[0],bufsz);
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
  int foo   = 15;
  int bufsz = 9;
  bfam_long_real_t x0[4] = {0,0.25,0,0.25};
  bfam_long_real_t y0[4] = {0,0,0.5,0.25};
  bfam_gloidx_t    N0[2] = {foo*1,foo*2};

  bfam_long_real_t x1[4] = {0.25,1,0.25,0.5};
  bfam_long_real_t y1[4] = {0,0,0.25,0.5};
  bfam_gloidx_t    N1[2] = {foo*3,foo*2};

  bfam_long_real_t x2[4] = {0.25,0.5,0,0};
  bfam_long_real_t y2[4] = {0.25,0.5,0.5,1};
  bfam_gloidx_t    N2[2] = {foo*3,foo*1};


  int nb = 3;
  bfam_long_real_t *x[] = {x0,x1,x2};
  bfam_long_real_t *y[] = {y0,y1,y2};
  bfam_gloidx_t    *N[] = {N0,N1,N2};
  const char *names[] = {"sub0","sub1","sub2"};

  /* figure out which part of the domain I handle */
  if(size <= nb)
  {
    for(int b = 0; b < nb; b++)
    {
      bfam_locidx_t Nl[2] = {-1,-1};
      bfam_gloidx_t gx[2] = {-1,-1};
      bfam_locidx_t Nb[4] = {0,0,0,0};
      bfam_locidx_t pd[2] = {1,size};
      simple_load_balance(Nl,gx,Nb,N[b],pd,rank,bufsz);
      BFAM_VERBOSE("s%d: Nl: %d %d gx: %d %d %d %d %d %d",
          b,Nl[0],Nl[1],gx[0],gx[1],(int)Nb[0],(int)Nb[1],(int)Nb[2],(int) Nb[3]);
      bfam_subdomain_sbp_t *sub =
        bfam_subdomain_sbp_new(0,rank,size,names[b],2,
            N[b],Nl,Nb,gx,x[b],y[b],NULL);
      bfam_domain_add_subdomain(&domain,(bfam_subdomain_t*)sub);
    }
  }
  else
  {
    /* get cost per block */
    bfam_real_t cost[nb];
    bfam_locidx_t num_procs[nb];
    for(int b = 0; b < nb; b++)
    {
      num_procs[b] = 1;
      cost[b]      = N[b][0]*N[b][1];
    }

    /* figured out number of procs per block */
    for(int p = nb; p < size; p++)
    {
      bfam_locidx_t max_blck = 0;
      bfam_real_t   max_cost = cost[max_blck]/num_procs[max_blck];
      for(int b = 1; b < nb; b++)
      {
        bfam_real_t test = cost[b]/num_procs[b];
        if(test > max_cost)
        {
          max_cost = test;
          max_blck = b;
        }
      }
      num_procs[max_blck]++;
    }

    /* figured out which block I handle and which piece of that block */
    bfam_locidx_t blck = 0; // my block
    bfam_locidx_t tp = 0;   // procs before my block
    for(;blck < nb-1;blck++)
      if(rank < tp+num_procs[blck]) break;
      else tp += num_procs[blck];

    /* determine block topology */
    int mpi_dims[2] = {0,0};
    BFAM_MPI_CHECK(MPI_Dims_create(num_procs[blck], 2, mpi_dims));

    /* load balance */
    bfam_locidx_t Nl[2] = {-1,-1};
    bfam_gloidx_t gx[2] = {-1,-1};
    bfam_locidx_t Nb[4] = {0,0,0,0};
    simple_load_balance(Nl,gx,Nb,N[blck],mpi_dims,rank-tp,bufsz);
    BFAM_VERBOSE("s%d: Nl: %d %d gx: %d %d %d %d %d %d",
        blck,Nl[0],Nl[1],gx[0],gx[1],
        (int)Nb[0],(int)Nb[1],(int)Nb[2],(int)Nb[3]);
    bfam_subdomain_sbp_t *sub =
      bfam_subdomain_sbp_new(0,rank-tp,num_procs[blck],names[blck],2,
          N[blck],Nl,Nb,gx,x[blck],y[blck],NULL);
    bfam_domain_add_subdomain(&domain,(bfam_subdomain_t*)sub);

  }

  const char *tags[] = {NULL};
  const char *scalars[] = {"_grid_x","_grid_y","_grid_z",NULL};
  const char *vectors[] = {NULL};
  const char *components[] = {NULL};
  bfam_vtk_write_struc_file(&domain,BFAM_DOMAIN_AND,
      tags,"sbp_fields",scalars,vectors,components,0,1);


  /* clean up */
  bfam_domain_free(&domain);
  BFAM_MPI_CHECK(MPI_Finalize());
}
