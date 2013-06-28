#include <bfam.h>
#include <bfam_subdomain_sbp.h>

void
simple_partition_1d(bfam_locidx_t *Nl, bfam_gloidx_t *gx, bfam_locidx_t *Nb,
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
simple_partition(bfam_locidx_t *Nl, bfam_gloidx_t *gx, bfam_locidx_t *Nb,
    const bfam_gloidx_t *N, bfam_locidx_t size_in, bfam_locidx_t rank,
    bfam_locidx_t bufsz,int dim)
{
  int pd[dim];
  for(int d = 0; d < dim; d++) pd[d] = 0;
  BFAM_MPI_CHECK(MPI_Dims_create(size_in, dim, pd));

  bfam_locidx_t size = size_in;
  if(dim == 3)
  {
    simple_partition_1d(&Nl[2],&gx[2],&Nb[2*2],N[2],pd[2],
        rank/(pd[0]*pd[1]),bufsz);
    size = size % (pd[0]*pd[1]);
  }
  simple_partition_1d(&Nl[1],&gx[1],&Nb[1*2],N[1],pd[1],rank/pd[0],bufsz);
  simple_partition_1d(&Nl[0],&gx[0],&Nb[0*2],N[0],pd[0],rank%pd[0],bufsz);
}

void
simple_load_balance(bfam_locidx_t *procs,const bfam_gloidx_t *N,
    int mpi_size,int num_blocks,int dim)
{
  if(mpi_size <= num_blocks)
  {
    for(int b = 0; b < num_blocks; b++)
    {
      procs[2*b  ] = 0;      // starting rank
      procs[2*b+1] = mpi_size-1; // ending rank
    }
  }
  else
  {
    /* get cost per block */
    bfam_real_t cost[num_blocks];
    bfam_locidx_t num_procs[num_blocks];
    for(int b = 0; b < num_blocks; b++)
    {
      num_procs[b] = 1;
      cost[b] = 1;
      for(int d = 0;d < dim;d++) cost[b] *= N[dim*b+1];
    }

    /* figured out number of procs per block */
    for(int p = num_blocks; p < mpi_size; p++)
    {
      bfam_locidx_t max_blck = 0;
      bfam_real_t   max_cost = cost[max_blck]/num_procs[max_blck];
      for(int b = 1; b < num_blocks; b++)
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

    bfam_locidx_t lp = -1;
    for(int b = 0; b < num_blocks; b++)
    {
      procs[2*b  ] = lp+1;
      procs[2*b+1] = lp+num_procs[b];
      lp           = procs[2*b+1];
    }
  }
}

int
main (int argc, char *argv[])
{
  int rank, mpi_size;

  /* startup MPI */
  BFAM_MPI_CHECK(MPI_Init(&argc,&argv));
  BFAM_MPI_CHECK(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));
  BFAM_MPI_CHECK(MPI_Comm_rank(MPI_COMM_WORLD, &rank));
  bfam_log_init(rank,stdout,BFAM_LL_VERBOSE);

  /* setup the domain*/
  bfam_domain_t domain;
  bfam_domain_init(&domain,MPI_COMM_WORLD);

  /* connectivity */
  int num_blocks = 3;
  const char *names[] = {"sub0","sub1","sub2"};
  bfam_locidx_t bufsz = 5;

  int foo   = 15;
  int dim   =  2;
  bfam_gloidx_t   EtoV[] = {0,1,3,4,
                            1,2,4,5,
                            3,4,6,5};
  bfam_gloidx_t   EtoE[] = {0,1,0,2,
                            0,1,1,2,
                            2,1,0,2};
  bfam_locidx_t   EtoF[] = {0,0,2,2,
                            1,1,2,1,
                            2,3,3,3};
  bfam_locidx_t   Eto0[] = {0,0,2,2,
                            1,1,2,1,
                            2,3,3,3};
  bfam_gloidx_t      N[] = {1*foo,2*foo,
                            3*foo,2*foo,
                            1*foo,3*foo};
  bfam_long_real_t  Vx[] = {0,0.25,1,  0,0.25,0.5,0};
  bfam_long_real_t  Vy[] = {0,   0,0,0.5,0.25,0.5,1};
  bfam_long_real_t *Vz   = NULL;
  bfam_locidx_t  procs[2*num_blocks];

  simple_load_balance(procs,N,mpi_size,num_blocks,dim);

  for(int b = 0; b < num_blocks;b++)
  {
    if(procs[2*b] <= rank && rank <= procs[2*b+1])
    {
      bfam_locidx_t Nl[2] = {-1,-1};
      bfam_gloidx_t gx[2] = {-1,-1};
      bfam_locidx_t Nb[4] = {0,0,0,0};
      bfam_locidx_t l_rank = rank-procs[2*b];
      bfam_locidx_t l_size = procs[2*b+1]-procs[2*b]+1;

      simple_partition(Nl,gx,Nb,&N[2*b],l_size,l_rank,bufsz,dim);
      BFAM_VERBOSE("s%-3d: Nl: %3d %3d gx: %3d %3d buf: %3d %3d %3d %3d",
          b,Nl[0],Nl[1],gx[0],gx[1],(int)Nb[0],(int)Nb[1],(int)Nb[2],(int) Nb[3]);

      bfam_gloidx_t V[] = {EtoV[4*b+0],EtoV[4*b+1],EtoV[4*b+2],EtoV[4*b+3]};
      bfam_long_real_t  x[] = {Vx[V[0]],Vx[V[1]],Vx[V[2]],Vx[V[3]]};
      bfam_long_real_t  y[] = {Vy[V[0]],Vy[V[1]],Vy[V[2]],Vy[V[3]]};
      bfam_long_real_t *z   = Vz;

      bfam_subdomain_sbp_t *sub =
        bfam_subdomain_sbp_new(0,l_rank,l_size,names[b],dim,
            &N[2*b],Nl,Nb,gx,x,y,z);

      bfam_domain_add_subdomain(&domain,(bfam_subdomain_t*)sub);
    }
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
