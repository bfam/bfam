#include <bfam.h>
#include <bfam_subdomain_sbp.h>

static void
poly1_field(bfam_locidx_t npoints, bfam_real_t time, bfam_real_t *restrict x,
    bfam_real_t *restrict y, bfam_real_t *restrict z,
    struct bfam_subdomain *s, void *arg, bfam_real_t *restrict field)
{
  BFAM_ASSUME_ALIGNED(x, 32);
  BFAM_ASSUME_ALIGNED(y, 32);
  BFAM_ASSUME_ALIGNED(z, 32);
  BFAM_ASSUME_ALIGNED(field, 32);

  if(z == NULL)
    for(bfam_locidx_t n=0; n < npoints; ++n)
      field[n] = -x[n] - y[n];
  else
    for(bfam_locidx_t n=0; n < npoints; ++n)
      field[n] = -x[n] - y[n] - z[n];
}

static void
poly2_field(bfam_locidx_t npoints, bfam_real_t time, bfam_real_t *restrict x,
    bfam_real_t *restrict y, bfam_real_t *restrict z,
    struct bfam_subdomain *s, void *arg, bfam_real_t *restrict field)
{
  BFAM_ASSUME_ALIGNED(x, 32);
  BFAM_ASSUME_ALIGNED(y, 32);
  BFAM_ASSUME_ALIGNED(z, 32);
  BFAM_ASSUME_ALIGNED(field, 32);

  if(z == NULL)
    for(bfam_locidx_t n=0; n < npoints; ++n)
      field[n] = x[n] + y[n];
  else
    for(bfam_locidx_t n=0; n < npoints; ++n)
      field[n] = x[n] + y[n] + z[n];
}

static void
poly3_field(bfam_locidx_t npoints, bfam_real_t time, bfam_real_t *restrict x,
    bfam_real_t *restrict y, bfam_real_t *restrict z,
    struct bfam_subdomain *s, void *arg, bfam_real_t *restrict field)
{
  BFAM_ASSUME_ALIGNED(x, 32);
  BFAM_ASSUME_ALIGNED(y, 32);
  BFAM_ASSUME_ALIGNED(z, 32);
  BFAM_ASSUME_ALIGNED(field, 32);

  BFAM_ASSERT(z != NULL);
  for(bfam_locidx_t n=0; n < npoints; ++n)
    field[n] = x[n] + y[n] - z[n];
}


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


void
setup_subdomains(bfam_domain_t *domain,
    int dim, int num_blocks, bfam_locidx_t bufsz,int rank,const char **names,
    const bfam_locidx_t *procs,const bfam_gloidx_t *EtoV,
    const bfam_gloidx_t *N, const bfam_long_real_t* Vx,
    const bfam_long_real_t *Vy, const bfam_long_real_t *Vz)
{
  int ncorn = (int)pow(2,dim);
  bfam_long_real_t *x = NULL;
  bfam_long_real_t *y = NULL;
  bfam_long_real_t *z = NULL;
  if(dim > 0) x = bfam_malloc_aligned(ncorn*sizeof(bfam_long_real_t));
  if(dim > 1) y = bfam_malloc_aligned(ncorn*sizeof(bfam_long_real_t));
  if(dim > 2) z = bfam_malloc_aligned(ncorn*sizeof(bfam_long_real_t));
  for(int b = 0; b < num_blocks;b++)
  {
    if(procs[2*b] <= rank && rank <= procs[2*b+1])
    {
      bfam_locidx_t Nl[  dim];
      bfam_gloidx_t gx[  dim];
      bfam_locidx_t Nb[2*dim];
      for(int d = 0;d<dim;d++)
      {
        Nl[d] = -1;
        gx[d] = -1;
        Nb[2*d] = 0;
        Nb[2*d+1] = 0;
      }

      bfam_locidx_t l_rank = rank-procs[2*b];
      bfam_locidx_t l_size = procs[2*b+1]-procs[2*b]+1;

      simple_partition(Nl,gx,Nb,&N[dim*b],l_size,l_rank,bufsz,dim);

      for(int ix = 0; ix < ncorn; ix++)
      {
        int V = EtoV[ncorn*b+ix];
        if(x != NULL) x[ix] = Vx[V];
        if(y != NULL) y[ix] = Vy[V];
        if(z != NULL) z[ix] = Vz[V];
      }

      bfam_subdomain_sbp_t *sub =
        bfam_subdomain_sbp_new(b,l_rank,l_size,names[b],dim,
            &N[dim*b],Nl,Nb,gx,x,y,z);

      bfam_subdomain_add_tag((bfam_subdomain_t*)sub,"_volume");

      bfam_domain_add_subdomain(domain,(bfam_subdomain_t*)sub);
    }
  }
  if(x != NULL) bfam_free_aligned(x);
  if(y != NULL) bfam_free_aligned(y);
  if(z != NULL) bfam_free_aligned(z);
}

void
test_2d(int rank, int mpi_size)
{
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
  // bfam_gloidx_t   EtoE[] = {0,1,0,2,
  //                           0,1,1,2,
  //                           2,1,0,2};
  // bfam_locidx_t   EtoF[] = {0,0,2,2,
  //                           1,1,2,1,
  //                           2,3,3,3};
  bfam_gloidx_t      N[] = {1*foo,2*foo,
                            3*foo,2*foo,
                            1*foo,3*foo};
  bfam_long_real_t  Vx[] = {0,0.25,1,  0,0.25,0.5,0};
  bfam_long_real_t  Vy[] = {0,   0,0,0.5,0.25,0.5,1};
  bfam_long_real_t *Vz   = NULL;
  bfam_locidx_t  procs[2*num_blocks];

  /* load balance */
  simple_load_balance(procs,N,mpi_size,num_blocks,dim);

  /* setup subdomains */
  setup_subdomains(&domain,dim,num_blocks,bufsz,rank,
      names,procs,EtoV,N,Vx,Vy,Vz);

  /* put in some fields */
  const char *volume[] = {"_volume", NULL};
  bfam_domain_add_field(&domain, BFAM_DOMAIN_OR, volume, "p1");
  bfam_domain_add_field(&domain, BFAM_DOMAIN_OR, volume, "p2");

  bfam_domain_init_field(&domain, BFAM_DOMAIN_OR, volume, "p1",
      0, poly1_field, NULL);
  bfam_domain_init_field(&domain, BFAM_DOMAIN_OR, volume, "p2",
      0, poly2_field, NULL);

  /* dump the entire mesh */
  const char *tags[] = {NULL};
  const char *scalars[] = {"p1","p2",NULL};
  const char *vectors[] = {NULL};
  const char *components[] = {NULL};
  bfam_vtk_write_struc_file(&domain,BFAM_DOMAIN_AND,
      tags,"sbp_fields_2d",scalars,vectors,components,0,1);

  /* clean up */
  bfam_domain_free(&domain);
}

void
test_3d(int rank, int mpi_size)
{
  /* setup the domain*/
  bfam_domain_t domain;
  bfam_domain_init(&domain,MPI_COMM_WORLD);

  int foo   = 15;
  int dim   =  3;

  const char *names[] = {"sub0","sub1","sub2","sub3","sub4","sub5"};
  bfam_locidx_t bufsz = 5;

  /* connectivity from p8est_connectivity_new_rotcubes */
  bfam_locidx_t num_blocks = 6;
  const bfam_gloidx_t EtoV[6 * 8] = {
    0, 17, 3, 4, 15, 11, 13, 14,
    7, 2, 6, 17, 9, 12, 8, 11,
    2, 12, 5, 10, 17, 11, 4, 14,
    19, 13, 18, 14, 16, 15, 1, 11,
    14, 11, 21, 25, 18, 1, 22, 23,
    21, 20, 25, 24, 14, 10, 11, 12,
  };
  const bfam_gloidx_t N[6 * 3] = {
     1+foo, 2+foo, 3+foo,
     4+foo, 5+foo, 6+foo,
     7+foo, 8+foo, 9+foo,
    10+foo,11+foo,12+foo,
    13+foo,14+foo,15+foo,
    16+foo,17+foo,18+foo
  };

  const bfam_long_real_t Vx[26] =
    {0,1,2,0,1,2,1,2,1,2,2,1,2,0,1,0,0,1,1,0,2.5,2,2,2,2.5,2};
  const bfam_long_real_t Vy[26] =
    {0,0,0,1,1,1,-1,-1,-1,-1,1,0,0,1,1,0,0,0,1,1,1.5,1.5,1.5,.5,.5,.5};
  const bfam_long_real_t Vz[26] =
    {0,2,0,0,0,0,0,0,1,1,1,1,1,1,1,1,2,0,2,2,2,2,2.5,2.5,2,2};

  bfam_locidx_t  procs[2*num_blocks];

  /* load balance */
  simple_load_balance(procs,N,mpi_size,num_blocks,dim);

  for(int b = 0; b < num_blocks;b++)
    BFAM_ROOT_VERBOSE("b%-3d: %3d %3d",b,procs[2*b],procs[2*b+1]);

  /* setup subdomains */
  setup_subdomains(&domain,dim,num_blocks,bufsz,rank,
      names,procs,EtoV,N,Vx,Vy,Vz);

  /* put in some fields */
  const char *volume[] = {"_volume", NULL};
  bfam_domain_add_field(&domain, BFAM_DOMAIN_OR, volume, "p1");
  bfam_domain_add_field(&domain, BFAM_DOMAIN_OR, volume, "p2");
  bfam_domain_add_field(&domain, BFAM_DOMAIN_OR, volume, "p3");

  bfam_domain_init_field(&domain, BFAM_DOMAIN_OR, volume, "p1",
      0, poly1_field, NULL);
  bfam_domain_init_field(&domain, BFAM_DOMAIN_OR, volume, "p2",
      0, poly2_field, NULL);
  bfam_domain_init_field(&domain, BFAM_DOMAIN_OR, volume, "p3",
      0, poly3_field, NULL);

  /* dump the entire mesh */
  const char *tags[] = {NULL};
  const char *scalars[] = {"p1","p2","p3",NULL};
  const char *vectors[] = {NULL};
  const char *components[] = {NULL};
  bfam_vtk_write_struc_file(&domain,BFAM_DOMAIN_AND,
      tags,"sbp_fields_3d",scalars,vectors,components,1,1);


  /* clean up */
  bfam_domain_free(&domain);
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

  /* test 3d */
  test_2d(rank,mpi_size);

  /* test 3d */
  test_3d(rank,mpi_size);

  /* stop MPI */
  BFAM_MPI_CHECK(MPI_Finalize());
}
