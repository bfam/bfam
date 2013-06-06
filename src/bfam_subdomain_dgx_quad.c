#include <bfam_subdomain_dgx_quad.h>
#include <bfam_jacobi.h>

bfam_subdomain_dgx_quad_t*
bfam_subdomain_dgx_quad_new(const char             *name,
                            const int               N,
                            const bfam_locidx_t     Nv,
                            const bfam_long_real_t *VX,
                            const bfam_long_real_t *VY,
                            const bfam_long_real_t *VZ,
                            const bfam_locidx_t     K,
                            const bfam_locidx_t    *EToV,
                            const bfam_locidx_t    *EToE,
                            const int8_t           *EToF)
{
  bfam_subdomain_dgx_quad_t* newSubdomain =
    bfam_malloc(sizeof(bfam_subdomain_dgx_quad_t));

  bfam_subdomain_dgx_quad_init(newSubdomain, name, N, Nv, VX, VY, VZ, K, EToV,
                               EToE, EToF);

  return newSubdomain;
}

void
bfam_subdomain_dgx_quad_init(bfam_subdomain_dgx_quad_t       *subdomain,
                             const char                      *name,
                             const int                        N,
                             const bfam_locidx_t              Nv,
                             const bfam_long_real_t          *VX,
                             const bfam_long_real_t          *VY,
                             const bfam_long_real_t          *VZ,
                             const bfam_locidx_t              K,
                             const bfam_locidx_t             *EToV,
                             const bfam_locidx_t             *EToE,
                             const int8_t                    *EToF)
{
  bfam_subdomain_init(&subdomain->base, name);
  bfam_subdomain_add_tag(&subdomain->base, "_subdomain_dgx_quad");
  subdomain->base.free = bfam_subdomain_dgx_quad_free;

  const int Np = (N+1)*(N+1);
  const int Nfp = N+1;
  const int Nfaces = 4;
  const int Ncorners = 4;

  const int Nrp = N+1;
  bfam_long_real_t *lr, *lw;
  lr = bfam_malloc_aligned(Nrp*sizeof(bfam_long_real_t));
  lw = bfam_malloc_aligned(Nrp*sizeof(bfam_long_real_t));

  bfam_jacobi_gauss_lobatto_quadrature(0, 0, N, lr, lw);

  bfam_long_real_t *lx, *ly, *lz;
  lx = bfam_malloc_aligned(K*Np*sizeof(bfam_long_real_t));
  ly = bfam_malloc_aligned(K*Np*sizeof(bfam_long_real_t));
  lz = bfam_malloc_aligned(K*Np*sizeof(bfam_long_real_t));


  for(bfam_locidx_t k = 0; k < K; ++k)
  {
    bfam_locidx_t va = EToV[Ncorners * k + 0];
    bfam_locidx_t vb = EToV[Ncorners * k + 1];
    bfam_locidx_t vc = EToV[Ncorners * k + 2];
    bfam_locidx_t vd = EToV[Ncorners * k + 3];

    for(int n = 0; n < Nrp; ++n)
    {
      for(int m = 0; m < Nrp; ++m)
      {
        bfam_long_real_t wa, wb, wc, wd;

        wa = (1-lr[m])*(1-lr[n]);
        wb = (1+lr[m])*(1-lr[n]);
        wc = (1-lr[m])*(1+lr[n]);
        wd = (1+lr[m])*(1+lr[n]);

        lx[Np*k + Nrp*n + m] = (wa*VX[va]+wb*VX[vb]+wc*VX[vc]+wd*VX[vd])/4;
        ly[Np*k + Nrp*n + m] = (wa*VY[va]+wb*VY[vb]+wc*VY[vc]+wd*VY[vd])/4;
        lz[Np*k + Nrp*n + m] = (wa*VZ[va]+wb*VZ[vb]+wc*VZ[vc]+wd*VZ[vd])/4;
      }
    }
  }

  /*
   * Set subdomain values
   */
  subdomain->N = N;
  subdomain->Np = Np;
  subdomain->Nfp = Nfp;
  subdomain->Nfaces = Nfaces;
  subdomain->Ncorners = Ncorners;

  subdomain->r = bfam_malloc_aligned(Nrp*sizeof(bfam_real_t));
  subdomain->w = bfam_malloc_aligned(Nrp*sizeof(bfam_real_t));
  for(int n = 0; n<Nrp; ++n)
  {
    subdomain->r[n] = (bfam_real_t) lr[n];
    subdomain->w[n] = (bfam_real_t) lw[n];
  }

  subdomain->K = K;

  subdomain->x = bfam_malloc_aligned(K*Np*sizeof(bfam_real_t));
  subdomain->y = bfam_malloc_aligned(K*Np*sizeof(bfam_real_t));
  subdomain->z = bfam_malloc_aligned(K*Np*sizeof(bfam_real_t));
  for(int n = 0; n < K*Np; ++n)
  {
    subdomain->x[n] = (bfam_real_t) lx[n];
    subdomain->y[n] = (bfam_real_t) ly[n];
    subdomain->z[n] = (bfam_real_t) lz[n];
  }

  bfam_free_aligned(lr);
  bfam_free_aligned(lw);

  bfam_free_aligned(lx);
  bfam_free_aligned(ly);
  bfam_free_aligned(lz);
}

void
bfam_subdomain_dgx_quad_free(bfam_subdomain_t *thisSubdomain)
{
  bfam_subdomain_dgx_quad_t *sub = (bfam_subdomain_dgx_quad_t*) thisSubdomain;

  bfam_subdomain_free(thisSubdomain);

  bfam_free_aligned(sub->r);
  bfam_free_aligned(sub->w);

  bfam_free_aligned(sub->x);
  bfam_free_aligned(sub->y);
  bfam_free_aligned(sub->z);
}
