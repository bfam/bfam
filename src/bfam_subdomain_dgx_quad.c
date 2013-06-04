#include <bfam_subdomain_dgx_quad.h>
#include <bfam_jacobi.h>

bfam_subdomain_dgx_quad_t*
bfam_subdomain_dgx_quad_new(const char             *name,
                            const int               N,
                            const bfam_locidx_t     Nv,
                            const bfam_long_real_t *VX,
                            const bfam_long_real_t *VY,
                            const bfam_locidx_t     K,
                            const bfam_locidx_t    *EToV,
                            const bfam_locidx_t    *EToE,
                            const int8_t           *EToF)
{
  bfam_subdomain_dgx_quad_t* newSubdomain =
    bfam_malloc(sizeof(bfam_subdomain_dgx_quad_t));

  bfam_subdomain_dgx_quad_init(newSubdomain, name, N, Nv, VX, VY, K, EToV,
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
                             const bfam_locidx_t              K,
                             const bfam_locidx_t             *EToV,
                             const bfam_locidx_t             *EToE,
                             const int8_t                    *EToF)
{
  bfam_subdomain_init(&subdomain->base, name);

  subdomain->Np  = (N+1)*(N+1);
  subdomain->Nfp =  N+1;

  const int Nrp = N+1;
  bfam_long_real_t *lr, *lw;
  lr = bfam_malloc(Nrp*sizeof(bfam_long_real_t));
  lw = bfam_malloc(Nrp*sizeof(bfam_long_real_t));

  bfam_jacobi_gauss_lobatto_quadrature(0, 0, N, lr, lw);

  bfam_free(lr);
  bfam_free(lw);
}

void
bfam_subdomain_dgx_quad_free(bfam_subdomain_dgx_quad_t *thisSubdomain)
{
  bfam_subdomain_free(&thisSubdomain->base);
}