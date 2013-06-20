#include <bfam_subdomain_sbp.h>
#include <bfam_log.h>
#include <bfam_vtk.h>

bfam_subdomain_sbp_t*
bfam_subdomain_sbp_new(const bfam_locidx_t     id,
                            const char             *name,
                            const int               dim,
                            const bfam_gloidx_t    *N,
                            const bfam_locidx_t    *Nl,
                            const bfam_locidx_t    *Nb,
                            const bfam_gloidx_t    *gx)
{
  bfam_subdomain_sbp_t *newSub = bfam_malloc(sizeof(bfam_subdomain_sbp_t));
  bfam_subdomain_sbp_init(newSub,id,name,dim,N,Nl,Nb,gx);
  return newSub;
}

void
bfam_subdomain_sbp_init(bfam_subdomain_sbp_t *subdomain,
                            const bfam_locidx_t     id,
                            const char             *name,
                            const int               dim,
                            const bfam_gloidx_t    *N,
                            const bfam_locidx_t    *Nl,
                            const bfam_locidx_t    *Nb,
                            const bfam_gloidx_t    *gx)
{
  bfam_subdomain_init(&subdomain->base, id, name);
  bfam_subdomain_add_tag(&subdomain->base, "_subdomain_sbp");
  subdomain->base.free = bfam_subdomain_sbp_free;

  subdomain->dim = dim;

  subdomain->N  = bfam_malloc(  dim*sizeof(bfam_gloidx_t));
  subdomain->Nl = bfam_malloc(  dim*sizeof(bfam_locidx_t));
  subdomain->Nb = bfam_malloc(2*dim*sizeof(bfam_locidx_t));
  subdomain->gx = bfam_malloc(  dim*sizeof(bfam_gloidx_t));

  for(int d = 0; d < dim; d++)
  {
    subdomain->N [  d  ] = N [  d  ];
    subdomain->Nl[  d  ] = Nl[  d  ];
    subdomain->Nb[2*d  ] = Nb[2*d  ];
    subdomain->Nb[2*d+1] = Nb[2*d+1];
    subdomain->gx[  d  ] = gx[  d  ];

    BFAM_ABORT_IF(N [d] < 0 || Nl[d] < 0 || Nb[2*d] < 0 || Nb[2*d+1] < 0 ||
                  gx[d] < 0 || gx[d] > N[d] || gx[d]+Nl[d] > N[d],
                  "%s: problem with a dimension %d", d, subdomain->base.name);
  }
}
