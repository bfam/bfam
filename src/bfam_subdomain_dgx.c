#include <bfam_subdomain_dgx.h>
#include <bfam_jacobi.h>
#include <bfam_kron.h>
#include <bfam_log.h>
#include <bfam_util.h>
#include <bfam_vtk.h>

#ifndef BFAM_DGX_DIMENSION
#define BFAM_DGX_DIMENSION
#define USE_GENERIC_DGX_DIMENSION
#else
#define DIM (BFAM_DGX_DIMENSION)
#endif

/* calculate the integer power of a function using exponential by squaring
 * method, taken from
 * http://stackoverflow.com/questions/101439/the-most-efficient-way-to-implement-an-integer-based-power-function-powint-int
 * accessed on 11/19/13
 */
static inline int
ipow(int base, int exp)
{
  BFAM_ASSERT(exp >= 0);

  int result = 1;
  while (exp)
  {
    /* check is odd */
    if (exp & 1)
      result *= base;

    /* divide by 2 */
    exp >>= 1;

    base *= base;
  }

  return result;
}


void
BFAM_APPEND_EXPAND(bfam_subdomain_dgx_init_,BFAM_DGX_DIMENSION)(
                              bfam_subdomain_dgx_t *subdomain,
                        const bfam_locidx_t         id,
                        const char                 *name,
                        const int                   N,
                        const bfam_locidx_t         Nv,
                        const int                   num_Vi,
                        const bfam_long_real_t    **Vi,
                        const bfam_locidx_t         K,
                        const bfam_locidx_t        *EToV,
                        const bfam_locidx_t        *EToE,
                        const int8_t               *EToF,
                        const int                   inDIM)
{
#ifdef USE_GENERIC_DGX_DIMENSION
  BFAM_WARNING("Using generic bfam_subdomain_dgx_init");
  const int DIM = inDIM;
#endif
  BFAM_ASSERT(DIM == inDIM);
  BFAM_ABORT_IF(DIM < 0, "dimension %d is not possible in bfam",DIM);
  BFAM_ABORT_IF(DIM == 0 && N != 0,
                "if DIM < 1 then N must be zero (i.e., constant");

  bfam_subdomain_init(&subdomain->base, id, name);
  bfam_subdomain_add_tag(&subdomain->base, "_subdomain_dgx");
  char dim_str[BFAM_BUFSIZ];
  snprintf(dim_str,BFAM_BUFSIZ,"_dimension_%d",DIM);
  bfam_subdomain_add_tag(&subdomain->base, dim_str);
  subdomain->dim = DIM;

  subdomain->base.free =
              BFAM_APPEND_EXPAND(bfam_subdomain_dgx_free_,BFAM_DGX_DIMENSION);
  // subdomain->base.vtk_write_vtu_piece =
  //   bfam_subdomain_dgx_quad_vtk_write_vtu_piece;
  // subdomain->base.field_add = bfam_subdomain_dgx_quad_field_add;
  // subdomain->base.field_face_add = bfam_subdomain_dgx_quad_field_face_add;
  // subdomain->base.field_init = bfam_subdomain_dgx_quad_field_init;

  subdomain->numg = DIM;
  const int numg = subdomain->numg;

  int *Ng  = NULL;
  if(numg > 0) Ng  = bfam_malloc_aligned(sizeof(int)*numg);
  subdomain->Ng = Ng;

  int *Ngp  = NULL;
  if(numg > 0) Ngp = bfam_malloc_aligned(sizeof(int)*numg);
  subdomain->Ngp = Ngp;

  const int Nrp = N+1;

  /* this could probably be made generic for arbitrary dimensions, but until
   * that's needed... */
  switch(DIM)
  {
    case 0:
      subdomain->Np = 1;
      break;

    case 1:
      subdomain->Np = Nrp;

      /* just corners */
      Ng [0]  = 2;
      Ngp[0]  = 1;
      break;

    case 2:
      subdomain->Np = Nrp*Nrp;

      /* edges */
      Ng [0]  = 4;
      Ngp[0]  = Nrp;

      /* corners */
      Ng [1]  = 4;
      Ngp[1]  = 1;
      break;

    case 3:
      subdomain->Np = Nrp*Nrp*Nrp;

      /* faces */
      Ng [0]  = 8;
      Ngp[0]  = Nrp*Nrp;

      /* edges */
      Ng [1]  = 12;
      Ngp[1]  = Nrp;

      /* corners */
      Ng [2]  = 8;
      Ngp[2]  = 1;
      break;

    default:
      BFAM_ABORT("cannot handle dim = %d",DIM);
  }

  const int Np = subdomain->Np;

  subdomain->gmask = NULL;
  if(numg > 0) subdomain->gmask = bfam_malloc_aligned(numg * sizeof(int**));
  for(int g = 0; g < numg; g++)
  {
    subdomain->gmask[g] = bfam_malloc_aligned(Ng[g] * sizeof(int*));
    for(int i = 0; i < Ng[g]; i++)
      subdomain->gmask[g][i] = bfam_malloc_aligned(Ngp[g] * sizeof(int));
  }
}

bfam_subdomain_dgx_t*
BFAM_APPEND_EXPAND(bfam_subdomain_dgx_new_,BFAM_DGX_DIMENSION)(
                       const bfam_locidx_t      id,
                       const char              *name,
                       const int                N,
                       const bfam_locidx_t      Nv,
                       const int                num_Vi,
                       const bfam_long_real_t **Vi,
                       const bfam_locidx_t      K,
                       const bfam_locidx_t     *EToV,
                       const bfam_locidx_t     *EToE,
                       const int8_t            *EToF,
                       const int                inDIM)
{
#ifdef USE_GENERIC_DGX_DIMENSION
  BFAM_WARNING("Using generic bfam_subdomain_dgx_new");
  const int DIM = inDIM;
#endif
  BFAM_ASSERT(DIM == inDIM);

  bfam_subdomain_dgx_t* newSubdomain =
    bfam_malloc(sizeof(bfam_subdomain_dgx_t));

  BFAM_APPEND_EXPAND(bfam_subdomain_dgx_init_,BFAM_DGX_DIMENSION)(
                         newSubdomain, id, name, N, Nv, num_Vi, Vi, K, EToV,
                         EToE, EToF,DIM);
  return newSubdomain;
}

void
BFAM_APPEND_EXPAND(bfam_subdomain_dgx_free_,BFAM_DGX_DIMENSION)(
    bfam_subdomain_t *thisSubdomain)
{
  bfam_subdomain_dgx_t *sub = (bfam_subdomain_dgx_t*) thisSubdomain;

  // bfam_dictionary_allprefixed_ptr(&sub->base.fields,"",
  //     &bfam_subdomain_dgx_free_fields,NULL);
  // bfam_dictionary_allprefixed_ptr(&sub->base.fields_p,"",
  //     &bfam_subdomain_dgx_free_fields,NULL);
  // bfam_dictionary_allprefixed_ptr(&sub->base.fields_m,"",
  //     &bfam_subdomain_dgx_free_fields,NULL);
  // bfam_dictionary_allprefixed_ptr(&sub->base.fields_face,"",
  //     &bfam_subdomain_dgx_free_fields,NULL);

  bfam_subdomain_free(thisSubdomain);

  // bfam_free_aligned(sub->Dr);
  // bfam_free_aligned(sub->V);

  // bfam_free_aligned(sub->r);
  // bfam_free_aligned(sub->w);
  // bfam_free_aligned(sub->wi);

  // for(int f = 0; f < sub->Nfaces; ++f)
  //    bfam_free_aligned(sub->fmask[f]);
  // bfam_free_aligned(sub->fmask);

  if(sub->gmask)
  {
    for(int g = 0; g < sub->numg; g++)
    {
      for(int i = 0; i < sub->Ng[g]; i++)
        bfam_free_aligned(sub->gmask[g][i]);
      bfam_free_aligned(sub->gmask[g]);
    }
    bfam_free_aligned(sub->gmask);
    sub->gmask = NULL;
  }
  if(sub->Ng ) bfam_free_aligned(sub->Ng ); sub->Ng  = NULL;
  if(sub->Ngp) bfam_free_aligned(sub->Ngp); sub->Ngp = NULL;

  // bfam_free_aligned(sub->vmapP);
  // bfam_free_aligned(sub->vmapM);
}
