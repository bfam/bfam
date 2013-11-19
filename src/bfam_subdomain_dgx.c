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

  bfam_subdomain_init(&subdomain->base, id, name);
  bfam_subdomain_add_tag(&subdomain->base, "_subdomain_dgx");
  char dim_str[BFAM_BUFSIZ];
  snprintf(dim_str,BFAM_BUFSIZ,"_dimension_%d",DIM);
  bfam_subdomain_add_tag(&subdomain->base, dim_str);

  subdomain->base.free =
              BFAM_APPEND_EXPAND(bfam_subdomain_dgx_free_,BFAM_DGX_DIMENSION);
  // subdomain->base.vtk_write_vtu_piece =
  //   bfam_subdomain_dgx_quad_vtk_write_vtu_piece;
  // subdomain->base.field_add = bfam_subdomain_dgx_quad_field_add;
  // subdomain->base.field_face_add = bfam_subdomain_dgx_quad_field_face_add;
  // subdomain->base.field_init = bfam_subdomain_dgx_quad_field_init;

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
  // bfam_subdomain_dgx_t *sub = (bfam_subdomain_dgx_t*) thisSubdomain;

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
  //   bfam_free_aligned(sub->fmask[f]);
  // bfam_free_aligned(sub->fmask);

  // bfam_free_aligned(sub->vmapP);
  // bfam_free_aligned(sub->vmapM);
}
