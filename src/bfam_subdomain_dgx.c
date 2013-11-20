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

static inline int***
BFAM_APPEND_EXPAND(bfam_subdomain_dgx_gmask_set_,BFAM_DGX_DIMENSION)
         (const int numg, const int N, int *Np, int *Ng, int *Ngp, int inDIM)
{
#ifdef USE_GENERIC_DGX_DIMENSION
  BFAM_WARNING("Using generic bfam_subdomain_dgx_gmask_set");
  const int DIM = inDIM;
#endif
  BFAM_ABORT_IF(DIM > 3 || DIM < 0,
      "bfam_subdomain_dgx_gmask_set cannot handle dim = %d",DIM);

  if(DIM == 0)
  {
    *Np = 1;
    return NULL;
  }


  /* this could probably be made generic for arbitrary dimensions, but until
   * that's needed... */
  switch(DIM)
  {
    case 1:
      *Np = N+1;

      /* just corners */
      Ng [0]  = 2;
      Ngp[0]  = 1;
      break;

    case 2:
      *Np = (N+1)*(N+1);

      /* edges */
      Ng [0]  = 4;
      Ngp[0]  = N+1;

      /* corners */
      Ng [1]  = 4;
      Ngp[1]  = 1;
      break;

    case 3:
      *Np = (N+1)*(N+1)*(N+1);

      /* faces */
      Ng [0]  = 6;
      Ngp[0]  = (N+1)*(N+1);

      /* edges */
      Ng [1]  = 12;
      Ngp[1]  = N+1;

      /* corners */
      Ng [2]  = 8;
      Ngp[2]  = 1;
      break;

    default:
      BFAM_ABORT("cannot handle dim = %d",DIM);
  }

  int ***gmask = bfam_malloc_aligned(numg * sizeof(int**));
  for(int g = 0; g < numg; g++)
  {
    gmask[g] = bfam_malloc_aligned(Ng[g] * sizeof(int*));
    for(int i = 0; i < Ng[g]; i++)
      gmask[g][i] = bfam_malloc_aligned(Ngp[g] * sizeof(int));
  }

  switch(DIM)
  {
    case 1:
      gmask[0][0][0] = 0;
      gmask[0][1][0] = N;
      break;

    case 2:
      /* edges */
      for(int i = 0; i < N+1; ++i) gmask[0][0][i] = i*(N+1);
      for(int i = 0; i < N+1; ++i) gmask[0][1][i] = (i+1)*(N+1)-1;
      for(int i = 0; i < N+1; ++i) gmask[0][2][i] = i;
      for(int i = 0; i < N+1; ++i) gmask[0][3][i] = (N+1)*N + i;

      /* corners */
      for(int j = 0; j < 2; ++j)
        for(int i = 0; i < 2; ++i)
          gmask[1][i+j*2][0] = i*N + j*(N+1);
      break;

    case 3:
      /* This could all probably be cleaned up... */

      /* faces */
      {
        int n,i,j,k,f=-1;

        n = 0; i = 0; f++;
        for(k = 0; k < N+1; k++) for(j = 0; j < N+1; j++)
          gmask[0][f][n++] = i+j*(N+1)+k*(N+1)*(N+1);

        n = 0; i = N; f++;
        for(k = 0; k < N+1; k++) for(j = 0; j < N+1; j++)
          gmask[0][f][n++] = i+j*(N+1)+k*(N+1)*(N+1);

        n = 0; j = 0; f++;
        for(k = 0; k < N+1; k++) for(i = 0; i < N+1; i++)
          gmask[0][f][n++] = i+j*(N+1)+k*(N+1)*(N+1);

        n = 0; j = N; f++;
        for(k = 0; k < N+1; k++) for(i = 0; i < N+1; i++)
          gmask[0][f][n++] = i+j*(N+1)+k*(N+1)*(N+1);

        n = 0; k = 0; f++;
        for(j = 0; j < N+1; j++) for(i = 0; i < N+1; i++)
          gmask[0][f][n++] = i+j*(N+1)+k*(N+1)*(N+1);

        n = 0; k = N; f++;
        for(j = 0; j < N+1; j++) for(i = 0; i < N+1; i++)
          gmask[0][f][n++] = i+j*(N+1)+k*(N+1)*(N+1);
      }

      /* edges */
      {
        int n,i,j,k,e = 0;

        for(k = 0; k < N+1;k+=N)
          for(j = 0; j < N+1;j+=N)
          {
            n = 0;
            for(i = 0; i < N+1; i++)
              gmask[1][e][n++] = i+j*(N+1)+k*(N+1)*(N+1);
            e++;
          }
        for(k = 0; k < N+1;k+=N)
          for(i = 0; i < N+1;i+=N)
          {
            n = 0;
            for(j = 0; j < N+1; j++)
              gmask[1][e][n++] = i+j*(N+1)+k*(N+1)*(N+1);
            e++;
          }
        for(j = 0; j < N+1;j+=N)
          for(i = 0; i < N+1;i+=N)
          {
            n = 0;
            for(k = 0; k < N+1; k++)
              gmask[1][e][n++] = i+j*(N+1)+k*(N+1)*(N+1);
            e++;
          }
      }

      /* corners */
      for(int k = 0, c = 0; k < N+1; k+=N)
        for(int j = 0; j < N+1;j+=N)
          for(int i = 0; i < N+1;i+=N)
            gmask[2][c++][0] = i+j*(N+1)+k*(N+1)*(N+1);

      break;

    default:
      BFAM_ABORT("cannot handle dim = %d",DIM);
  }

  return gmask;
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

  subdomain->gmask =
    BFAM_APPEND_EXPAND(bfam_subdomain_dgx_gmask_set_,BFAM_DGX_DIMENSION)(
                       numg, N, &subdomain->Np, Ng, Ngp, DIM);

  const int Np = subdomain->Np;

  if(DIM > 0)
  {
    const int Nrp = N+1;
    bfam_long_real_t *lr, *lw;
    lr = bfam_malloc_aligned(Nrp*sizeof(bfam_long_real_t));
    lw = bfam_malloc_aligned(Nrp*sizeof(bfam_long_real_t));

    bfam_jacobi_gauss_lobatto_quadrature(0, 0, N, lr, lw);

    bfam_long_real_t **lxi =
      bfam_malloc_aligned(num_Vi*sizeof(bfam_long_real_t*));

    for(int i = 0;i < num_Vi;i++)
      lxi[i] = bfam_malloc_aligned(K*Np*sizeof(bfam_long_real_t));

    /* Loop over all the elements and set up the grid*/
    int Ncorners = Ng[numg-1];
    for(bfam_locidx_t k = 0; k < K; ++k)
    {
      const bfam_locidx_t *v = EToV+Ncorners*k;
      bfam_long_real_t w[Ncorners];

      if(DIM == 1)
        for(int n = 0; n < Nrp; ++n)
        {
          w[0] = 1-lr[n];
          w[1] = 1+lr[n];
          for(int i = 0; i < num_Vi;i++)
          {
            lxi[i][Np*k + n] = 0;
            for(int c = 0; c < Ncorners; c++)
              lxi[i][Np*k + n] += w[c]*Vi[i][v[c]];
            lxi[i][Np*k + n] /= Ncorners;
          }
        }
      else if(DIM == 2)
        for(int n = 0; n < Nrp; ++n)
          for(int m = 0; m < Nrp; ++m)
          {
            w[0] = (1-lr[m])*(1-lr[n]);
            w[1] = (1+lr[m])*(1-lr[n]);
            w[2] = (1-lr[m])*(1+lr[n]);
            w[3] = (1+lr[m])*(1+lr[n]);
            for(int i = 0; i < num_Vi;i++)
            {
              lxi[i][Np*k + n] = 0;
              for(int c = 0; c < Ncorners; c++)
                lxi[i][Np*k + n] += w[c]*Vi[i][v[c]];
              lxi[i][Np*k + n] /= Ncorners;
            }
        }
      else if(DIM == 3)
        for(int n = 0; n < Nrp; ++n)
          for(int m = 0; m < Nrp; ++m)
            for(int l = 0; l < Nrp; ++l)
            {
              w[0] = (1-lr[m])*(1-lr[n])*(1-lr[l]);
              w[1] = (1+lr[m])*(1-lr[n])*(1-lr[l]);
              w[2] = (1-lr[m])*(1+lr[n])*(1-lr[l]);
              w[3] = (1+lr[m])*(1+lr[n])*(1-lr[l]);
              w[4] = (1-lr[m])*(1-lr[n])*(1-lr[l]);
              w[5] = (1+lr[m])*(1-lr[n])*(1-lr[l]);
              w[6] = (1-lr[m])*(1+lr[n])*(1-lr[l]);
              w[7] = (1+lr[m])*(1+lr[n])*(1-lr[l]);
              for(int i = 0; i < num_Vi;i++)
              {
                lxi[i][Np*k + n] = 0;
                for(int c = 0; c < Ncorners; c++)
                  lxi[i][Np*k + n] += w[c]*Vi[i][v[c]];
                lxi[i][Np*k + n] /= Ncorners;
              }
            }
      else BFAM_ABORT("not setup of dim = %d",DIM);

    }

    subdomain->V = bfam_malloc_aligned(Nrp*Nrp*sizeof(bfam_long_real_t));
    bfam_long_real_t *restrict V = subdomain->V;

    bfam_jacobi_p_vandermonde(0, 0, N, Nrp, lr, V);

    bfam_long_real_t *restrict D =
      bfam_malloc_aligned(Nrp*Nrp*sizeof(bfam_long_real_t));

    bfam_jacobi_p_differentiation(0, 0, N, Nrp, lr, V, D);

    bfam_long_real_t *restrict M =
      bfam_malloc_aligned(Nrp*Nrp*sizeof(bfam_long_real_t));

    bfam_jacobi_p_mass(0, 0, N, V, M);


    bfam_long_real_t **lJrx =
      bfam_malloc_aligned(num_Vi*DIM*sizeof(bfam_long_real_t*));
    for(int n = 0; n < num_Vi*DIM; n++)
      lJrx[n] = bfam_malloc_aligned(K*Np*sizeof(bfam_long_real_t));

    bfam_long_real_t *lJ = bfam_malloc_aligned(K*Np*sizeof(bfam_long_real_t));

    /* Ng[0] = number of faces, Ngp[0] = number of face points */
    bfam_long_real_t **ln =
      bfam_malloc_aligned(num_Vi*Ng[0]*sizeof(bfam_long_real_t*));
    for(int n = 0; n < num_Vi*Ng[0]; n++)
      ln[n] = bfam_malloc_aligned(K*Ng[0]*Ngp[0]*sizeof(bfam_long_real_t));

    bfam_long_real_t *lsJ =
      bfam_malloc_aligned(K*Ng[0]*Ngp[0]*sizeof(bfam_long_real_t));

    /*
    bfam_subdomain_dgx_geo(N, K, subdomain->gmask, lxi, D, lJrx, lJ, ln, lsJ,
        DIM);
    */

    subdomain->r = bfam_malloc_aligned(Nrp*sizeof(bfam_real_t));
    subdomain->w = bfam_malloc_aligned(Nrp*sizeof(bfam_real_t));
    subdomain->wi = bfam_malloc_aligned(Nrp*sizeof(bfam_real_t));

    for(int n = 0; n<Nrp; ++n)
    {
      subdomain->r[n]  = (bfam_real_t) lr[n];
      subdomain->w[n]  = (bfam_real_t) lw[n];
      subdomain->wi[n] = (bfam_real_t) (1.0l/lw[n]);
    }

    subdomain->K = K;

    /* free stuff */
    bfam_free_aligned(lsJ);
    for(int n = 0; n < num_Vi*Ng[0]; n++)
      bfam_free_aligned(ln[n]);
    bfam_free_aligned(ln);

    bfam_free_aligned(lJ);

    for(int n = 0; n < num_Vi*DIM; n++)
      bfam_free_aligned(lJrx[n]);
    bfam_free_aligned(lJrx);

    bfam_free_aligned(D);
    bfam_free_aligned(M);

    bfam_free_aligned(lr);
    bfam_free_aligned(lw);

    for(int i = 0;i < num_Vi;i++) bfam_free_aligned(lxi[i]);
    bfam_free_aligned(lxi);
  }
  else
  {
    subdomain->V  = NULL;
    subdomain->r  = NULL;
    subdomain->w  = NULL;
    subdomain->wi = NULL;
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
  if(sub->V) bfam_free_aligned(sub->V); sub->V = NULL;

  if(sub->r ) bfam_free_aligned(sub->r ); sub->r  = NULL;
  if(sub->w ) bfam_free_aligned(sub->w ); sub->w  = NULL;
  if(sub->wi) bfam_free_aligned(sub->wi); sub->wi = NULL;

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
