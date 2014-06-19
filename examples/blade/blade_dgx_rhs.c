#ifdef BLADE_DGX_DIMENSION

/* START DIM 2 */
#if    BLADE_DGX_DIMENSION==2
#include "blade_dgx_rhs_2.h"

#define DIM 2
#define BLADE_D3_AP(A1,A2) (A1)
#define BLADE_D3_OP(A) BFAM_NOOP()


#define BLADE_DR1(A,B)     BFAM_KRON_IXA    (N+1, Dr , A, B)
#define BLADE_DR1_PE(A,B)  BFAM_KRON_IXA_PE (N+1, Dr , A, B)
#define BLADE_DR1T(A,B)    BFAM_KRON_IXAT   (N+1, Dr , A, B)
#define BLADE_DR1T_PE(A,B) BFAM_KRON_IXAT_PE(N+1, Dr , A, B)

#define BLADE_DR2(A,B)     BFAM_KRON_AXI    (N+1, Dr , A, B)
#define BLADE_DR2_PE(A,B)  BFAM_KRON_AXI_PE (N+1, Dr , A, B)
#define BLADE_DR2T(A,B)    BFAM_KRON_ATXI   (N+1, Dr , A, B)
#define BLADE_DR2T_PE(A,B) BFAM_KRON_ATXI_PE(N+1, Dr , A, B)

#define BLADE_DR3(A,B)     BFAM_NOOP()
#define BLADE_DR3_PE(A,B)  BFAM_NOOP()
#define BLADE_DR3T(A,B)    BFAM_NOOP()
#define BLADE_DR3T_PE(A,B) BFAM_NOOP()

/* START DIM 2 Generic */
#ifndef NORDER
#define NORDER
#define USE_GENERIC
#define GENERIC_INIT(INN,NFM)\
  const int N   = INN; \
  const int Np  = (INN+1)*(INN+1); \
  const int Nfp = (INN+1); \
  const int Nfaces  = 4; \
  BFAM_WARNING(\
      "using generic NFM for 2d of order %d with "\
      "Np = %d, Nfp = %d, and Nfaces = %d",\
      N, Np, Nfp, Nfaces);
#else
#define GENERIC_INIT(INN,NFM) BFAM_NOOP()
#define N    NORDER
#define Np  (NORDER+1)*(NORDER+1)
#define Np_BACK  (NORDER+1)*(NORDER+1)
#define Nfp (NORDER+1)
#define Nfaces  4
#endif
/* END DIM 2 Generic */

/* END DIM 2 */
/* START DIM 3 */
#elif  BLADE_DGX_DIMENSION==3
#include "blade_dgx_rhs_3.h"

#define DIM 3
#define BLADE_D3_AP(A1,A2) (A1 A2)
#define BLADE_D3_OP(A) A

#define BLADE_DR1(A,B)     BFAM_KRON_IXIXA    (N+1, Dr , A, B)
#define BLADE_DR1_PE(A,B)  BFAM_KRON_IXIXA_PE (N+1, Dr , A, B)
#define BLADE_DR1T(A,B)    BFAM_KRON_IXIXAT   (N+1, Dr , A, B)
#define BLADE_DR1T_PE(A,B) BFAM_KRON_IXIXAT_PE(N+1, Dr , A, B)

#define BLADE_DR2(A,B)     BFAM_KRON_IXAXI    (N+1, Dr , A, B)
#define BLADE_DR2_PE(A,B)  BFAM_KRON_IXAXI_PE (N+1, Dr , A, B)
#define BLADE_DR2T(A,B)    BFAM_KRON_IXATXI   (N+1, Dr , A, B)
#define BLADE_DR2T_PE(A,B) BFAM_KRON_IXATXI_PE(N+1, Dr , A, B)

#define BLADE_DR3(A,B)     BFAM_KRON_AXIXI   (N+1, Dr , A, B)
#define BLADE_DR3_PE(A,B)  BFAM_KRON_AXIXI_PE(N+1, Dr , A, B)
#define BLADE_DR3T(A,B)    BFAM_KRON_ATXIXI   (N+1, Dr , A, B)
#define BLADE_DR3T_PE(A,B) BFAM_KRON_ATXIXI_PE(N+1, Dr , A, B)

#ifndef NORDER
#define NORDER
#define USE_GENERIC
#define GENERIC_INIT(INN,NFM)\
  const int N   = INN; \
  const int Np  = (INN+1)*(INN+1)*(INN+1); \
  const int Nfp = (INN+1)*(INN+1); \
  const int Nfaces  = 6; \
  BFAM_WARNING(\
      "using generic NFM for 3d of order %d with "\
      "Np = %d, Nfp = %d, and Nfaces = %d",\
      N, Np, Nfp, Nfaces);
#else
#define GENERIC_INIT(INN,NFM) BFAM_NOOP()
#define N    NORDER
#define Np  (NORDER+1)*(NORDER+1)*(NORDER+1)
#define Np_BACK  (NORDER+1)*(NORDER+1)*(NORDER+1)
#define Nfp (NORDER+1)*(NORDER+1)
#define Nfaces  6
#endif

/* END DIM 3 */
/* START BAD DIMENSION */
#else
#error "bad dimension"
#endif
/* END BAD DIMENSION */

/* Some useful field macros */
#define BFAM_LOAD_FIELD_RESTRICT_ALIGNED(field,prefix,base,dictionary)         \
bfam_real_t *restrict field;                                                   \
{                                                                              \
  char bfam_load_field_name[BFAM_BUFSIZ];                                      \
  snprintf(bfam_load_field_name,BFAM_BUFSIZ,"%s%s",(prefix),(base));           \
  field = bfam_dictionary_get_value_ptr(dictionary, bfam_load_field_name);     \
  BFAM_ASSERT(field != NULL);                                                  \
}                                                                              \
BFAM_ASSUME_ALIGNED(field,32);

#define BFAM_LOAD_FIELD_ALIGNED(field,prefix,base,dictionary)                  \
bfam_real_t *field;                                                            \
{                                                                              \
  char bfam_load_field_name[BFAM_BUFSIZ];                                      \
  snprintf(bfam_load_field_name,BFAM_BUFSIZ,"%s%s",(prefix),(base));           \
  field = bfam_dictionary_get_value_ptr(dictionary, bfam_load_field_name);     \
  BFAM_ASSERT(field != NULL);                                                  \
}                                                                              \
BFAM_ASSUME_ALIGNED(field,32);


/* Function macros */
#define BLADE_APPEND_4(a,b,c,d) a ## b ## c ##d
#define BLADE_APPEND_EXPAND_4(a,b,c,d) BLADE_APPEND_4(a,b,c,d)
#define blade_dgx_energy \
  BLADE_APPEND_EXPAND_4(blade_dgx_energy_,DIM,_,NORDER)

void blade_dgx_energy(
    int inN, bfam_real_t *energy_sq,
    bfam_subdomain_dgx_t *sub, const char *field_prefix)
{
  GENERIC_INIT(inN,blade_dgx_energy);

  /* get the fields we will need */
  bfam_dictionary_t *fields = &sub->base.fields;
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(q ,field_prefix,"q" ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(J  ,"","_grid_J",fields);

  bfam_locidx_t K  = sub->K;

  bfam_real_t *w  = sub->w;
  BFAM_ASSUME_ALIGNED(w ,32);

  /* loop through all the elements */
  for(bfam_locidx_t e = 0; e < K;e++)
  {
    bfam_locidx_t n = e*Np;
#if DIM==3
    for(bfam_locidx_t k = 0; k < N+1; k++)
    {
      const bfam_real_t wk = w[k];
#endif
      for(bfam_locidx_t j = 0; j < N+1; j++)
      {
        const bfam_real_t wj = w[j];
        for(bfam_locidx_t i = 0; i < N+1; i++,n++)
        {
          const bfam_real_t wi = w[i];

          energy_sq[0] += BLADE_D3_AP(wi*wj,*wk)*J[n]*q[n]*q[n];
        }
      }
#if DIM==3
    }
#endif
  }
}

#endif
