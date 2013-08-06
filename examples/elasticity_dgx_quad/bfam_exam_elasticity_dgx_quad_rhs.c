#include "bfam_exam_elasticity_dgx_quad_rhs.h"

#ifndef NORDER
#define NORDER
#define USE_GENERIC
#else
#define N    NORDER
#define Np  (NORDER+1)*(NORDER+1)
#define Nfp (NORDER+1)
#endif

void BFAM_APPEND_EXPAND(bfam_elasticity_dgx_quad_print_order_,NORDER)(int inN)
{

#ifdef USE_GENERIC
  const int N = inN;
  BFAM_WARNING("Using generic order print function");
#endif

  BFAM_INFO("N = %d", N);
}

void BFAM_APPEND_EXPAND(bfam_elasticity_dgx_quad_intra_rhs_,NORDER)(int inN,
    bfam_subdomain_dgx_quad_t *sub, const char *rate_prefix,
    const char *field_prefix)
{

#ifdef USE_GENERIC
  const int N   = inN;
  const int Np  = (inN+1)*(inN+1);
  const int Nfp = inN+1;
  BFAM_WARNING("Using generic intra rhs function");
#endif

  char tmp_name[BFAM_BUFSIZ];

  /* get the fields we will need */
  bfam_dictionary_t *fields = &sub->base.fields;

  snprintf(tmp_name,BFAM_BUFSIZ,"%sv1",field_prefix);
  bfam_real_t *restrict v1 = bfam_dictionary_get_value_ptr(fields, tmp_name);

  snprintf(tmp_name,BFAM_BUFSIZ,"%sv2",field_prefix);
  bfam_real_t *restrict v2 = bfam_dictionary_get_value_ptr(fields, tmp_name);

  snprintf(tmp_name,BFAM_BUFSIZ,"%sv3",field_prefix);
  bfam_real_t *restrict v3 = bfam_dictionary_get_value_ptr(fields, tmp_name);

  snprintf(tmp_name,BFAM_BUFSIZ,"%sS11",field_prefix);
  bfam_real_t *restrict S11 = bfam_dictionary_get_value_ptr(fields, tmp_name);

  snprintf(tmp_name,BFAM_BUFSIZ,"%sS12",field_prefix);
  bfam_real_t *restrict S12 = bfam_dictionary_get_value_ptr(fields, tmp_name);

  snprintf(tmp_name,BFAM_BUFSIZ,"%sS13",field_prefix);
  bfam_real_t *restrict S13 = bfam_dictionary_get_value_ptr(fields, tmp_name);

  snprintf(tmp_name,BFAM_BUFSIZ,"%sS22",field_prefix);
  bfam_real_t *restrict S22 = bfam_dictionary_get_value_ptr(fields, tmp_name);

  snprintf(tmp_name,BFAM_BUFSIZ,"%sS23",field_prefix);
  bfam_real_t *restrict S23 = bfam_dictionary_get_value_ptr(fields, tmp_name);

  /*snprintf(tmp_name,BFAM_BUFSIZ,"%sS33",field_prefix);
  bfam_real_t *restrict S33 = bfam_dictionary_get_value_ptr(fields, tmp_name);*/

  /* get the rates we will need */
  snprintf(tmp_name,BFAM_BUFSIZ,"%sv1",rate_prefix);
  bfam_real_t *restrict dv1 = bfam_dictionary_get_value_ptr(fields, tmp_name);

  snprintf(tmp_name,BFAM_BUFSIZ,"%sv2",rate_prefix);
  bfam_real_t *restrict dv2 = bfam_dictionary_get_value_ptr(fields, tmp_name);

  snprintf(tmp_name,BFAM_BUFSIZ,"%sv3",rate_prefix);
  bfam_real_t *restrict dv3 = bfam_dictionary_get_value_ptr(fields, tmp_name);

  snprintf(tmp_name,BFAM_BUFSIZ,"%sS11",rate_prefix);
  bfam_real_t *restrict dS11 = bfam_dictionary_get_value_ptr(fields, tmp_name);

  snprintf(tmp_name,BFAM_BUFSIZ,"%sS12",rate_prefix);
  bfam_real_t *restrict dS12 = bfam_dictionary_get_value_ptr(fields, tmp_name);

  snprintf(tmp_name,BFAM_BUFSIZ,"%sS13",rate_prefix);
  bfam_real_t *restrict dS13 = bfam_dictionary_get_value_ptr(fields, tmp_name);

  snprintf(tmp_name,BFAM_BUFSIZ,"%sS22",rate_prefix);
  bfam_real_t *restrict dS22 = bfam_dictionary_get_value_ptr(fields, tmp_name);

  snprintf(tmp_name,BFAM_BUFSIZ,"%sS23",rate_prefix);
  bfam_real_t *restrict dS23 = bfam_dictionary_get_value_ptr(fields, tmp_name);

  snprintf(tmp_name,BFAM_BUFSIZ,"%sS33",rate_prefix);
  bfam_real_t *restrict dS33 = bfam_dictionary_get_value_ptr(fields, tmp_name);

  /* get the material properties and metric terms */
  /* bfam_real_t *restrict rho = bfam_dictionary_get_value_ptr(fields,"rho"); */
  bfam_real_t *restrict lam = bfam_dictionary_get_value_ptr(fields,"lam");
  bfam_real_t *restrict mu  = bfam_dictionary_get_value_ptr(fields,"mu" );
  bfam_real_t *restrict J   = bfam_dictionary_get_value_ptr(fields,"_grid_J");
  bfam_real_t *restrict Jrx = bfam_dictionary_get_value_ptr(fields,"_grid_Jrx");
  bfam_real_t *restrict Jry = bfam_dictionary_get_value_ptr(fields,"_grid_Jry");
  bfam_real_t *restrict Jsx = bfam_dictionary_get_value_ptr(fields,"_grid_Jsx");
  bfam_real_t *restrict Jsy = bfam_dictionary_get_value_ptr(fields,"_grid_Jsy");

  bfam_locidx_t K  = sub->K;
  BFAM_ALIGN(32) bfam_real_t aux1[Np];
  BFAM_ALIGN(32) bfam_real_t aux2[Np];
  BFAM_ALIGN(32) bfam_real_t rate[Np];

  BFAM_ALIGN(32) bfam_real_t MJv1[Np];
  BFAM_ALIGN(32) bfam_real_t MJv2[Np];
  BFAM_ALIGN(32) bfam_real_t MJv3[Np];

  bfam_real_t *Dr = sub->Dr;
  bfam_real_t *w  = sub->w;

  bfam_locidx_t num_fp = 4*Nfp;
  /* bfam_locidx_t *vmapM = sub->vmapM; */
  /* bfam_locidx_t *vmapP = sub->vmapP; */


  /* loop through all the elements */
  for(bfam_locidx_t e = 0; e < K;e++)
  {
    bfam_locidx_t off = e*Np;

    /* do differential equation */
    /* Note: Matrices on the LHS will be handled in add rates routine */

    /* M*J*rho*v1 += M*(Dr*(Jrx*S11+Jry*S12)) +  M*(Ds*(Jsx*S11+Jsy*S12));*/
    BFAM_DOT_AX    (Np , Jrx+off, S11+off, aux1);
    BFAM_DOT_AX_PE (Np , Jry+off, S12+off, aux1);
    BFAM_KRON_AXI   (N+1, Dr   , aux1  , rate);

    BFAM_DOT_AX    (Np , Jsx+off, S11+off, aux1);
    BFAM_DOT_AX_PE (Np , Jsy+off, S12+off, aux1);
    BFAM_KRON_IXA_PE(N+1, Dr   , aux1  , rate);

    BFAM_KRON_AB_DOT_C_PE(N+1,w,w,J+off,rate,dv1+off);

    /* M*J*rho*v2 += M*(Dr*(Jrx*S12+Jry*S22)) +  M*(Ds*(Jsx*S12+Jsy*S22));*/
    BFAM_DOT_AX    (Np , Jrx+off, S12+off, aux1);
    BFAM_DOT_AX_PE (Np , Jry+off, S22+off, aux1);
    BFAM_KRON_AXI   (N+1, Dr   , aux1  , rate);

    BFAM_DOT_AX    (Np , Jsx+off, S12+off, aux1);
    BFAM_DOT_AX_PE (Np , Jsy+off, S22+off, aux1);
    BFAM_KRON_IXA_PE(N+1, Dr   , aux1  , rate);

    BFAM_KRON_AB_DOT_C_PE(N+1,w,w,J+off,rate,dv2+off);

    /* M*J*rho*v3 += M*(Dr*(Jrx*S13+Jry*S23)) +  M*(Ds*(Jsx*S13+Jsy*S23));*/
    BFAM_DOT_AX    (Np , Jrx+off, S13+off, aux1);
    BFAM_DOT_AX_PE (Np , Jry+off, S23+off, aux1);
    BFAM_KRON_AXI   (N+1, Dr   , aux1  , rate);

    BFAM_DOT_AX    (Np , Jsx+off, S13+off, aux1);
    BFAM_DOT_AX_PE (Np , Jsy+off, S23+off, aux1);
    BFAM_KRON_IXA_PE(N+1, Dr   , aux1  , rate);

    BFAM_KRON_AB_DOT_C_PE(N+1,w,w,J+off,rate,dv3+off);

    /* we use these a lot so store them */
    BFAM_KRON_AB_DOT_C(N+1,w,w,J+off,v1+off,MJv1);
    BFAM_KRON_AB_DOT_C(N+1,w,w,J+off,v2+off,MJv2);
    BFAM_KRON_AB_DOT_C(N+1,w,w,J+off,v3+off,MJv3);

    /* M*J*S11 += -(lam+mu)*(Jrx*Dr'+Jsx*Ds')*v1 - lam*(Jry*Dr'+Jsy*Ds')*v2 */
    /* M*J*S22 += -lam*(Jrx*Dr'+Jsx*Ds')*v1 - (lam+mu)*(Jry*Dr'+Jsy*Ds')*v2 */
    /* M*J*S33 += -lam*(Jrx*Dr'+Jsx*Ds')*v1 - lam*(Jry*Dr'+Jsy*Ds')*v2 */

    BFAM_KRON_ATXI(N+1, Dr     , MJv1, aux1); /* a1=Dr'*v1 */
    BFAM_DOT_AX   (Np , Jrx+off, aux1, aux2); /* a2=Jrx*Dr'*v1 */
    BFAM_KRON_IXAT(N+1, Dr     , MJv1, aux1); /* a1=Ds'v1 */
    BFAM_DOT_AX_PE(Np , Jsx+off, aux1, aux2); /* a2=Jrx*Dr'*v1+Jsx*Ds'*v1 */
    /* once computed aux2 can be used for all diagonal components of S */
    BFAM_DOT_AX_ME(Np , mu +off, aux2, dS11);
    BFAM_DOT_AX_ME(Np , lam+off, aux2, dS11);
    BFAM_DOT_AX_ME(Np , lam+off, aux2, dS22);
    BFAM_DOT_AX_ME(Np , lam+off, aux2, dS33);

    BFAM_KRON_ATXI(N+1, Dr     , MJv2, aux1); /* a1=Dr'*v2 */
    BFAM_DOT_AX   (Np , Jry+off, aux1, aux2); /* a2= (..)*v1 + Jry*Dr'*v2 */
    BFAM_KRON_IXAT(N+1, Dr     , MJv2, aux1); /* a1=Ds'*v2 */
    BFAM_DOT_AX_PE(Np , Jsy+off, aux1, aux2); /* a2= (..)*v1 + (..)v2 */
    /* once computed aux2 can be used for all diagonal components of S */
    BFAM_DOT_AX_ME(Np , mu +off, aux2, dS22);
    BFAM_DOT_AX_ME(Np , lam+off, aux2, dS11);
    BFAM_DOT_AX_ME(Np , lam+off, aux2, dS22);
    BFAM_DOT_AX_ME(Np , lam+off, aux2, dS33);

    /* M*J*S12 += -mu*((Jrx*Dr'+Jsx*Ds')*v2 -mu*((Jry*Dr'+Jsy*Ds')*v1 */
    BFAM_KRON_ATXI(N+1, Dr     , MJv2, aux1); /* a1=Dr'*v2 */
    BFAM_DOT_AX   (Np , Jrx+off, aux2, aux2); /* a2=Jrx*Dr'*v2 */
    BFAM_KRON_IXAT(N+1, Dr     , MJv2, aux1); /* a1=Ds'v2 */
    BFAM_DOT_AX_PE(Np , Jsx+off, aux2, aux2); /* a2=Jrx*Dr'*v2+Jsx*Ds'*v2 */
    BFAM_DOT_AX_ME(Np , mu +off, aux2, dS12);

    BFAM_KRON_ATXI(N+1, Dr     , MJv1, aux1); /* a1=Dr'*v1 */
    BFAM_DOT_AX   (Np , Jry+off, aux2, aux2); /* a2=Jry*Dr'*v1 */
    BFAM_KRON_IXAT(N+1, Dr     , MJv1, aux1); /* a1=Ds'v1 */
    BFAM_DOT_AX_PE(Np , Jsy+off, aux2, aux2); /* a2=Jry*Dr'*v1+Jsy*Ds'*v1 */
    BFAM_DOT_AX_ME(Np , mu +off, aux2, dS12);

    /* M*J*S13 += -mu*((Jrx*Dr'+Jsx*Ds')*v3 */
    BFAM_KRON_ATXI(N+1, Dr     , MJv3, aux1); /* a1=Dr'*v3 */
    BFAM_DOT_AX   (Np , Jrx+off, aux2, aux2); /* a2=Jrx*Dr'*v3 */
    BFAM_KRON_IXAT(N+1, Dr     , MJv3, aux1); /* a1=Ds'v3 */
    BFAM_DOT_AX_PE(Np , Jsx+off, aux2, aux2); /* a2=Jrx*Dr'*v3+Jsx*Ds'*v3 */
    BFAM_DOT_AX_ME(Np , mu +off, aux2, dS13);

    /* M*J*S23 += -mu*((Jry*Dr'+Jsy*Ds')*v3 */
    BFAM_KRON_ATXI(N+1, Dr     , MJv3, aux1); /* a1=Dr'*v3 */
    BFAM_DOT_AX   (Np , Jry+off, aux2, aux2); /* a2=Jry*Dr'*v3 */
    BFAM_KRON_IXAT(N+1, Dr     , MJv3, aux1); /* a1=Ds'v3 */
    BFAM_DOT_AX_PE(Np , Jsy+off, aux2, aux2); /* a2=Jry*Dr'*v3+Jsy*Ds'*v3 */
    BFAM_DOT_AX_ME(Np , mu +off, aux2, dS23);


    /* loop over faces */
    for(bfam_locidx_t f = e*num_fp; f < (e+1)*num_fp; f++)
    {
    }
  }
}
