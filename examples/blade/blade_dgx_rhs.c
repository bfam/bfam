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

#define blade_dgx_scale_rates_advection \
  BLADE_APPEND_EXPAND_4(blade_dgx_scale_rates_advection_,DIM,_,NORDER)

#define blade_dgx_intra_rhs_advection \
  BLADE_APPEND_EXPAND_4(blade_dgx_intra_rhs_advection_,DIM,_,NORDER)

#define blade_dgx_energy \
  BLADE_APPEND_EXPAND_4(blade_dgx_energy_,DIM,_,NORDER)

void blade_dgx_scale_rates_advection(
    int inN, bfam_subdomain_dgx_t *sub, const char *rate_prefix,
    const bfam_long_real_t a)
{
  GENERIC_INIT(inN,blade_dgx_scale_rates_advection);

  const bfam_locidx_t num_pts = sub->K * Np;
  bfam_dictionary_t *fields = &sub->base.fields;

  const char *f_names[] = {"q",NULL};

  for(bfam_locidx_t f = 0; f_names[f]!=NULL; ++f)
  {
    BFAM_LOAD_FIELD_RESTRICT_ALIGNED(rate,rate_prefix,f_names[f],fields);
    for(bfam_locidx_t n = 0; n < num_pts; n++) rate[n] *= a;
  }
}

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

void blade_dgx_intra_rhs_advection(
    int inN, bfam_subdomain_dgx_t *sub, const char *rate_prefix,
    const char *field_prefix, const bfam_long_real_t t)
{
  GENERIC_INIT(inN,blade_dgx_intra_rhs_advection);

  /* get the fields we will need */
  bfam_dictionary_t *fields = &sub->base.fields;
  bfam_dictionary_t *fields_face = &sub->base.fields_face;

  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(q  ,field_prefix,"q" ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(dq ,rate_prefix ,"q" ,fields);

  /* get the material properties and metric terms */
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(ux ,field_prefix,"ux" ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(uy ,field_prefix,"uy" ,fields);
  BLADE_D3_OP(BFAM_LOAD_FIELD_RESTRICT_ALIGNED(uz,field_prefix,"uz",fields));

  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(J    ,"","_grid_J"    ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(JI   ,"","_grid_JI"   ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Jr1x1,"","_grid_Jr0x0",fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Jr1x2,"","_grid_Jr0x1",fields);
  BLADE_D3_OP(BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Jr1x3,"","_grid_Jr0x2",fields));
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Jr2x1,"","_grid_Jr1x0",fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Jr2x2,"","_grid_Jr1x1",fields);
  BLADE_D3_OP(BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Jr2x3,"","_grid_Jr1x2",fields));

  BLADE_D3_OP(BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Jr3x1,"","_grid_Jr2x0",fields));
  BLADE_D3_OP(BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Jr3x2,"","_grid_Jr2x1",fields));
  BLADE_D3_OP(BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Jr3x3,"","_grid_Jr2x2",fields));

  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(n1,"","_grid_nx0",fields_face);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(n2,"","_grid_nx1",fields_face);
  BLADE_D3_OP(BFAM_LOAD_FIELD_RESTRICT_ALIGNED(n3,"","_grid_nx2",fields_face));
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(sJ  ,"","_grid_sJ",fields_face);

  bfam_locidx_t K  = sub->K;

  BFAM_ALIGN(32) bfam_real_t aux0[Np];
  BFAM_ALIGN(32) bfam_real_t aux1[Np];
  BFAM_ALIGN(32) bfam_real_t aux2[Np];
  BFAM_ALIGN(32) bfam_real_t aux3[Np];
//JK  BFAM_ALIGN(32) bfam_real_t aux2[Np];
//JK  BLADE_D3_OP(BFAM_ALIGN(32) bfam_real_t aux3[Np]);
//JK
//JK  BFAM_ALIGN(32) bfam_real_t Dr1Tv1[Np];
//JK  BFAM_ALIGN(32) bfam_real_t Dr2Tv1[Np];
//JK  BLADE_D3_OP(BFAM_ALIGN(32) bfam_real_t Dr3Tv1[Np]);
//JK  BFAM_ALIGN(32) bfam_real_t Dr1Tv2[Np];
//JK  BFAM_ALIGN(32) bfam_real_t Dr2Tv2[Np];
//JK  BLADE_D3_OP(BFAM_ALIGN(32) bfam_real_t Dr3Tv2[Np]);
//JK  BFAM_ALIGN(32) bfam_real_t Dr1Tv3[Np];
//JK  BFAM_ALIGN(32) bfam_real_t Dr2Tv3[Np];
//JK  BLADE_D3_OP(BFAM_ALIGN(32) bfam_real_t Dr3Tv3[Np]);
//JK
//JK  BFAM_ALIGN(32) bfam_real_t Mv1[Np];
//JK  BFAM_ALIGN(32) bfam_real_t Mv2[Np];
//JK  BFAM_ALIGN(32) bfam_real_t Mv3[Np];
//JK
  bfam_real_t *Dr = sub->Dr;
  bfam_real_t *w  = sub->w;
  bfam_real_t *wi  = sub->wi;
  BFAM_ASSUME_ALIGNED(Dr,32);
  BFAM_ASSUME_ALIGNED(w ,32);
  BFAM_ASSUME_ALIGNED(wi ,32);

  bfam_locidx_t *vmapM = sub->vmapM;
  bfam_locidx_t *vmapP = sub->vmapP;


  /* loop through all the elements */
  for(bfam_locidx_t e = 0; e < K;e++)
  {
    bfam_locidx_t off = e*Np;

    /* do differential equation */
    /* Note: Matrices on the LHS will be handled in add rates routine */

    /****** X-DERIVATIVE ******/
    /* q += JI*(Jrx*Dr*(ux*q) + Jsx*Ds*(ux*q) + Jtx*Dt*(ux*q)) */
    for(bfam_locidx_t i = 0; i < Np; i++) aux0[i] = ux[off+i] * q[off+i];
    BLADE_DR1(aux0, aux1); /* Dr */
    BLADE_DR2(aux0, aux2); /* Ds */
    BLADE_DR3(aux0, aux3); /* Dt */
    for(bfam_locidx_t i = 0; i < Np; i++)
          dq[off+i] += BFAM_REAL(0.5)*JI[off+i]
                        *BLADE_D3_AP( Jr1x1[off+i]*aux1[i]
                                    + Jr2x1[off+i]*aux2[i],
                                    + Jr3x1[off+i]*aux3[i]);

    /* q += MI*JI*ux*(Dr'*M*Jrx*q + Ds'*M*Jsx*q + Dt'*M*Jtx*q) */
    bfam_locidx_t n = 0;
#if DIM==3
    for(bfam_locidx_t k = 0; k < N+1; k++)
#endif
      for(bfam_locidx_t j = 0; j < N+1; j++)
        for(bfam_locidx_t i = 0; i < N+1; i++,n++)
        {
          aux1[n] = BLADE_D3_AP(w[i]*w[j],*w[k])*Jr1x1[off+i]*q[off+i];
          aux2[n] = BLADE_D3_AP(w[i]*w[j],*w[k])*Jr2x1[off+i]*q[off+i];
          BLADE_D3_OP(aux3[n] = w[i]*w[j] *w[k] *Jr3x1[off+i]*q[off+i]);
        }
    BLADE_DR1T   (aux1, aux0);
    BLADE_DR2T_PE(aux2, aux0);
    BLADE_DR3T_PE(aux3, aux0);

    n = 0;
#if DIM==3
    for(bfam_locidx_t k = 0; k < N+1; k++)
#endif
      for(bfam_locidx_t j = 0; j < N+1; j++)
        for(bfam_locidx_t i = 0; i < N+1; i++,n++)
        {
          const bfam_real_t wi_ijk = BLADE_D3_AP(wi[i]*wi[j],*wi[k]);
          dq[off+n] -= BFAM_REAL(0.5)*wi_ijk*JI[off+n]*ux[off+n]*aux0[n];
        }


    /****** Y-DERIVATIVE ******/
    /* q += JI*(Jry*Dr*(uy*q) + Jsy*Ds*(uy*q) + Jty*Dt*(uy*q)) */
    for(bfam_locidx_t i = 0; i < Np; i++) aux0[i] = uy[off+i] * q[off+i];
    BLADE_DR1(aux0, aux1); /* Dr */
    BLADE_DR2(aux0, aux2); /* Ds */
    BLADE_DR3(aux0, aux3); /* Dt */
    for(bfam_locidx_t i = 0; i < Np; i++)
          dq[off+i] += BFAM_REAL(0.5)*JI[off+i]
                        *BLADE_D3_AP( Jr1x2[off+i]*aux1[i]
                                    + Jr2x2[off+i]*aux2[i],
                                    + Jr3x2[off+i]*aux3[i]);

    /* q += MI*JI*uy*(Dr'*M*Jry*q + Ds'*M*Jsy*q + Dt'*M*Jty*q) */
    n = 0;
#if DIM==3
    for(bfam_locidx_t k = 0; k < N+1; k++)
#endif
      for(bfam_locidx_t j = 0; j < N+1; j++)
        for(bfam_locidx_t i = 0; i < N+1; i++,n++)
        {
          aux1[n] = BLADE_D3_AP(w[i]*w[j],*w[k])*Jr1x2[off+i]*q[off+i];
          aux2[n] = BLADE_D3_AP(w[i]*w[j],*w[k])*Jr2x2[off+i]*q[off+i];
          BLADE_D3_OP(aux3[n] = w[i]*w[j] *w[k] *Jr3x2[off+i]*q[off+i]);
        }
    BLADE_DR1T   (aux1, aux0);
    BLADE_DR2T_PE(aux2, aux0);
    BLADE_DR3T_PE(aux3, aux0);

    n = 0;
#if DIM==3
    for(bfam_locidx_t k = 0; k < N+1; k++)
#endif
      for(bfam_locidx_t j = 0; j < N+1; j++)
        for(bfam_locidx_t i = 0; i < N+1; i++,n++)
        {
          const bfam_real_t wi_ijk = BLADE_D3_AP(wi[i]*wi[j],*wi[k]);
          dq[off+n] -= BFAM_REAL(0.5)*wi_ijk*JI[off+n]*uy[off+n]*aux0[n];
        }

    /****** Z-DERIVATIVE ******/
#if DIM==3
    /* q += JI*(Jrz*Dr*(uz*q) + Jsz*Ds*(uz*q) + Jtz*Dt*(uz*q)) */
    for(bfam_locidx_t i = 0; i < Np; i++) aux0[i] = uz[off+i] * q[off+i];
    BLADE_DR1(aux0, aux1); /* Dr */
    BLADE_DR2(aux0, aux2); /* Ds */
    BLADE_DR3(aux0, aux3); /* Dt */
    for(bfam_locidx_t i = 0; i < Np; i++)
          dq[off+i] += BFAM_REAL(0.5)*JI[off+i]
                        *BLADE_D3_AP( Jr1x3[off+i]*aux1[i]
                                    + Jr2x3[off+i]*aux2[i],
                                    + Jr3x3[off+i]*aux3[i]);

    /* q += MI*JI*uz*(Dr'*M*Jrz*q + Ds'*M*Jsz*q + Dt'*M*Jtz*q) */
    n = 0;
    for(bfam_locidx_t k = 0; k < N+1; k++)
      for(bfam_locidx_t j = 0; j < N+1; j++)
        for(bfam_locidx_t i = 0; i < N+1; i++,n++)
        {
          aux1[n] = BLADE_D3_AP(w[i]*w[j],*w[k])*Jr1x3[off+i]*q[off+i];
          aux2[n] = BLADE_D3_AP(w[i]*w[j],*w[k])*Jr2x3[off+i]*q[off+i];
          BLADE_D3_OP(aux3[n] = w[i]*w[j] *w[k] *Jr3x3[off+i]*q[off+i]);
        }
    BLADE_DR1T   (aux1, aux0);
    BLADE_DR2T_PE(aux2, aux0);
    BLADE_DR3T_PE(aux3, aux0);

    n = 0;
    for(bfam_locidx_t k = 0; k < N+1; k++)
      for(bfam_locidx_t j = 0; j < N+1; j++)
        for(bfam_locidx_t i = 0; i < N+1; i++,n++)
        {
          const bfam_real_t wi_ijk = BLADE_D3_AP(wi[i]*wi[j],*wi[k]);
          dq[off+n] -= BFAM_REAL(0.5)*wi_ijk*JI[off+n]*uz[off+n]*aux0[n];
        }
#endif

//JK    /* loop over faces */
//JK    for(bfam_locidx_t face = 0; face < Nfaces; face++)
//JK    {
//JK      for(bfam_locidx_t pnt = 0; pnt < Nfp; pnt++)
//JK      {
//JK        const bfam_locidx_t f = pnt + Nfp*(face + Nfaces*e);
//JK        const bfam_locidx_t iM = vmapM[f];
//JK        const bfam_locidx_t iP = vmapP[f];
//JK
//JK        /* Setup stuff for the minus side */
//JK        const bfam_real_t ZsM = Zs[iM];
//JK        const bfam_real_t ZpM = Zp[iM];
//JK
//JK        const bfam_real_t nM[] = {n1[f],n2[f],BLADE_D3_AP(0,+n3[f])};
//JK
//JK        bfam_real_t TpM[] = {
//JK          BLADE_D3_AP(nM[0]*S11[iM] + nM[1]*S12[iM], + nM[2]*S13[iM]),
//JK          BLADE_D3_AP(nM[0]*S12[iM] + nM[1]*S22[iM], + nM[2]*S23[iM]),
//JK          BLADE_D3_AP(nM[0]*S13[iM] + nM[1]*S23[iM], + nM[2]*S33[iM]),
//JK        };
//JK        const bfam_real_t TnM = BLADE_D3_AP(TpM[0]*nM[0]
//JK                                           +TpM[1]*nM[1],
//JK                                           +TpM[2]*nM[2]);
//JK        TpM[0] = TpM[0]-TnM*nM[0];
//JK        TpM[1] = TpM[1]-TnM*nM[1];
//JK        BLADE_D3_OP(TpM[2] = TpM[2]-TnM*nM[2]);
//JK
//JK        bfam_real_t vpM[] = {v1[iM],v2[iM],v3[iM]};
//JK        const bfam_real_t vnM = BLADE_D3_AP(nM[0]*vpM[0]
//JK                                           +nM[1]*vpM[1],
//JK                                           +nM[2]*vpM[2]);
//JK        vpM[0] = vpM[0]-vnM*nM[0];
//JK        vpM[1] = vpM[1]-vnM*nM[1];
//JK        BLADE_D3_OP(vpM[2] = vpM[2]-vnM*nM[2]);
//JK
//JK        /* Setup stuff for the plus side */
//JK        bfam_real_t ZsP = Zs[iP];
//JK        bfam_real_t ZpP = Zp[iP];
//JK
//JK        bfam_real_t nP[] = {-nM[0],-nM[1],-nM[2]};
//JK
//JK        bfam_real_t TpP[] = {
//JK          BLADE_D3_AP(nP[0]*S11[iP] + nP[1]*S12[iP], + nP[2]*S13[iP]),
//JK          BLADE_D3_AP(nP[0]*S12[iP] + nP[1]*S22[iP], + nP[2]*S23[iP]),
//JK          BLADE_D3_AP(nP[0]*S13[iP] + nP[1]*S23[iP], + nP[2]*S33[iP]),
//JK        };
//JK        const bfam_real_t TnP = BLADE_D3_AP(TpP[0]*nP[0]
//JK                                           +TpP[1]*nP[1],
//JK                                           +TpP[2]*nP[2]);
//JK        TpP[0] = TpP[0]-TnP*nP[0];
//JK        TpP[1] = TpP[1]-TnP*nP[1];
//JK        BLADE_D3_OP(TpP[2] = TpP[2]-TnP*nP[2]);
//JK
//JK        bfam_real_t vpP[] = {v1[iP],v2[iP],v3[iP]};
//JK        const bfam_real_t vnP = BLADE_D3_AP(nP[0]*vpP[0]
//JK                                           +nP[1]*vpP[1],
//JK                                           +nP[2]*vpP[2]);
//JK        vpP[0] = vpP[0]-vnP*nP[0];
//JK        vpP[1] = vpP[1]-vnP*nP[1];
//JK        BLADE_D3_OP(vpP[2] = vpP[2]-vnP*nP[2]);
//JK
//JK        bfam_real_t TnS;
//JK        bfam_real_t TpS[3];
//JK        bfam_real_t vnS;
//JK        bfam_real_t vpS[3];
//JK
//JK        BLADE_STATE(&TnS,TpS,&vnS,vpS,
//JK            TnM, TnP, TpM, TpP, vnM, vnP, vpM, vpP, ZpM, ZpP, ZsM, ZsP);
//JK
//JK        TnS    -= TnM;
//JK        TpS[0] -= TpM[0];
//JK        TpS[1] -= TpM[1];
//JK        TpS[2] -= TpM[2];
//JK
//JK        /* intra */
//JK        blade_dgx_add_flux(1, TnS,TpS,vnS,vpS,iM,
//JK            dv1,dv2,dv3, dS11,dS22,dS33,dS12,dS13,dS23,
//JK            lam[iM],mu[iM],rhoi[iM],nM,sJ[f],JI[iM],wi[0]);
//JK      }
//JK    }
  }
}

#endif
