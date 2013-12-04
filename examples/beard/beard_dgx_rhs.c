#ifdef BEARD_DGX_DIMENSION

#if    BEARD_DGX_DIMENSION==2
#include "beard_dgx_rhs_2.h"

#define DIM 2
#define BEARD_D3_AP(A1,A2) (A1)
#define BEARD_D3_OP(A) BFAM_NOOP()


#define BEARD_DR1(A,B)     BFAM_KRON_IXA    (N+1, Dr , A, B)
#define BEARD_DR1_PE(A,B)  BFAM_KRON_IXA_PE (N+1, Dr , A, B)
#define BEARD_DR1T(A,B)    BFAM_KRON_IXAT   (N+1, Dr , A, B)
#define BEARD_DR1T_PE(A,B) BFAM_KRON_IXAT_PE(N+1, Dr , A, B)

#define BEARD_DR2(A,B)     BFAM_KRON_AXI    (N+1, Dr , A, B)
#define BEARD_DR2_PE(A,B)  BFAM_KRON_AXI_PE (N+1, Dr , A, B)
#define BEARD_DR2T(A,B)    BFAM_KRON_ATXI   (N+1, Dr , A, B)
#define BEARD_DR2T_PE(A,B) BFAM_KRON_ATXI_PE(N+1, Dr , A, B)

#define BEARD_DR3(A,B)     BFAM_NOOP()
#define BEARD_DR3_PE(A,B)  BFAM_NOOP()
#define BEARD_DR3T(A,B)    BFAM_NOOP()
#define BEARD_DR3T_PE(A,B) BFAM_NOOP()

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
#define Nfp (NORDER+1)
#define Nfaces  4
#endif

#elif  BEARD_DGX_DIMENSION==3
#include "beard_dgx_rhs_3.h"

#define DIM 3
#define BEARD_D3_AP(A1,A2) (A1 A2)
#define BEARD_D3_OP(A) A

#define BEARD_DR1(A,B)     BFAM_KRON_IXIXA    (N+1, Dr , A, B)
#define BEARD_DR1_PE(A,B)  BFAM_KRON_IXIXA_PE (N+1, Dr , A, B)
#define BEARD_DR1T(A,B)    BFAM_KRON_IXIXAT   (N+1, Dr , A, B)
#define BEARD_DR1T_PE(A,B) BFAM_KRON_IXIXAT_PE(N+1, Dr , A, B)

#define BEARD_DR2(A,B)     BFAM_KRON_IXAXI    (N+1, Dr , A, B)
#define BEARD_DR2_PE(A,B)  BFAM_KRON_IXAXI_PE (N+1, Dr , A, B)
#define BEARD_DR2T(A,B)    BFAM_KRON_IXATXI   (N+1, Dr , A, B)
#define BEARD_DR2T_PE(A,B) BFAM_KRON_IXATXI_PE(N+1, Dr , A, B)

#define BEARD_DR3(A,B)     BFAM_KRON_AXIXI   (N+1, Dr , A, B)
#define BEARD_DR3_PE(A,B)  BFAM_KRON_AXIXI_PE(N+1, Dr , A, B)
#define BEARD_DR3T(A,B)    BFAM_KRON_ATXIXI   (N+1, Dr , A, B)
#define BEARD_DR3T_PE(A,B) BFAM_KRON_ATXIXI_PE(N+1, Dr , A, B)

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
#define Nfp (NORDER+1)*(NORDER+1)
#define Nfaces  6
#endif

#else
#error "bad dimension"
#endif

#define BEARD_APPEND_4(a,b,c,d) a ## b ## c ##d
#define BEARD_APPEND_EXPAND_4(a,b,c,d) BEARD_APPEND_4(a,b,c,d)
#define beard_dgx_intra_rhs_elastic \
  BEARD_APPEND_EXPAND_4(beard_dgx_intra_rhs_elastic_,DIM,_,NORDER)

#define beard_dgx_scale_rates_elastic \
  BEARD_APPEND_EXPAND_4(beard_dgx_scale_rates_elastic_,DIM,_,NORDER)

#define beard_dgx_add_rates_elastic \
  BEARD_APPEND_EXPAND_4(beard_dgx_add_rates_elastic_,DIM,_,NORDER)

#define beard_dgx_inter_rhs_boundary \
  BEARD_APPEND_EXPAND_4(beard_dgx_inter_rhs_boundary_,DIM,_,NORDER)

#define beard_dgx_inter_rhs_interface \
  BEARD_APPEND_EXPAND_4(beard_dgx_inter_rhs_interface_,DIM,_,NORDER)

#define beard_dgx_energy \
  BEARD_APPEND_EXPAND_4(beard_dgx_energy_,DIM,_,NORDER)

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


#define BEARD_STATE beard_dgx_upwind_state_m

static inline void
beard_dgx_upwind_state_m(
          bfam_real_t *Tns,       bfam_real_t *Tps,
          bfam_real_t *vns,       bfam_real_t *vps,
          bfam_real_t  Tnm,       bfam_real_t  Tnp,
    const bfam_real_t *Tpm, const bfam_real_t *Tpp,
    const bfam_real_t  vnm, const bfam_real_t  vnp,
    const bfam_real_t *vpm, const bfam_real_t *vpp,
    const bfam_real_t  Zpm, const bfam_real_t  Zpp,
    const bfam_real_t  Zsm, const bfam_real_t  Zsp)
{
  /* upwind normal velocity and traction */
  /* wnm = Tnm - Zpm*vnm = TnS - Zpm*vnS */
  /* wnp = Tnp - Zpp*vnm = TnS + Zpp*vnS */
  bfam_real_t wnm = Tnm - Zpm*vnm;
  bfam_real_t wnp = Tnp - Zpp*vnp;
  vns[0] = (    wnp -     wnm)/(Zpp+Zpm);
  Tns[0] = (Zpm*wnp + Zpp*wnm)/(Zpp+Zpm);

  /* upwind perpendiculat velocities and tractions */
  /* wpm = Tpm - Zsm*vpm = TpS - Zsm*vpS */
  /* wpp = Tpp - Zsp*vpp =-TpS - Zsp*vpS */
  for(bfam_locidx_t i = 0; i < 3; i++)
  {
    bfam_real_t wpm = Tpm[i] - Zsm*vpm[i];
    bfam_real_t wpp = Tpp[i] - Zsp*vpp[i];
    vps[i] =-(wpm    +wpp    )/(Zsp+Zsm);
    Tps[i] = (wpm*Zsp-wpp*Zsm)/(Zsp+Zsm);
  }
}

static inline void
beard_dgx_central_state_m(
          bfam_real_t *Tns,       bfam_real_t *Tps,
          bfam_real_t *vns,       bfam_real_t *vps,
    const bfam_real_t  Tnm, const bfam_real_t  Tnp,
    const bfam_real_t *Tpm, const bfam_real_t *Tpp,
    const bfam_real_t  vnm, const bfam_real_t  vnp,
    const bfam_real_t *vpm, const bfam_real_t *vpp,
    const bfam_real_t  Zpm, const bfam_real_t  Zpp,
    const bfam_real_t  Zsm, const bfam_real_t  Zsp)
{
  vns[0] = 0.5*(vnm-vnp);
  Tns[0] = 0.5*(Tnm+Tnp);

  for(bfam_locidx_t i = 0; i < 3; i++)
  {
    vps[i] = 0.5*(vpm[i]+vpp[i]);
    Tps[i] = 0.5*(Tpm[i]-Tpp[i]);
  }
}

static inline void
beard_dgx_add_flux(const bfam_real_t scale,
    const bfam_real_t delta_TnS, const bfam_real_t* delta_TpS,
    const bfam_real_t vnS, const bfam_real_t *vpS,
    const bfam_locidx_t iM,
    bfam_real_t *dv1,  bfam_real_t *dv2,  bfam_real_t *dv3,
    bfam_real_t *dS11, bfam_real_t *dS22, bfam_real_t *dS33,
    bfam_real_t *dS12, bfam_real_t *dS13, bfam_real_t *dS23,
    const bfam_real_t lam, const bfam_real_t mu,
    const bfam_real_t rhoi, const bfam_real_t* nm,
    const bfam_real_t sJ, const bfam_real_t JI, const bfam_real_t wi)
{

  /* velocities for the flux */
  bfam_real_t vS[] = {vpS[0] + nm[0]*vnS,
                      vpS[1] + nm[1]*vnS,
                      BEARD_D3_AP(vpS[2], + nm[2]*vnS)};

  /* add the flux back in */
  bfam_real_t JI_wi_sJ = scale * wi * sJ * JI;

  dv1[iM] += rhoi*JI_wi_sJ*(           delta_TpS[0] + delta_TnS*nm[0]);
  dv2[iM] += rhoi*JI_wi_sJ*(           delta_TpS[1] + delta_TnS*nm[1]);
  dv3[iM] += rhoi*JI_wi_sJ*BEARD_D3_AP(delta_TpS[2],+ delta_TnS*nm[2]);

  const bfam_real_t M = lam+2*mu;

  dS11[iM] += JI_wi_sJ*
              BEARD_D3_AP(M  *vS[0]*nm[0] + lam*vS[1]*nm[1],+ lam*vS[2]*nm[2]);
  dS22[iM] += JI_wi_sJ*
              BEARD_D3_AP(lam*vS[0]*nm[0] + M  *vS[1]*nm[1],+ lam*vS[2]*nm[2]);
  dS33[iM] += JI_wi_sJ*
              BEARD_D3_AP(lam*vS[0]*nm[0] + lam*vS[1]*nm[1],+ M  *vS[2]*nm[2]);

  dS12[iM] += JI_wi_sJ*mu*           (vS[0]*nm[1] + vS[1]*nm[0]);
  dS13[iM] += JI_wi_sJ*mu*BEARD_D3_AP(vS[2]*nm[0],+ vS[0]*nm[2]);
  dS23[iM] += JI_wi_sJ*mu*BEARD_D3_AP(vS[2]*nm[1],+ vS[1]*nm[2]);
}

void beard_dgx_intra_rhs_elastic(
    int inN, bfam_subdomain_dgx_t *sub, const char *rate_prefix,
    const char *field_prefix, const bfam_long_real_t t)
{
  GENERIC_INIT(inN,beard_dgx_intra_rhs_elastic);

  /* get the fields we will need */
  bfam_dictionary_t *fields = &sub->base.fields;
  bfam_dictionary_t *fields_face = &sub->base.fields_face;

  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(v1 ,field_prefix,"v1" ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(v2 ,field_prefix,"v2" ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(v3 ,field_prefix,"v3" ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S11,field_prefix,"S11",fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S22,field_prefix,"S22",fields);
  BEARD_D3_OP(BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S33,field_prefix,"S33",fields));
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S12,field_prefix,"S12",fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S13,field_prefix,"S13",fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S23,field_prefix,"S23",fields);

  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(dv1 ,rate_prefix,"v1" ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(dv2 ,rate_prefix,"v2" ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(dv3 ,rate_prefix,"v3" ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(dS11,rate_prefix,"S11",fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(dS22,rate_prefix,"S22",fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(dS33,rate_prefix,"S33",fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(dS12,rate_prefix,"S12",fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(dS13,rate_prefix,"S13",fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(dS23,rate_prefix,"S23",fields);


  /* get the material properties and metric terms */
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(rhoi ,"","rho_inv"    ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(lam  ,"","lam"        ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(mu   ,"","mu"         ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Zs   ,"","Zs"         ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Zp   ,"","Zp"         ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(J    ,"","_grid_J"    ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(JI   ,"","_grid_JI"   ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Jr1x1,"","_grid_Jr0x0",fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Jr1x2,"","_grid_Jr0x1",fields);
  BEARD_D3_OP(BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Jr1x3,"","_grid_Jr0x2",fields));
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Jr2x1,"","_grid_Jr1x0",fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Jr2x2,"","_grid_Jr1x1",fields);
  BEARD_D3_OP(BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Jr2x3,"","_grid_Jr1x2",fields));

  BEARD_D3_OP(BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Jr3x1,"","_grid_Jr2x0",fields));
  BEARD_D3_OP(BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Jr3x2,"","_grid_Jr2x1",fields));
  BEARD_D3_OP(BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Jr3x3,"","_grid_Jr2x2",fields));

  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(n1,"","_grid_nx0",fields_face);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(n2,"","_grid_nx1",fields_face);
  BEARD_D3_OP(BFAM_LOAD_FIELD_RESTRICT_ALIGNED(n3,"","_grid_nx2",fields_face));
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(sJ  ,"","_grid_sJ",fields_face);

  bfam_locidx_t K  = sub->K;
  BFAM_ALIGN(32) bfam_real_t aux1[Np];
  BFAM_ALIGN(32) bfam_real_t aux2[Np];
  BEARD_D3_OP(BFAM_ALIGN(32) bfam_real_t aux3[Np]);

  BFAM_ALIGN(32) bfam_real_t Dr1Tv1[Np];
  BFAM_ALIGN(32) bfam_real_t Dr2Tv1[Np];
  BEARD_D3_OP(BFAM_ALIGN(32) bfam_real_t Dr3Tv1[Np]);
  BFAM_ALIGN(32) bfam_real_t Dr1Tv2[Np];
  BFAM_ALIGN(32) bfam_real_t Dr2Tv2[Np];
  BEARD_D3_OP(BFAM_ALIGN(32) bfam_real_t Dr3Tv2[Np]);
  BFAM_ALIGN(32) bfam_real_t Dr1Tv3[Np];
  BFAM_ALIGN(32) bfam_real_t Dr2Tv3[Np];
  BEARD_D3_OP(BFAM_ALIGN(32) bfam_real_t Dr3Tv3[Np]);
  BFAM_ALIGN(32) bfam_real_t rate[Np];

  BFAM_ALIGN(32) bfam_real_t Mv1[Np];
  BFAM_ALIGN(32) bfam_real_t Mv2[Np];
  BFAM_ALIGN(32) bfam_real_t Mv3[Np];

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

    /* v1 +=
     *      JI*rhoi*(Dr*(Jrx*S11 + Jry*S12 + Jrz*S13))
     *     +JI*rhoi*(Ds*(Jsx*S11 + Jsy*S12 + Jsz*S13))
     *     +JI*rhoi*(Dt*(Jtx*S11 + Jty*S12 + Jtz*S13))
     */
    for(bfam_locidx_t i = 0; i < Np; i++)
    {
      aux1[i] =
        BEARD_D3_AP(Jr1x1[off+i] * S11[off+i] + Jr1x2[off+i] * S12[off+i],
                  + Jr1x3[off+i] * S13[off+i]);
      aux2[i] =
        BEARD_D3_AP(Jr2x1[off+i] * S11[off+i] + Jr2x2[off+i] * S12[off+i],
                  + Jr2x3[off+i] * S13[off+i]);
      BEARD_D3_OP(
          aux3[i] = Jr3x1[off+i] * S11[off+i] + Jr3x2[off+i] * S12[off+i]
                  + Jr3x3[off+i] * S13[off+i]);
    }
    BEARD_DR1   (aux1   , rate); /* Dr */
    BEARD_DR2_PE(aux2   , rate); /* Ds */
    BEARD_DR3_PE(aux3   , rate); /* Dt */

    for(bfam_locidx_t i = 0; i < Np; i++)
      dv1[off+i] += JI[off+i]*rhoi[off+i]*rate[i];

    /* v2 +=
     *      JI*rhoi*(Dr*(Jrx*S12 + Jry*S22 + Jrz*S23))
     *     +JI*rhoi*(Ds*(Jsx*S12 + Jsy*S22 + Jsz*S23))
     *     +JI*rhoi*(Dt*(Jtx*S12 + Jty*S22 + Jtz*S23))
     */
    for(bfam_locidx_t i = 0; i < Np; i++)
    {
      aux1[i] =
        BEARD_D3_AP(Jr1x1[off+i] * S12[off+i] + Jr1x2[off+i] * S22[off+i],
                  + Jr1x3[off+i] * S23[off+i]);
      aux2[i] =
        BEARD_D3_AP(Jr2x1[off+i] * S12[off+i] + Jr2x2[off+i] * S22[off+i],
                  + Jr2x3[off+i] * S23[off+i]);
      BEARD_D3_OP(
          aux3[i] = Jr3x1[off+i] * S12[off+i] + Jr3x2[off+i] * S22[off+i]
                  + Jr3x3[off+i] * S23[off+i]);
    }
    BEARD_DR1   (aux1   , rate); /* Dr */
    BEARD_DR2_PE(aux2   , rate); /* Ds */
    BEARD_DR3_PE(aux3   , rate); /* Dt */

    for(bfam_locidx_t i = 0; i < Np; i++)
      dv2[off+i] += JI[off+i]*rhoi[off+i]*rate[i];


    /* v3 +=
     *      JI*rhoi*(Dr*(Jrx*S13 + Jry*S23 + Jrz*S33))
     *     +JI*rhoi*(Ds*(Jsx*S13 + Jsy*S23 + Jsz*S33))
     *     +JI*rhoi*(Dt*(Jtx*S13 + Jty*S23 + Jtz*S33))
     */
    for(bfam_locidx_t i = 0; i < Np; i++)
    {
      aux1[i] =
        BEARD_D3_AP(Jr1x1[off+i] * S13[off+i] + Jr1x2[off+i] * S23[off+i],
                  + Jr1x3[off+i] * S33[off+i]);
      aux2[i] =
        BEARD_D3_AP(Jr2x1[off+i] * S13[off+i] + Jr2x2[off+i] * S23[off+i],
                  + Jr2x3[off+i] * S33[off+i]);
      BEARD_D3_OP(
          aux3[i] = Jr3x1[off+i] * S13[off+i] + Jr3x2[off+i] * S23[off+i]
                  + Jr3x3[off+i] * S33[off+i]);
    }
    BEARD_DR1   (aux1   , rate); /* Dr */
    BEARD_DR2_PE(aux2   , rate); /* Ds */
    BEARD_DR3_PE(aux3   , rate); /* Dt */

    for(bfam_locidx_t i = 0; i < Np; i++)
      dv3[off+i] += JI[off+i]*rhoi[off+i]*rate[i];


   /* we use these a lot so store them */
    bfam_locidx_t n = 0;
#if DIM==3
    for(bfam_locidx_t k = 0; k < N+1; k++)
#endif
      for(bfam_locidx_t j = 0; j < N+1; j++)
        for(bfam_locidx_t i = 0; i < N+1; i++,n++)
        {
          Mv1[n] = BEARD_D3_AP(w[i]*w[j],*w[k])*v1[off+n];
          Mv2[n] = BEARD_D3_AP(w[i]*w[j],*w[k])*v2[off+n];
          Mv3[n] = BEARD_D3_AP(w[i]*w[j],*w[k])*v3[off+n];
        }

    /*
     * S11 += -(lam+2*mu)*MI*JI*((Jrx*Dr' + Jsx*Ds' + Jtx*Dt')*v1
     *        -lam       *MI*JI*((Jry*Dr' + Jsy*Ds' + Jty*Dt')*v2
     *        -lam       *MI*JI*((Jrz*Dr' + Jsz*Ds' + Jtz*Dt')*v3
     *
     * S22 += -lam       *MI*JI*((Jrx*Dr' + Jsx*Ds' + Jtx*Dt')*v1
     *        -(lam+2*mu)*MI*JI*((Jry*Dr' + Jsy*Ds' + Jty*Dt')*v2
     *        -lam       *MI*JI*((Jrz*Dr' + Jsz*Ds' + Jtz*Dt')*v3
     *
     * S33 += -lam l     *MI*JI*((Jrx*Dr' + Jsx*Ds' + Jtx*Dt')*v1
     *        -lam l     *MI*JI*((Jry*Dr' + Jsy*Ds' + Jty*Dt')*v2
     *        -(lam+2*mu)*MI*JI*((Jrz*Dr' + Jsz*Ds' + Jtz*Dt')*v3
     *
     * S12 += -mu*MI*JI*((Jrx*Dr' + Jsx*Ds' + Jtx*Dt')*v2
     *        -mu*MI*JI*((Jry*Dr' + Jsy*Ds' + Jty*Dt')*v1
     *
     * S13 += -mu*MI*JI*((Jrx*Dr' + Jsx*Ds' + Jtx*Dt')*v3
     *        -mu*MI*JI*((Jrz*Dr' + Jsz*Ds' + Jtz*Dt')*v1
     *
     * S23 += -mu*MI*JI*((Jrx*Dr' + Jsx*Ds' + Jtx*Dt')*v3
     *        -mu*MI*JI*((Jrz*Dr' + Jsz*Ds' + Jtz*Dt')*v2
     */
    BEARD_DR1T(Mv1, Dr1Tv1);
    BEARD_DR2T(Mv1, Dr2Tv1);
    BEARD_DR3T(Mv1, Dr3Tv1);
    BEARD_DR1T(Mv2, Dr1Tv2);
    BEARD_DR2T(Mv2, Dr2Tv2);
    BEARD_DR3T(Mv2, Dr3Tv2);
    BEARD_DR1T(Mv3, Dr1Tv3);
    BEARD_DR2T(Mv3, Dr2Tv3);
    BEARD_DR3T(Mv3, Dr3Tv3);

    n = 0;
#if DIM==3
    for(bfam_locidx_t k = 0; k < N+1; k++)
#endif
      for(bfam_locidx_t j = 0; j < N+1; j++)
        for(bfam_locidx_t i = 0; i < N+1; i++,n++)
        {
          const bfam_real_t wi_ijk = BEARD_D3_AP(w[i]*w[j],*w[k]);
          const bfam_locidx_t nid = off+n;

          const bfam_real_t lam_l = lam[nid];
          const bfam_real_t mu_l  = mu[nid];
          const bfam_real_t M_l   = lam_l+2*mu_l;

          dS11[nid] -= wi_ijk*JI[nid]*(
              +M_l  *BEARD_D3_AP(Jr1x1[nid]*Dr1Tv1[n]
                                +Jr2x1[nid]*Dr2Tv1[n],
                                +Jr3x1[nid]*Dr3Tv1[n])
              +lam_l*BEARD_D3_AP(Jr1x2[nid]*Dr1Tv2[n]
                                +Jr2x2[nid]*Dr2Tv2[n],
                                +Jr3x2[nid]*Dr3Tv2[n])
#if DIM==3
              +lam_l*BEARD_D3_AP(Jr1x3[nid]*Dr1Tv3[n]
                                +Jr2x3[nid]*Dr2Tv3[n],
                                +Jr3x3[nid]*Dr3Tv3[n])
#endif
              );

          dS22[nid] -= wi_ijk*JI[nid]*(
              +lam_l*BEARD_D3_AP(Jr1x1[nid]*Dr1Tv1[n]
                                +Jr2x1[nid]*Dr2Tv1[n],
                                +Jr3x1[nid]*Dr3Tv1[n])
              +M_l  *BEARD_D3_AP(Jr1x2[nid]*Dr1Tv2[n]
                                +Jr2x2[nid]*Dr2Tv2[n],
                                +Jr3x2[nid]*Dr3Tv2[n])
#if DIM==3
              +lam_l*BEARD_D3_AP(Jr1x3[nid]*Dr1Tv3[n]
                                +Jr2x3[nid]*Dr2Tv3[n],
                                +Jr3x3[nid]*Dr3Tv3[n])
#endif
              );

          dS33[nid] -= wi_ijk*JI[nid]*(
              +lam_l*BEARD_D3_AP(Jr1x1[nid]*Dr1Tv1[n]
                                +Jr2x1[nid]*Dr2Tv1[n],
                                +Jr3x1[nid]*Dr3Tv1[n])
              +lam_l*BEARD_D3_AP(Jr1x2[nid]*Dr1Tv2[n]
                                +Jr2x2[nid]*Dr2Tv2[n],
                                +Jr3x2[nid]*Dr3Tv2[n])
#if DIM==3
              +M_l  *BEARD_D3_AP(Jr1x3[nid]*Dr1Tv3[n]
                                 +Jr2x3[nid]*Dr2Tv3[n],
                                 +Jr3x3[nid]*Dr3Tv3[n])
#endif
              );



          dS12[nid] -= wi_ijk*JI[nid]*(
              +mu_l*BEARD_D3_AP(Jr1x1[nid]*Dr1Tv2[n]
                               +Jr2x1[nid]*Dr2Tv2[n],
                               +Jr3x1[nid]*Dr3Tv2[n])
              +mu_l*BEARD_D3_AP(Jr1x2[nid]*Dr1Tv1[n]
                               +Jr2x2[nid]*Dr2Tv1[n],
                               +Jr3x2[nid]*Dr3Tv1[n])
              );

          dS13[nid] -= wi_ijk*JI[nid]*(
              +mu_l*BEARD_D3_AP(Jr1x1[nid]*Dr1Tv3[n]
                               +Jr2x1[nid]*Dr2Tv3[n],
                               +Jr3x1[nid]*Dr3Tv3[n])
#if DIM==3
              +mu_l*BEARD_D3_AP(Jr1x3[nid]*Dr1Tv1[n]
                               +Jr2x3[nid]*Dr2Tv1[n],
                               +Jr3x3[nid]*Dr3Tv1[n])
#endif
              );


          dS23[nid] -= wi_ijk*JI[nid]*(
              +mu_l*BEARD_D3_AP(Jr1x2[nid]*Dr1Tv3[n]
                               +Jr2x2[nid]*Dr2Tv3[n],
                               +Jr3x2[nid]*Dr3Tv3[n])
#if DIM==3
              +mu_l*BEARD_D3_AP(Jr1x3[nid]*Dr1Tv2[n]
                               +Jr2x3[nid]*Dr2Tv2[n],
                               +Jr3x3[nid]*Dr3Tv2[n])
#endif
              );
      }

    /* loop over faces */
    for(bfam_locidx_t face = 0; face < Nfaces; face++)
    {
      for(bfam_locidx_t pnt = 0; pnt < Nfp; pnt++)
      {
        const bfam_locidx_t f = pnt + Nfp*(face + Nfaces*e);
        const bfam_locidx_t iM = vmapM[f];
        const bfam_locidx_t iP = vmapP[f];

        /* Setup stuff for the minus side */
        const bfam_real_t ZsM = Zs[iM];
        const bfam_real_t ZpM = Zp[iM];

        const bfam_real_t nM[] = {n1[f],n2[f],BEARD_D3_AP(0,+n3[f])};

        bfam_real_t TpM[] = {
          BEARD_D3_AP(nM[0]*S11[iM] + nM[1]*S12[iM], + nM[2]*S13[iM]),
          BEARD_D3_AP(nM[0]*S12[iM] + nM[1]*S22[iM], + nM[2]*S23[iM]),
          BEARD_D3_AP(nM[0]*S13[iM] + nM[1]*S23[iM], + nM[2]*S33[iM]),
        };
        const bfam_real_t TnM = BEARD_D3_AP(TpM[0]*nM[0]
                                           +TpM[1]*nM[1],
                                           +TpM[2]*nM[2]);
        TpM[0] = TpM[0]-TnM*nM[0];
        TpM[1] = TpM[1]-TnM*nM[1];
        BEARD_D3_OP(TpM[2] = TpM[2]-TnM*nM[2]);

        bfam_real_t vpM[] = {v1[iM],v2[iM],v3[iM]};
        const bfam_real_t vnM = BEARD_D3_AP(nM[0]*vpM[0]
                                           +nM[1]*vpM[1],
                                           +nM[2]*vpM[2]);
        vpM[0] = vpM[0]-vnM*nM[0];
        vpM[1] = vpM[1]-vnM*nM[1];
        BEARD_D3_OP(vpM[2] = vpM[2]-vnM*nM[2]);

        /* Setup stuff for the plus side */
        bfam_real_t ZsP = Zs[iP];
        bfam_real_t ZpP = Zp[iP];

        bfam_real_t nP[] = {-nM[0],-nM[1],-nM[2]};

        bfam_real_t TpP[] = {
          BEARD_D3_AP(nP[0]*S11[iP] + nP[1]*S12[iP], + nP[2]*S13[iP]),
          BEARD_D3_AP(nP[0]*S12[iP] + nP[1]*S22[iP], + nP[2]*S23[iP]),
          BEARD_D3_AP(nP[0]*S13[iP] + nP[1]*S23[iP], + nP[2]*S33[iP]),
        };
        const bfam_real_t TnP = BEARD_D3_AP(TpP[0]*nP[0]
                                           +TpP[1]*nP[1],
                                           +TpP[2]*nP[2]);
        TpP[0] = TpP[0]-TnP*nP[0];
        TpP[1] = TpP[1]-TnP*nP[1];
        BEARD_D3_OP(TpP[2] = TpP[2]-TnP*nP[2]);

        bfam_real_t vpP[] = {v1[iP],v2[iP],v3[iP]};
        const bfam_real_t vnP = BEARD_D3_AP(nP[0]*vpP[0]
                                           +nP[1]*vpP[1],
                                           +nP[2]*vpP[2]);
        vpP[0] = vpP[0]-vnP*nP[0];
        vpP[1] = vpP[1]-vnP*nP[1];
        BEARD_D3_OP(vpP[2] = vpP[2]-vnP*nP[2]);

        bfam_real_t TnS;
        bfam_real_t TpS[3];
        bfam_real_t vnS;
        bfam_real_t vpS[3];

        BEARD_STATE(&TnS,TpS,&vnS,vpS,
            TnM, TnP, TpM, TpP, vnM, vnP, vpM, vpP, ZpM, ZpP, ZsM, ZsP);

        TnS    -= TnM;
        TpS[0] -= TpM[0];
        TpS[1] -= TpM[1];
        TpS[2] -= TpM[2];

        /* intra */
        beard_dgx_add_flux(1, TnS,TpS,vnS,vpS,iM,
            dv1,dv2,dv3, dS11,dS22,dS33,dS12,dS13,dS23,
            lam[iM],mu[iM],rhoi[iM],nM,sJ[f],JI[iM],wi[0]);
      }
    }
  }
}

void beard_dgx_scale_rates_elastic(
    int inN, bfam_subdomain_dgx_t *sub, const char *rate_prefix,
    const bfam_long_real_t a)
{
  GENERIC_INIT(inN,beard_dgx_scale_rates_elastic);

  const bfam_locidx_t num_pts = sub->K * Np;
  bfam_dictionary_t *fields = &sub->base.fields;
  const char *f_names[] =
    {"v1","v2","v3","S11","S22","S33","S12","S13","S23",NULL};
  for(bfam_locidx_t f = 0; f_names[f]!=NULL; ++f)
  {
    BFAM_LOAD_FIELD_RESTRICT_ALIGNED(rate,rate_prefix,f_names[f],fields);
    for(bfam_locidx_t n = 0; n < num_pts; n++) rate[n] *= a;
  }
}

void beard_dgx_add_rates_elastic(
    int inN, bfam_subdomain_dgx_t *sub, const char *field_prefix_lhs,
    const char *field_prefix_rhs, const char *rate_prefix,
    const bfam_long_real_t a)
{
  GENERIC_INIT(inN,beard_dgx_add_rates_elastic);

  const char *f_names[] = {"v1", "v2", "v3",
    "S11", "S22", "S33", "S12", "S13", "S23", NULL};

  const bfam_locidx_t num_pts = sub->K * Np;

  /* get the fields we will need */
  bfam_dictionary_t *fields = &sub->base.fields;

  for(bfam_locidx_t f = 0; f_names[f] != NULL; f++)
  {
    BFAM_LOAD_FIELD_ALIGNED(         lhs ,field_prefix_lhs,f_names[f],fields);
    BFAM_LOAD_FIELD_ALIGNED(         rhs ,field_prefix_rhs,f_names[f],fields);
    BFAM_LOAD_FIELD_RESTRICT_ALIGNED(rate,rate_prefix     ,f_names[f],fields);
    for(bfam_locidx_t n = 0; n < num_pts; n++) lhs[n] = rhs[n] + a*rate[n];
  }
}

void beard_dgx_inter_rhs_boundary(
    int inN, bfam_subdomain_dgx_t *sub_g, const char *rate_prefix,
    const char *field_prefix, const bfam_long_real_t t, const bfam_real_t R)
{
  GENERIC_INIT(inN,beard_dgx_inter_rhs_boundary);

}

void beard_dgx_inter_rhs_interface(
    int inN, bfam_subdomain_dgx_t *sub_g, const char *rate_prefix,
    const char *field_prefix, const bfam_long_real_t t)
{
  GENERIC_INIT(inN,beard_dgx_inter_rhs_interface);

}

void beard_dgx_energy(
    int inN, bfam_real_t *energy_sq,
    bfam_subdomain_dgx_t *sub, const char *field_prefix)
{
  GENERIC_INIT(inN,beard_dgx_energy);

  /* get the fields we will need */
  bfam_dictionary_t *fields = &sub->base.fields;
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(v1 ,field_prefix,"v1" ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(v2 ,field_prefix,"v2" ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(v3 ,field_prefix,"v3" ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S11,field_prefix,"S11",fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S22,field_prefix,"S22",fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S33,field_prefix,"S33",fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S12,field_prefix,"S12",fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S13,field_prefix,"S13",fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S23,field_prefix,"S23",fields);

  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(rho,"","rho"    ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(lam,"","lam"    ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(mu ,"","mu"     ,fields);
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

          /* setup up deviatoric stress tensor */
          bfam_real_t mean_stress = (S11[n]+S22[n]+S33[n])/3.0;
          bfam_real_t s11 = S11[n] - mean_stress;
          bfam_real_t s22 = S22[n] - mean_stress;
          bfam_real_t s33 = S33[n] - mean_stress;
          bfam_real_t s12 = S12[n];
          bfam_real_t s13 = S13[n];
          bfam_real_t s23 = S23[n];

          /* bulk modulus */
          bfam_real_t K = lam[n] + 2.0*mu[n]/3.0;

          energy_sq[0] += BEARD_D3_AP(wi*wj,*wk)*J[n]*(
              rho[n]*(v1[n]*v1[n] + v2[n]*v2[n] + v3[n]*v3[n])/2
              + (s11*s11 + s22*s22 + s33*s33
                + 2*s12*s12 + 2*s13*s13 + 2*s23*s23)/(4*mu[n])
              + mean_stress*mean_stress/(2*K)
              );
        }
      }
#if DIM==3
    }
#endif
  }
}

#endif
