#include "beard_dgx_rhs.h"

#ifndef NORDER
#define NORDER
#define USE_GENERIC
#else
#define N    NORDER
#define Np  (NORDER+1)*(NORDER+1)
#define Nfp (NORDER+1)
#endif

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


void BFAM_APPEND_EXPAND(beard_dgx_print_order_,NORDER)(int inN)
{

#ifdef USE_GENERIC
  const int N = inN;
  BFAM_WARNING("Using generic order print function");
#endif

  BFAM_INFO("N = %d", N);
}

static inline void
project_flux(bfam_real_t *Tns,       bfam_real_t *Tps,
             bfam_real_t *vns,       bfam_real_t *vps,
        bfam_locidx_t  Nfp_in,     bfam_locidx_t  Npg,
       const bfam_real_t *Tng, const bfam_real_t *Tpg,
       const bfam_real_t *vng, const bfam_real_t *vpg,
       const bfam_real_t * MP, const bfam_real_t  *wi)
{
#ifdef USE_GENERIC
  const int Nfp = Nfp_in;
#endif
  /* zero out the components */
  for(int i = 0; i < Nfp; i++)
  {
    for(int k = 0; k < 3; k++)
    {
      Tps[3*i+k] = 0;
      vps[3*i+k] = 0;
    }
    Tns[i] = 0;
    vns[i] = 0;
  }
  for(int j = 0; j < Npg; j++)
  {
    for(int i = 0; i < Nfp; i++)
    {
      for(int k = 0; k < 3; k++)
      {
        Tps[3*i+k] += wi[i]*MP[i+j*Nfp]*Tpg[3*j+k];
        vps[3*i+k] += wi[i]*MP[i+j*Nfp]*vpg[3*j+k];
      }
      Tns[i] += wi[i]*MP[i+j*Nfp]*Tng[j];
      vns[i] += wi[i]*MP[i+j*Nfp]*vng[j];
    }
  }
}

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
beard_dgx_upwind_state_friction_m(
          bfam_real_t *Tps,       bfam_real_t *vpsm,        bfam_real_t *Vps,
    const bfam_real_t    T,
    const bfam_real_t *Tpm, const bfam_real_t *Tpp  , const bfam_real_t *Tp0,
    const bfam_real_t *vpm, const bfam_real_t *vpp  ,
    const bfam_real_t  Zsm, const bfam_real_t  Zsp  )
{
  /* upwind perpendiculat velocities and tractions */
  /* wpm = Tpm - Zsm*vpm = TpS - Zsm*vpS */
  /* wpp = Tpp - Zsp*vpp =-TpS - Zsp*vpS */
  bfam_real_t phi[3];
  bfam_real_t mag = 0;
  for(bfam_locidx_t i = 0; i < 3; i++)
  {
    bfam_real_t wpm = Tpm[i] + Tp0[i] - Zsm*vpm[i];
    bfam_real_t wpp = Tpp[i] - Tp0[i] - Zsp*vpp[i];
    phi[i] = (wpm*Zsp-wpp*Zsm)/(Zsp+Zsm);
    mag += phi[i]*phi[i];
  }
  mag = BFAM_REAL_SQRT(mag);
  for(bfam_locidx_t i = 0; i < 3; i++)
  {
    Tps[i]           = T*phi[i]/mag - Tp0[i];
    vpsm[i]          = (Tps[i]-Tpm[i])/Zsm + vpm[i];
    bfam_real_t vpsp =-(Tps[i]+Tpp[i])/Zsp + vpp[i];
    Vps[i]           = vpsm[i]-vpsp;
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
  bfam_real_t vS[] = {vpS[0] + nm[0]*vnS,vpS[1] + nm[1]*vnS,vpS[2]};

  /* add the flux back in */
  bfam_real_t JI_wi_sJ = scale * wi * sJ * JI;

  dv1[iM] += rhoi*JI_wi_sJ*(delta_TpS[0]+delta_TnS*nm[0]);
  dv2[iM] += rhoi*JI_wi_sJ*(delta_TpS[1]+delta_TnS*nm[1]);
  dv3[iM] += rhoi*JI_wi_sJ*(delta_TpS[2]);

  dS11[iM] += JI_wi_sJ*((lam+2*mu)*vS[0]*nm[0] +  lam      *vS[1]*nm[1]);
  dS22[iM] += JI_wi_sJ*( lam      *vS[0]*nm[0] + (lam+2*mu)*vS[1]*nm[1]);
  dS33[iM] += JI_wi_sJ*( lam      *vS[0]*nm[0] +  lam      *vS[1]*nm[1]);
  dS12[iM] += JI_wi_sJ*mu*(vS[0]*nm[1] + vS[1]*nm[0]);
  dS13[iM] += JI_wi_sJ*mu *vS[2]*nm[0];
  dS23[iM] += JI_wi_sJ*mu *vS[2]*nm[1];
}

static inline void
beard_dgx_remove_flux( const int Nfp_in,
    bfam_locidx_t face, bfam_locidx_t e, const bfam_locidx_t *vmapM,
    const bfam_real_t *n1, const bfam_real_t *n2,
    const bfam_real_t *Zs, const bfam_real_t *Zp,
    const bfam_real_t *mu, const bfam_real_t *rhoi, const bfam_real_t *lam,
    const bfam_real_t *sJ, const bfam_real_t *JI,   const bfam_real_t *wi,
    const bfam_real_t *v1,  const bfam_real_t *v2,  const bfam_real_t *v3,
    const bfam_real_t *S11, const bfam_real_t *S22, const bfam_real_t *S33,
    const bfam_real_t *S12, const bfam_real_t *S13, const bfam_real_t *S23,
          bfam_real_t *dv1,  bfam_real_t *dv2,  bfam_real_t *dv3,
          bfam_real_t *dS11, bfam_real_t *dS22, bfam_real_t *dS33,
          bfam_real_t *dS12, bfam_real_t *dS13, bfam_real_t *dS23
    )
{
#ifdef USE_GENERIC
  const int Nfp = Nfp_in;
#endif
  for(bfam_locidx_t pnt = 0; pnt < Nfp; pnt++)
  {
    bfam_locidx_t f = pnt + Nfp*(face + 4*e);
    bfam_locidx_t iM = vmapM[f];

    /* Setup stuff for the minus side */
    bfam_real_t Zsm = Zs[iM];
    bfam_real_t Zpm = Zp[iM];

    bfam_real_t nm[] = {n1[f],n2[f],0};

    bfam_real_t Tpm[] = {
      nm[0]*S11[iM]+nm[1]*S12[iM],
      nm[0]*S12[iM]+nm[1]*S22[iM],
      nm[0]*S13[iM]+nm[1]*S23[iM],
    };
    bfam_real_t Tnm = Tpm[0]*nm[0]+Tpm[1]*nm[1];
    Tpm[0] = Tpm[0]-Tnm*nm[0];
    Tpm[1] = Tpm[1]-Tnm*nm[1];

    bfam_real_t vpm[] = {v1[iM],v2[iM],v3[iM]};
    bfam_real_t vnm = nm[0]*vpm[0]+nm[1]*vpm[1];
    vpm[0] = vpm[0]-vnm*nm[0];
    vpm[1] = vpm[1]-vnm*nm[1];

    /* First remove what we already did */
    // bfam_real_t Zsp = Zs[iM];
    // bfam_real_t Zpp = Zp[iM];
    bfam_real_t Zsp = Zsm;
    bfam_real_t Zpp = Zpm;

    bfam_real_t Tpp[] = {-Tpm[0],-Tpm[1],-Tpm[2]};
    bfam_real_t Tnp   = Tnm;

    bfam_real_t vpp[] = { vpm[0], vpm[1], vpm[2]};
    bfam_real_t vnp   = -vnm;

    bfam_real_t TnS;
    bfam_real_t TpS[3];
    bfam_real_t vnS;
    bfam_real_t vpS[3];

    beard_dgx_upwind_state_m(&TnS,TpS,&vnS,vpS,
        Tnm, Tnp, Tpm, Tpp, vnm, vnp, vpm, vpp, Zpm, Zpp, Zsm, Zsp);

    TnS -= Tnm;
    TpS[0] -= Tpm[0];
    TpS[1] -= Tpm[1];
    TpS[2] -= Tpm[2];

    beard_dgx_add_flux(-1, TnS,TpS,vnS,vpS,iM,
        dv1,dv2,dv3, dS11,dS22,dS33,dS12,dS13,dS23,
        lam[iM],mu[iM],rhoi[iM],nm,sJ[f],JI[iM],wi[0]);
  }
}



void BFAM_APPEND_EXPAND(beard_dgx_intra_rhs_elastic_,NORDER)(
    int inN, bfam_subdomain_dgx_quad_t *sub, const char *rate_prefix,
    const char *field_prefix, const bfam_long_real_t t)
{

#ifdef USE_GENERIC
  const int N   = inN;
  const int Np  = (inN+1)*(inN+1);
  const int Nfp = inN+1;
  BFAM_WARNING("Using generic intra rhs function");
#endif

  /* get the fields we will need */
  bfam_dictionary_t *fields = &sub->base.fields;
  bfam_dictionary_t *fields_face = &sub->base.fields_face;

  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(v1 ,field_prefix,"v1" ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(v2 ,field_prefix,"v2" ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(v3 ,field_prefix,"v3" ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S11,field_prefix,"S11",fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S22,field_prefix,"S22",fields);
  /*
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S33,field_prefix,"S33",fields);
  */
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
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(rhoi,"","rho_inv"  ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(lam ,"","lam"      ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(mu  ,"","mu"       ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Zs  ,"","Zs"       ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Zp  ,"","Zp"       ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(J   ,"","_grid_J"  ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(JI  ,"","_grid_JI" ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Jrx ,"","_grid_Jrx",fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Jry ,"","_grid_Jry",fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Jsx ,"","_grid_Jsx",fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Jsy ,"","_grid_Jsy",fields);

  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(n1  ,"","_grid_nx",fields_face);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(n2  ,"","_grid_ny",fields_face);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(sJ  ,"","_grid_sJ",fields_face);

  bfam_locidx_t K  = sub->K;
  BFAM_ALIGN(32) bfam_real_t aux1[Np];
  BFAM_ALIGN(32) bfam_real_t aux2[Np];
  BFAM_ALIGN(32) bfam_real_t DrTv1[Np];
  BFAM_ALIGN(32) bfam_real_t DsTv1[Np];
  BFAM_ALIGN(32) bfam_real_t DrTv2[Np];
  BFAM_ALIGN(32) bfam_real_t DsTv2[Np];
  BFAM_ALIGN(32) bfam_real_t DrTv3[Np];
  BFAM_ALIGN(32) bfam_real_t DsTv3[Np];
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

    /* v1 += JI*rhoi*(Dr*(Jrx*S11+Jry*S12)) + JI*rhoi*(Ds*(Jsx*S11+Jsy*S12));*/
    for(bfam_locidx_t i = 0; i < Np; i++)
    {
      aux1[i] = Jrx[off+i] * S11[off+i] + Jry[off+i] * S12[off+i];
      aux2[i] = Jsx[off+i] * S11[off+i] + Jsy[off+i] * S12[off+i];
    }
    BFAM_KRON_IXA   (N+1, Dr     , aux1   , rate); /* Dr */
    BFAM_KRON_AXI_PE(N+1, Dr     , aux2   , rate); /* Ds */

    for(bfam_locidx_t i = 0; i < Np; i++)
      dv1[off+i] += JI[off+i]*rhoi[off+i]*rate[i];

    /* v2 += JI*rhoi*(Dr*(Jrx*S12+Jry*S22)) + JI*rhoi*(Ds*(Jsx*S12+Jsy*S22));*/
    for(bfam_locidx_t i = 0; i < Np; i++)
    {
      aux1[i] = Jrx[off+i] * S12[off+i] + Jry[off+i] * S22[off+i];
      aux2[i] = Jsx[off+i] * S12[off+i] + Jsy[off+i] * S22[off+i];
    }
    BFAM_KRON_IXA   (N+1, Dr     , aux1   , rate); /* Dr */
    BFAM_KRON_AXI_PE(N+1, Dr     , aux2   , rate); /* Ds */

    for(bfam_locidx_t i = 0; i < Np; i++)
      dv2[off+i] += JI[off+i]*rhoi[off+i]*rate[i];

    /* rho*v3 += JI*rhoi*(Dr*(Jrx*S13+Jry*S23)) +  JI*rhoi*(Ds*(Jsx*S13+Jsy*S23));*/
    for(bfam_locidx_t i = 0; i < Np; i++)
    {
      aux1[i] = Jrx[off+i] * S13[off+i] + Jry[off+i] * S23[off+i];
      aux2[i] = Jsx[off+i] * S13[off+i] + Jsy[off+i] * S23[off+i];
    }
    BFAM_KRON_IXA   (N+1, Dr     , aux1   , rate); /* Dr */
    BFAM_KRON_AXI_PE(N+1, Dr     , aux2   , rate); /* Ds */

    for(bfam_locidx_t i = 0; i < Np; i++)
      dv3[off+i] += JI[off+i]*rhoi[off+i]*rate[i];

    /* we use these a lot so store them */
    for(bfam_locidx_t i = 0; i < N+1; i++)
      for(bfam_locidx_t j = 0; j < N+1; j++)
      {
        bfam_locidx_t k = i*(N+1)+j;
        Mv1[k] = w[i]*w[j]*v1[off+k];
        Mv2[k] = w[i]*w[j]*v2[off+k];
        Mv3[k] = w[i]*w[j]*v3[off+k];
      }

    /*
     * S11 += -(lam+2*mu)*MI*JI*((Jrx*Dr'+Jsx*Ds')*v1
     *        -lam*MI*JI*((Jry*Dr'+Jsy*Ds')*v2
     *
     * S22 += -lam*MI*JI*((Jrx*Dr'+Jsx*Ds')*v1
     *        -(lam+2*mu)*MI*JI*((Jry*Dr'+Jsy*Ds')*v2
     *
     * S33 += -lam*MI*JI*((Jrx*Dr'+Jsx*Ds')*v1
     *        -lam*MI*JI*((Jry*Dr'+Jsy*Ds')*v2
     *
     * S12 += -mu*MI*JI*((Jrx*Dr'+Jsx*Ds')*v2
     *        -mu*MI*JI*((Jry*Dr'+Jsy*Ds')*v1
     *
     * S13 += -mu*MI*JI*((Jrx*Dr'+Jsx*Ds')*v3
     *
     * S23 += -mu*MI*JI*((Jry*Dr'+Jsy*Ds')*v3
     */
    BFAM_KRON_IXAT(N+1, Dr, Mv1, DrTv1);
    BFAM_KRON_ATXI(N+1, Dr, Mv1, DsTv1);
    BFAM_KRON_IXAT(N+1, Dr, Mv2, DrTv2);
    BFAM_KRON_ATXI(N+1, Dr, Mv2, DsTv2);
    BFAM_KRON_IXAT(N+1, Dr, Mv3, DrTv3);
    BFAM_KRON_ATXI(N+1, Dr, Mv3, DsTv3);
    for(bfam_locidx_t i = 0; i < N+1; i++)
      for(bfam_locidx_t j = 0; j < N+1; j++)
      {
        bfam_locidx_t k = i*(N+1)+j;

        dS11[off+k] -= wi[i]*wi[j]*JI[off+k]*(
            (lam[off+k]+2*mu[off+k])*(Jrx[off+k]*DrTv1[k]+Jsx[off+k]*DsTv1[k])
            +lam[off+k]             *(Jry[off+k]*DrTv2[k]+Jsy[off+k]*DsTv2[k]));

        dS22[off+k] -= wi[i]*wi[j]*JI[off+k]*(
            + lam[off+k]             *(Jrx[off+k]*DrTv1[k]+Jsx[off+k]*DsTv1[k])
            +(lam[off+k]+2*mu[off+k])*(Jry[off+k]*DrTv2[k]+Jsy[off+k]*DsTv2[k])
            );

        dS33[off+k] -= wi[i]*wi[j]*JI[off+k]*(
             lam[off+k]             *(Jrx[off+k]*DrTv1[k]+Jsx[off+k]*DsTv1[k])
            +lam[off+k]             *(Jry[off+k]*DrTv2[k]+Jsy[off+k]*DsTv2[k]));

        dS12[off+k] -= wi[i]*wi[j]*JI[off+k]*(
             mu[off+k]             *(Jrx[off+k]*DrTv2[k]+Jsx[off+k]*DsTv2[k])
            +mu[off+k]             *(Jry[off+k]*DrTv1[k]+Jsy[off+k]*DsTv1[k]));

        dS13[off+k] -= wi[i]*wi[j]*mu[off+k]*JI[off+k]
          *(Jrx[off+k]*DrTv3[k] + Jsx[off+k]*DsTv3[k]);

        dS23[off+k] -= wi[i]*wi[j]*mu[off+k]*JI[off+k]
          *(Jry[off+k]*DrTv3[k] + Jsy[off+k]*DsTv3[k]);
      }

    /* loop over faces */
    for(bfam_locidx_t face = 0; face < 4; face++)
    {
      for(bfam_locidx_t pnt = 0; pnt < Nfp; pnt++)
      {
        bfam_locidx_t f = pnt + Nfp*(face + 4*e);
        bfam_locidx_t iM = vmapM[f];
        bfam_locidx_t iP = vmapP[f];

        /* Setup stuff for the minus side */
        bfam_real_t Zsm = Zs[iM];
        bfam_real_t Zpm = Zp[iM];

        bfam_real_t nm[] = {n1[f],n2[f],0};

        bfam_real_t Tpm[] = {
          nm[0]*S11[iM]+nm[1]*S12[iM],
          nm[0]*S12[iM]+nm[1]*S22[iM],
          nm[0]*S13[iM]+nm[1]*S23[iM],
        };
        bfam_real_t Tnm = Tpm[0]*nm[0]+Tpm[1]*nm[1];
        Tpm[0] = Tpm[0]-Tnm*nm[0];
        Tpm[1] = Tpm[1]-Tnm*nm[1];

        bfam_real_t vpm[] = {v1[iM],v2[iM],v3[iM]};
        bfam_real_t vnm = nm[0]*vpm[0]+nm[1]*vpm[1];
        vpm[0] = vpm[0]-vnm*nm[0];
        vpm[1] = vpm[1]-vnm*nm[1];

        /* Setup stuff for the plus side */
        bfam_real_t Zsp = Zs[iP];
        bfam_real_t Zpp = Zp[iP];

        bfam_real_t np[] = {-nm[0],-nm[1],-nm[2]};

        bfam_real_t Tpp[] = {
          np[0]*S11[iP]+np[1]*S12[iP],
          np[0]*S12[iP]+np[1]*S22[iP],
          np[0]*S13[iP]+np[1]*S23[iP],
        };
        bfam_real_t Tnp = Tpp[0]*np[0]+Tpp[1]*np[1];
        Tpp[0] = Tpp[0]-Tnp*np[0];
        Tpp[1] = Tpp[1]-Tnp*np[1];

        bfam_real_t vpp[] = {v1[iP],v2[iP],v3[iP]};
        bfam_real_t vnp = np[0]*vpp[0]+np[1]*vpp[1];
        vpp[0] = vpp[0]-vnp*np[0];
        vpp[1] = vpp[1]-vnp*np[1];

        bfam_real_t TnS;
        bfam_real_t TpS[3];
        bfam_real_t vnS;
        bfam_real_t vpS[3];

        beard_dgx_upwind_state_m(&TnS,TpS,&vnS,vpS,
            Tnm, Tnp, Tpm, Tpp, vnm, vnp, vpm, vpp, Zpm, Zpp, Zsm, Zsp);

        TnS    -= Tnm;
        TpS[0] -= Tpm[0];
        TpS[1] -= Tpm[1];
        TpS[2] -= Tpm[2];

        /* intra */
        beard_dgx_add_flux(1, TnS,TpS,vnS,vpS,iM,
            dv1,dv2,dv3, dS11,dS22,dS33,dS12,dS13,dS23,
            lam[iM],mu[iM],rhoi[iM],nm,sJ[f],JI[iM],wi[0]);
      }
    }
  }
}

void BFAM_APPEND_EXPAND(beard_dgx_scale_rates_elastic_,NORDER)(
    int inN, bfam_subdomain_dgx_quad_t *sub, const char *rate_prefix,
    const bfam_long_real_t a)
{
#ifdef USE_GENERIC
  const int Np  = (inN+1)*(inN+1);
  BFAM_WARNING("Using generic scale rates function");
#endif
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

void BFAM_APPEND_EXPAND(beard_dgx_scale_rates_slip_weakening_,NORDER)(
    int inN, bfam_subdomain_dgx_quad_glue_t *sub, const char *rate_prefix,
    const bfam_long_real_t a)
{
#ifdef USE_GENERIC
  const int N   = inN;
  BFAM_WARNING("Using generic scale rates function");
#endif
  const bfam_locidx_t num_pts = sub->K * (N+1);
  bfam_dictionary_t *fields = &sub->base.fields;
  const char *f_names[] = {"Dp","Dn",NULL};
  for(bfam_locidx_t f = 0; f_names[f]!=NULL; ++f)
  {
    BFAM_LOAD_FIELD_RESTRICT_ALIGNED(rate,rate_prefix,f_names[f],fields);
    for(bfam_locidx_t n = 0; n < num_pts; n++) rate[n] *= a;
  }
}


void BFAM_APPEND_EXPAND(beard_dgx_add_rates_elastic_,NORDER)(
    int inN, bfam_subdomain_dgx_quad_t *sub, const char *field_prefix_lhs,
    const char *field_prefix_rhs, const char *rate_prefix,
    const bfam_long_real_t a)
{
#ifdef USE_GENERIC
  const int Np  = (inN+1)*(inN+1);
  BFAM_WARNING("Using generic add rates function");
#endif

  const char *f_names[] = {"v1", "v2", "v3",
    "S11", "S22", "S33", "S12", "S13", "S23", NULL};

  /* get the fields we will need */
  bfam_dictionary_t *fields = &sub->base.fields;

  for(bfam_locidx_t f = 0; f_names[f] != NULL; f++)
  {
    BFAM_LOAD_FIELD_ALIGNED(         lhs ,field_prefix_lhs,f_names[f],fields);
    BFAM_LOAD_FIELD_ALIGNED(         rhs ,field_prefix_rhs,f_names[f],fields);
    BFAM_LOAD_FIELD_RESTRICT_ALIGNED(rate,rate_prefix     ,f_names[f],fields);
    for(bfam_locidx_t p = 0; p < sub->K*Np; p++) lhs[p] = rhs[p] + a*rate[p];
  }
}

void BFAM_APPEND_EXPAND(beard_dgx_add_rates_slip_weakening_,NORDER)(
    int inN, bfam_subdomain_dgx_quad_glue_t *sub, const char *field_prefix_lhs,
    const char *field_prefix_rhs, const char *rate_prefix,
    const bfam_long_real_t a)
{
#ifdef USE_GENERIC
  const int N = inN;
  BFAM_WARNING("Using generic add rates function");
#endif

  const char *f_names[] = {"Dp","Dn",NULL};

  /* get the fields we will need */
  bfam_dictionary_t *fields = &sub->base.fields;

  for(bfam_locidx_t f = 0; f_names[f] != NULL; f++)
  {
    BFAM_LOAD_FIELD_ALIGNED(         lhs ,field_prefix_lhs,f_names[f],fields);
    BFAM_LOAD_FIELD_ALIGNED(         rhs ,field_prefix_rhs,f_names[f],fields);
    BFAM_LOAD_FIELD_RESTRICT_ALIGNED(rate,rate_prefix     ,f_names[f],fields);
    for(bfam_locidx_t p = 0; p < sub->K*(N+1); p++) lhs[p] = rhs[p] + a*rate[p];
  }
}

void BFAM_APPEND_EXPAND(beard_dgx_inter_rhs_boundary_,NORDER)(
    int inN, bfam_subdomain_dgx_quad_glue_t *sub_g, const char *rate_prefix,
    const char *field_prefix, const bfam_long_real_t t, const bfam_real_t R)
{
#ifdef USE_GENERIC
  /*const int N   = inN;*/
  /*const int Np  = (inN+1)*(inN+1);*/
  const int Nfp = inN+1;
  BFAM_WARNING("Using generic inter rhs function");
#endif

  bfam_subdomain_dgx_quad_t* sub_m = sub_g->sub_m;

  /* get the fields we will need */
  bfam_dictionary_t *fields = &sub_m->base.fields;
  bfam_dictionary_t *fields_face = &sub_m->base.fields_face;

  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(v1 ,field_prefix,"v1" ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(v2 ,field_prefix,"v2" ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(v3 ,field_prefix,"v3" ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S11,field_prefix,"S11",fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S22,field_prefix,"S22",fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S33,field_prefix,"S33",fields);
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
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(rhoi,"","rho_inv"  ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(lam ,"","lam"      ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(mu  ,"","mu"       ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Zs  ,"","Zs"       ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Zp  ,"","Zp"       ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(J   ,"","_grid_J"  ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(JI  ,"","_grid_JI" ,fields);

  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(n1  ,"","_grid_nx",fields_face);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(n2  ,"","_grid_ny",fields_face);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(sJ  ,"","_grid_sJ",fields_face);

  bfam_real_t *wi  = sub_m->wi;
  BFAM_ASSUME_ALIGNED(wi ,32);

  for(bfam_locidx_t le = 0; le < sub_g->K; le++)
  {
    bfam_locidx_t e = sub_g->EToEm[le];
    int8_t face = sub_g->EToFm[le];

    /* First remove the eroneous flux */
    beard_dgx_remove_flux(Nfp,face,e,sub_m->vmapM,n1,n2,Zs,Zp,
        mu,rhoi,lam,sJ,JI,wi,
        v1,v2,v3,S11,S22,S33,S12,S13,S23,
        dv1,dv2,dv3,dS11,dS22,dS33,dS12,dS13,dS23);

    for(bfam_locidx_t pnt = 0; pnt < Nfp; pnt++)
    {
      bfam_locidx_t f = pnt + Nfp*(face + 4*e);
      bfam_locidx_t iM = sub_m->vmapM[f];

      /* Setup stuff for the minus side */
      bfam_real_t Zsm = Zs[iM];
      bfam_real_t Zpm = Zp[iM];

      bfam_real_t nm[] = {n1[f],n2[f],0};

      bfam_real_t Tpm[] = {
        nm[0]*S11[iM]+nm[1]*S12[iM],
        nm[0]*S12[iM]+nm[1]*S22[iM],
        nm[0]*S13[iM]+nm[1]*S23[iM],
      };
      bfam_real_t Tnm = Tpm[0]*nm[0]+Tpm[1]*nm[1];
      Tpm[0] = Tpm[0]-Tnm*nm[0];
      Tpm[1] = Tpm[1]-Tnm*nm[1];

      bfam_real_t vpm[] = {v1[iM],v2[iM],v3[iM]};
      bfam_real_t vnm = nm[0]*vpm[0]+nm[1]*vpm[1];
      vpm[0] = vpm[0]-vnm*nm[0];
      vpm[1] = vpm[1]-vnm*nm[1];

      /* now add the real flux */
      /* Setup stuff for the plus side */
      bfam_real_t Zsp = Zsm;
      bfam_real_t Zpp = Zpm;

      bfam_real_t Tpp[] = {R*Tpm[0],R*Tpm[1],R*Tpm[2]};
      bfam_real_t Tnp   =-R*Tnm;

      bfam_real_t vpp[] = {R*vpm[0],R*vpm[1],R*vpm[2]};
      bfam_real_t vnp    =-R*vnm;

      bfam_real_t TnS;
      bfam_real_t TpS[3];
      bfam_real_t vnS;
      bfam_real_t vpS[3];

      beard_dgx_upwind_state_m(&TnS,TpS,&vnS,vpS,
          Tnm, Tnp, Tpm, Tpp, vnm, vnp, vpm, vpp, Zpm, Zpp, Zsm, Zsp);

      TnS -= Tnm;
      TpS[0] -= Tpm[0];
      TpS[1] -= Tpm[1];
      TpS[2] -= Tpm[2];

      beard_dgx_add_flux(1, TnS,TpS,vnS,vpS,iM,
          dv1,dv2,dv3, dS11,dS22,dS33,dS12,dS13,dS23,
          lam[iM],mu[iM],rhoi[iM],nm,sJ[f],JI[iM],wi[0]);
    }
  }
}

void BFAM_APPEND_EXPAND(beard_dgx_inter_rhs_interface_,NORDER)(
    int inN, bfam_subdomain_dgx_quad_glue_t *sub_g, const char *rate_prefix,
    const char *field_prefix, const bfam_long_real_t t)
{
#ifdef USE_GENERIC
  /*const int N   = inN;*/
  /*const int Np  = (inN+1)*(inN+1);*/
  const int Nfp = inN+1;
  BFAM_WARNING("Using generic inter rhs function");
#endif

  bfam_subdomain_dgx_quad_t* sub_m = sub_g->sub_m;

  /* get the fields we will need */
  bfam_dictionary_t *fields_m    = &sub_g->base.fields_m;
  bfam_dictionary_t *fields_p    = &sub_g->base.fields_p;
  bfam_dictionary_t *fields      = &sub_m->base.fields;
  bfam_dictionary_t *fields_face = &sub_m->base.fields_face;

  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(v1 ,field_prefix,"v1" ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(v2 ,field_prefix,"v2" ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(v3 ,field_prefix,"v3" ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S11,field_prefix,"S11",fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S22,field_prefix,"S22",fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S33,field_prefix,"S33",fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S12,field_prefix,"S12",fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S13,field_prefix,"S13",fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S23,field_prefix,"S23",fields);

  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(v1_m ,field_prefix,"v1" ,fields_m);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(v2_m ,field_prefix,"v2" ,fields_m);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(v3_m ,field_prefix,"v3" ,fields_m);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S11_m,field_prefix,"S11",fields_m);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S22_m,field_prefix,"S22",fields_m);
  /*
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S33_m,field_prefix,"S33",fields_m);
  */
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S12_m,field_prefix,"S12",fields_m);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S13_m,field_prefix,"S13",fields_m);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S23_m,field_prefix,"S23",fields_m);

  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(v1_p ,field_prefix,"v1" ,fields_p);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(v2_p ,field_prefix,"v2" ,fields_p);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(v3_p ,field_prefix,"v3" ,fields_p);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S11_p,field_prefix,"S11",fields_p);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S22_p,field_prefix,"S22",fields_p);
  /*
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S33_p,field_prefix,"S33",fields_p);
  */
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S12_p,field_prefix,"S12",fields_p);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S13_p,field_prefix,"S13",fields_p);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S23_p,field_prefix,"S23",fields_p);

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
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(rhoi,"","rho_inv"  ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(lam ,"","lam"      ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(mu  ,"","mu"       ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Zs  ,"","Zs"       ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Zp  ,"","Zp"       ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(J   ,"","_grid_J"  ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(JI  ,"","_grid_JI" ,fields);

  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(n1  ,"","_grid_nx",fields_face);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(n2  ,"","_grid_ny",fields_face);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(sJ  ,"","_grid_sJ",fields_face);

  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Zs_m  ,"","Zs"       ,fields_m);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Zp_m  ,"","Zp"       ,fields_m);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Zs_p  ,"","Zs"       ,fields_p);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Zp_p  ,"","Zp"       ,fields_p);

  bfam_real_t *wi  = sub_m->wi;
  BFAM_ASSUME_ALIGNED(wi ,32);

  for(bfam_locidx_t le = 0; le < sub_g->K; le++)
  {
    bfam_locidx_t e = sub_g->EToEm[le];
    int8_t face = sub_g->EToFm[le];

    /* Assumes conforming straight sided elements */
    bfam_real_t nm[] = {n1[Nfp*(face+4*e)],n2[Nfp*(face+4*e)],0};

    if(sub_g->EToHm[le] < 2)
      beard_dgx_remove_flux(Nfp,face,e,sub_m->vmapM,n1,n2,Zs,Zp,
          mu,rhoi,lam,sJ,JI,wi,
          v1,v2,v3,S11,S22,S33,S12,S13,S23,
          dv1,dv2,dv3,dS11,dS22,dS33,dS12,dS13,dS23);


#ifndef USE_GENERIC
#undef Np
#endif
    bfam_locidx_t Np_g = sub_g->Np;
#ifndef USE_GENERIC
#define Np (N+1)*(N+1)
#endif

    bfam_real_t TpS_g[3*Np_g];
    bfam_real_t TnS_g[  Np_g];
    bfam_real_t vpS_g[3*Np_g];
    bfam_real_t vnS_g[  Np_g];

    for(bfam_locidx_t pnt = 0; pnt < Np_g; pnt++)
    {
      bfam_locidx_t iG = pnt+le*Np_g;

      /* Setup stuff for the minus side */
      bfam_real_t Zsm = Zs_m[iG];
      bfam_real_t Zpm = Zp_m[iG];


      bfam_real_t Tpm[] = {
        nm[0]*S11_m[iG]+nm[1]*S12_m[iG],
        nm[0]*S12_m[iG]+nm[1]*S22_m[iG],
        nm[0]*S13_m[iG]+nm[1]*S23_m[iG],
      };
      bfam_real_t Tnm = Tpm[0]*nm[0]+Tpm[1]*nm[1];
      Tpm[0] = Tpm[0]-Tnm*nm[0];
      Tpm[1] = Tpm[1]-Tnm*nm[1];

      bfam_real_t vpm[] = {v1_m[iG],v2_m[iG],v3_m[iG]};
      bfam_real_t vnm = nm[0]*vpm[0]+nm[1]*vpm[1];
      vpm[0] = vpm[0]-vnm*nm[0];
      vpm[1] = vpm[1]-vnm*nm[1];

      /* now add the real flux */
      /* Setup stuff for the plus side */
      bfam_real_t Zsp = Zs_p[iG];
      bfam_real_t Zpp = Zp_p[iG];

      bfam_real_t np[] = {-nm[0],-nm[1],0};

      bfam_real_t Tpp[] = {
        np[0]*S11_p[iG]+np[1]*S12_p[iG],
        np[0]*S12_p[iG]+np[1]*S22_p[iG],
        np[0]*S13_p[iG]+np[1]*S23_p[iG],
      };
      bfam_real_t Tnp = Tpp[0]*np[0]+Tpp[1]*np[1];
      Tpp[0] = Tpp[0]-Tnp*np[0];
      Tpp[1] = Tpp[1]-Tnp*np[1];

      bfam_real_t vpp[] = {v1_p[iG], v2_p[iG], v3_p[iG]};
      bfam_real_t vnp = np[0]*vpp[0]+np[1]*vpp[1];
      vpp[0] = vpp[0]-vnp*np[0];
      vpp[1] = vpp[1]-vnp*np[1];

      beard_dgx_upwind_state_m(
          &TnS_g[pnt],&TpS_g[3*pnt],&vnS_g[pnt],&vpS_g[3*pnt],
          Tnm, Tnp, Tpm, Tpp, vnm, vnp, vpm, vpp, Zpm, Zpp, Zsm, Zsp);

      /* substract off the grid values */
      TpS_g[3*pnt+0] -= Tpm[0];
      TpS_g[3*pnt+1] -= Tpm[1];
      TpS_g[3*pnt+2] -= Tpm[2];
      TnS_g[pnt]     -= Tnm;
    }

    bfam_real_t *restrict TpS_m;
    bfam_real_t *restrict TnS_m;
    bfam_real_t *restrict vpS_m;
    bfam_real_t *restrict vnS_m;

    /* these will be used to the store the projected values if we need them */
    BFAM_ALIGN(32) bfam_real_t TpS_m_STORAGE[3*Nfp];
    BFAM_ALIGN(32) bfam_real_t TnS_m_STORAGE[  Nfp];
    BFAM_ALIGN(32) bfam_real_t vpS_m_STORAGE[3*Nfp];
    BFAM_ALIGN(32) bfam_real_t vnS_m_STORAGE[  Nfp];

    /* check to tee if projection */
    /* locked */
    if(0)
    {
      TpS_m = TpS_g;
      TnS_m = TnS_g;
      vpS_m = vpS_g;
      vnS_m = vnS_g;
    }
    else
    {
      /* set to the correct Mass times projection */
      TpS_m = TpS_m_STORAGE;
      TnS_m = TnS_m_STORAGE;
      vpS_m = vpS_m_STORAGE;
      vnS_m = vnS_m_STORAGE;
      project_flux(TnS_m, TpS_m, vnS_m, vpS_m, Nfp, Np_g, TnS_g, TpS_g, vnS_g,
          vpS_g, sub_g->massprojection[sub_g->EToHm[le]], wi);
    }

    for(bfam_locidx_t pnt = 0; pnt < Nfp; pnt++)
    {
      bfam_locidx_t f = pnt + Nfp*(face + 4*e);
      bfam_locidx_t iM = sub_m->vmapM[f];

      beard_dgx_add_flux(1,
          TnS_m[pnt],&TpS_m[3*pnt],vnS_m[pnt],&vpS_m[3*pnt],iM,
          dv1,dv2,dv3, dS11,dS22,dS33,dS12,dS13,dS23,
          lam[iM],mu[iM],rhoi[iM],nm,sJ[f],JI[iM],wi[0]);
    }
  }
}

void BFAM_APPEND_EXPAND(beard_dgx_inter_rhs_slip_weakening_interface_,NORDER)(
    int inN, bfam_subdomain_dgx_quad_glue_t *sub_g, const char *rate_prefix,
    const char *field_prefix, const bfam_long_real_t t)
{
#ifdef USE_GENERIC
  /*const int N   = inN;*/
  /*const int Np  = (inN+1)*(inN+1);*/
  const int Nfp = inN+1;
  BFAM_WARNING("Using generic inter rhs function");
#endif

  bfam_subdomain_dgx_quad_t* sub_m = sub_g->sub_m;

  /* get the fields we will need */
  bfam_dictionary_t *fields_m    = &sub_g->base.fields_m;
  bfam_dictionary_t *fields_p    = &sub_g->base.fields_p;
  bfam_dictionary_t *fields_g    = &sub_g->base.fields;
  bfam_dictionary_t *fields      = &sub_m->base.fields;
  bfam_dictionary_t *fields_face = &sub_m->base.fields_face;

  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(v1 ,field_prefix,"v1" ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(v2 ,field_prefix,"v2" ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(v3 ,field_prefix,"v3" ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S11,field_prefix,"S11",fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S22,field_prefix,"S22",fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S33,field_prefix,"S33",fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S12,field_prefix,"S12",fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S13,field_prefix,"S13",fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S23,field_prefix,"S23",fields);

  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(v1_m ,field_prefix,"v1" ,fields_m);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(v2_m ,field_prefix,"v2" ,fields_m);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(v3_m ,field_prefix,"v3" ,fields_m);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S11_m,field_prefix,"S11",fields_m);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S22_m,field_prefix,"S22",fields_m);
  /*
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S33_m,field_prefix,"S33",fields_m);
  */
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S12_m,field_prefix,"S12",fields_m);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S13_m,field_prefix,"S13",fields_m);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S23_m,field_prefix,"S23",fields_m);

  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(v1_p ,field_prefix,"v1" ,fields_p);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(v2_p ,field_prefix,"v2" ,fields_p);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(v3_p ,field_prefix,"v3" ,fields_p);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S11_p,field_prefix,"S11",fields_p);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S22_p,field_prefix,"S22",fields_p);
  /*
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S33_p,field_prefix,"S33",fields_p);
  */
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S12_p,field_prefix,"S12",fields_p);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S13_p,field_prefix,"S13",fields_p);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S23_p,field_prefix,"S23",fields_p);

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
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(rhoi,"","rho_inv"  ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(lam ,"","lam"      ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(mu  ,"","mu"       ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Zs  ,"","Zs"       ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Zp  ,"","Zp"       ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(J   ,"","_grid_J"  ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(JI  ,"","_grid_JI" ,fields);

  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(n1  ,"","_grid_nx",fields_face);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(n2  ,"","_grid_ny",fields_face);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(sJ  ,"","_grid_sJ",fields_face);

  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Zs_m  ,"","Zs"       ,fields_m);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Zp_m  ,"","Zp"       ,fields_m);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Zs_p  ,"","Zs"       ,fields_p);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Zp_p  ,"","Zp"       ,fields_p);

  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Tp1_0  ,"","Tp1_0" ,fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Tp2_0  ,"","Tp2_0" ,fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Tp3_0  ,"","Tp3_0" ,fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Tn_0   ,"","Tn_0"  ,fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Tp1    ,"","Tp1"   ,fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Tp2    ,"","Tp2"   ,fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Tp3    ,"","Tp3"   ,fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Tn     ,"","Tn"    ,fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(V      ,"","V"     ,fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Dc     ,"","Dc"    ,fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Dp     ,"","Dp"    ,fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(fs     ,"","fs"    ,fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(fd     ,"","fd"    ,fields_g);

  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(dDp,rate_prefix,"Dp",fields_g);

  bfam_real_t *wi  = sub_m->wi;
  BFAM_ASSUME_ALIGNED(wi ,32);

  for(bfam_locidx_t le = 0; le < sub_g->K; le++)
  {
    bfam_locidx_t e = sub_g->EToEm[le];
    int8_t face = sub_g->EToFm[le];

    /* Assumes conforming straight sided elements */
    bfam_real_t nm[] = {n1[Nfp*(face+4*e)],n2[Nfp*(face+4*e)],0};

    if(sub_g->EToHm[le] < 2)
      beard_dgx_remove_flux(Nfp,face,e,sub_m->vmapM,n1,n2,Zs,Zp,
          mu,rhoi,lam,sJ,JI,wi,
          v1,v2,v3,S11,S22,S33,S12,S13,S23,
          dv1,dv2,dv3,dS11,dS22,dS33,dS12,dS13,dS23);


#ifndef USE_GENERIC
#undef Np
#endif
    bfam_locidx_t Np_g = sub_g->Np;
#ifndef USE_GENERIC
#define Np (N+1)*(N+1)
#endif

    bfam_real_t TpS_g[3*Np_g];
    bfam_real_t TnS_g[  Np_g];
    bfam_real_t vpS_g[3*Np_g];
    bfam_real_t vnS_g[  Np_g];

    for(bfam_locidx_t pnt = 0; pnt < Np_g; pnt++)
    {
      bfam_locidx_t iG = pnt+le*Np_g;

      /* Setup stuff for the minus side */
      bfam_real_t Zsm = Zs_m[iG];
      bfam_real_t Zpm = Zp_m[iG];


      bfam_real_t Tpm[] = {
        nm[0]*S11_m[iG]+nm[1]*S12_m[iG],
        nm[0]*S12_m[iG]+nm[1]*S22_m[iG],
        nm[0]*S13_m[iG]+nm[1]*S23_m[iG],
      };
      bfam_real_t Tnm = Tpm[0]*nm[0]+Tpm[1]*nm[1];
      Tpm[0] = Tpm[0]-Tnm*nm[0];
      Tpm[1] = Tpm[1]-Tnm*nm[1];

      bfam_real_t vpm[] = {v1_m[iG],v2_m[iG],v3_m[iG]};
      bfam_real_t vnm = nm[0]*vpm[0]+nm[1]*vpm[1];
      vpm[0] = vpm[0]-vnm*nm[0];
      vpm[1] = vpm[1]-vnm*nm[1];

      /* now add the real flux */
      /* Setup stuff for the plus side */
      bfam_real_t Zsp = Zs_p[iG];
      bfam_real_t Zpp = Zp_p[iG];

      bfam_real_t np[] = {-nm[0],-nm[1],0};

      bfam_real_t Tpp[] = {
        np[0]*S11_p[iG]+np[1]*S12_p[iG],
        np[0]*S12_p[iG]+np[1]*S22_p[iG],
        np[0]*S13_p[iG]+np[1]*S23_p[iG],
      };
      bfam_real_t Tnp = Tpp[0]*np[0]+Tpp[1]*np[1];
      Tpp[0] = Tpp[0]-Tnp*np[0];
      Tpp[1] = Tpp[1]-Tnp*np[1];

      bfam_real_t vpp[] = {v1_p[iG], v2_p[iG], v3_p[iG]};
      bfam_real_t vnp = np[0]*vpp[0]+np[1]*vpp[1];
      vpp[0] = vpp[0]-vnp*np[0];
      vpp[1] = vpp[1]-vnp*np[1];

      /* compute the upwind state assuming that the fault is locked */
      beard_dgx_upwind_state_m(
          &TnS_g[pnt],&TpS_g[3*pnt],&vnS_g[pnt],&vpS_g[3*pnt],
          Tnm, Tnp, Tpm, Tpp, vnm, vnp, vpm, vpp, Zpm, Zpp, Zsm, Zsp);

      Tn[iG] = TnS_g[pnt]+Tn_0[iG];
      BFAM_ABORT_IF(Tn[iG] > 0, "fault opening not implemented");

      bfam_real_t Slock2 =
        + (TpS_g[3*pnt+0]+Tp1_0[iG])*(TpS_g[3*pnt+0]+Tp1_0[iG])
        + (TpS_g[3*pnt+1]+Tp2_0[iG])*(TpS_g[3*pnt+1]+Tp2_0[iG])
        + (TpS_g[3*pnt+2]+Tp3_0[iG])*(TpS_g[3*pnt+2]+Tp3_0[iG]);

      bfam_real_t Sfric =
        -Tn[iG]*(fs[iG]-(fs[iG]-fd[iG])*BFAM_MIN(Dp[iG],Dc[iG])/Dc[iG]);

      V[iG] = 0;
      if(Sfric*Sfric < Slock2)
      {
        bfam_real_t Vps[3];
        const bfam_real_t Tp0[] = {Tp1_0[iG],Tp2_0[iG],Tp3_0[iG]};
        beard_dgx_upwind_state_friction_m(&TpS_g[3*pnt], &vpS_g[3*pnt], Vps,
            Sfric, Tpm, Tpp, Tp0, vpm, vpp, Zsm, Zsp);
        V[iG] = BFAM_REAL_SQRT(Vps[0]*Vps[0] + Vps[1]*Vps[1] + Vps[2]*Vps[2]);
      }
      dDp[iG] += V[iG];
      Tp1[iG] = TpS_g[3*pnt+0];
      Tp2[iG] = TpS_g[3*pnt+1];
      Tp3[iG] = TpS_g[3*pnt+2];

      /* substract off the grid values */
      TpS_g[3*pnt+0] -= Tpm[0];
      TpS_g[3*pnt+1] -= Tpm[1];
      TpS_g[3*pnt+2] -= Tpm[2];
      TnS_g[pnt]     -= Tnm;
    }

    bfam_real_t *restrict TpS_m;
    bfam_real_t *restrict TnS_m;
    bfam_real_t *restrict vpS_m;
    bfam_real_t *restrict vnS_m;

    /* these will be used to the store the projected values if we need them */
    BFAM_ALIGN(32) bfam_real_t TpS_m_STORAGE[3*Nfp];
    BFAM_ALIGN(32) bfam_real_t TnS_m_STORAGE[  Nfp];
    BFAM_ALIGN(32) bfam_real_t vpS_m_STORAGE[3*Nfp];
    BFAM_ALIGN(32) bfam_real_t vnS_m_STORAGE[  Nfp];

    /* check to tee if projection */
    /* slip weakening */
    if(0)
    {
      TpS_m = TpS_g;
      TnS_m = TnS_g;
      vpS_m = vpS_g;
      vnS_m = vnS_g;
    }
    else
    {
      /* set to the correct Mass times projection */
      TpS_m = TpS_m_STORAGE;
      TnS_m = TnS_m_STORAGE;
      vpS_m = vpS_m_STORAGE;
      vnS_m = vnS_m_STORAGE;
      project_flux(TnS_m, TpS_m, vnS_m, vpS_m, Nfp, Np_g, TnS_g, TpS_g, vnS_g,
          vpS_g, sub_g->massprojection[sub_g->EToHm[le]], wi);
    }

    for(bfam_locidx_t pnt = 0; pnt < Nfp; pnt++)
    {
      bfam_locidx_t f = pnt + Nfp*(face + 4*e);
      bfam_locidx_t iM = sub_m->vmapM[f];

      beard_dgx_add_flux(1,
          TnS_m[pnt],&TpS_m[3*pnt],vnS_m[pnt],&vpS_m[3*pnt],iM,
          dv1,dv2,dv3, dS11,dS22,dS33,dS12,dS13,dS23,
          lam[iM],mu[iM],rhoi[iM],nm,sJ[f],JI[iM],wi[0]);
    }
  }
}

void BFAM_APPEND_EXPAND(beard_dgx_energy_,NORDER)(
    int inN, bfam_real_t *energy_sq,
    bfam_subdomain_dgx_quad_t *sub, const char *field_prefix)
{
#ifdef USE_GENERIC
  const int N   = inN;
  const int Np  = (inN+1)*(inN+1);
  BFAM_WARNING("Using generic compute energy function");
#endif

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
    bfam_locidx_t off = e*Np;
    for(bfam_locidx_t i = 0; i < N+1; i++)
      for(bfam_locidx_t j = 0; j < N+1; j++)
      {
        /* node id */
        bfam_locidx_t nid = i*(N+1)+j+off;

        /* setup up deviatoric stress tensor */
        bfam_real_t mean_stress = (S11[nid]+S22[nid]+S33[nid])/3.0;
        bfam_real_t s11 = S11[nid] - mean_stress;
        bfam_real_t s22 = S22[nid] - mean_stress;
        bfam_real_t s33 = S33[nid] - mean_stress;
        bfam_real_t s12 = S12[nid];
        bfam_real_t s13 = S13[nid];
        bfam_real_t s23 = S23[nid];

        /* bulk modulus */
        bfam_real_t K = lam[nid] + 2.0*mu[nid]/3.0;

        energy_sq[0] += w[i]*w[j]*J[nid]*(
            rho[nid]*(v1[nid]*v1[nid] + v2[nid]*v2[nid] + v3[nid]*v3[nid])/2
            + (s11*s11 + s22*s22 + s33*s33
              + 2*s12*s12 + 2*s13*s13 + 2*s23*s23)/(4*mu[nid])
            + mean_stress*mean_stress/(2*K)
          );
      }
  }
}
