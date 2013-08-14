#include "bfam_exam_elasticity_dgx_quad_rhs.h"

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
}                                                                              \
BFAM_ASSUME_ALIGNED(field,32);
#define BFAM_LOAD_FIELD_ALIGNED(field,prefix,base,dictionary)         \
bfam_real_t *field;                                                   \
{                                                                              \
  char bfam_load_field_name[BFAM_BUFSIZ];                                      \
  snprintf(bfam_load_field_name,BFAM_BUFSIZ,"%s%s",(prefix),(base));           \
  field = bfam_dictionary_get_value_ptr(dictionary, bfam_load_field_name);     \
}                                                                              \
BFAM_ASSUME_ALIGNED(field,32);


void BFAM_APPEND_EXPAND(bfam_elasticity_dgx_quad_print_order_,NORDER)(int inN)
{

#ifdef USE_GENERIC
  const int N = inN;
  BFAM_WARNING("Using generic order print function");
#endif

  BFAM_INFO("N = %d", N);
}

static inline void
bfam_elasticity_dgx_quad_upwind_state_m(
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

  /* upwind perpendiculat velocitie and tractions */
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


void BFAM_APPEND_EXPAND(bfam_elasticity_dgx_quad_intra_rhs_elastic_,NORDER)(
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

  BFAM_ALIGN(32) bfam_real_t MJv1[Np];
  BFAM_ALIGN(32) bfam_real_t MJv2[Np];
  BFAM_ALIGN(32) bfam_real_t MJv3[Np];

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
        MJv1[k] = w[i]*w[j]*v1[off+k];
        MJv2[k] = w[i]*w[j]*v2[off+k];
        MJv3[k] = w[i]*w[j]*v3[off+k];
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
    BFAM_KRON_IXAT(N+1, Dr, MJv1, DrTv1);
    BFAM_KRON_ATXI(N+1, Dr, MJv1, DsTv1);
    BFAM_KRON_IXAT(N+1, Dr, MJv2, DrTv2);
    BFAM_KRON_ATXI(N+1, Dr, MJv2, DsTv2);
    BFAM_KRON_IXAT(N+1, Dr, MJv3, DrTv3);
    BFAM_KRON_ATXI(N+1, Dr, MJv3, DsTv3);
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

        bfam_elasticity_dgx_quad_upwind_state_m(&TnS,TpS,&vnS,vpS,
            Tnm, Tnp, Tpm, Tpp, vnm, vnp, vpm, vpp, Zpm, Zpp, Zsm, Zsp);

        /* velocities for the flux */
        bfam_real_t vS[] = {vpS[0] + nm[0]*vnS,vpS[1] + nm[1]*vnS,vpS[2]};

        /* add the flux back in */
        bfam_real_t JI_wi_sJ = wi[0] * sJ[f] * JI[iM];

        bfam_real_t  lamM =  lam[iM];
        bfam_real_t   muM =   mu[iM];
        bfam_real_t rhoiM = rhoi[iM];

        dv1[iM] += rhoiM*JI_wi_sJ*(TpS[0]-Tpm[0]+(TnS-Tnm)*nm[0]);
        dv2[iM] += rhoiM*JI_wi_sJ*(TpS[1]-Tpm[1]+(TnS-Tnm)*nm[1]);
        dv3[iM] += rhoiM*JI_wi_sJ*(TpS[2]-Tpm[2]);

        dS11[iM] += JI_wi_sJ*((lamM+2*muM)*vS[0]*nm[0] +  lamM       *vS[1]*nm[1]);
        dS22[iM] += JI_wi_sJ*( lamM       *vS[0]*nm[0] + (lamM+2*muM)*vS[1]*nm[1]);
        dS33[iM] += JI_wi_sJ*( lamM       *vS[0]*nm[0] +  lamM       *vS[1]*nm[1]);
        dS12[iM] += JI_wi_sJ*muM*(vS[0]*nm[1] + vS[1]*nm[0]);
        dS13[iM] += JI_wi_sJ*muM *vS[2]*nm[0];
        dS23[iM] += JI_wi_sJ*muM *vS[2]*nm[1];
      }
    }
  }
}

void BFAM_APPEND_EXPAND(bfam_elasticity_dgx_quad_scale_rates_elastic_,NORDER)(
    int inN, bfam_subdomain_dgx_quad_t *sub, const char *rate_prefix,
    const bfam_long_real_t a)
{
#ifdef USE_GENERIC
  const int Np  = (inN+1)*(inN+1);
  BFAM_WARNING("Using generic intra rhs function");
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

void BFAM_APPEND_EXPAND(bfam_elasticity_dgx_quad_add_rates_elastic_,NORDER)(
    int inN, bfam_subdomain_dgx_quad_t *sub, const char *field_prefix_lhs,
    const char *field_prefix_rhs, const char *rate_prefix,
    const bfam_long_real_t a)
{
#ifdef USE_GENERIC
  const int Np  = (inN+1)*(inN+1);
  BFAM_WARNING("Using generic intra rhs function");
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

void BFAM_APPEND_EXPAND(bfam_elasticity_dgx_quad_inter_rhs_boundary_,NORDER)(
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

  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(n1  ,"","_grid_nx",fields_face);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(n2  ,"","_grid_ny",fields_face);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(sJ  ,"","_grid_sJ",fields_face);

  bfam_real_t *wi  = sub_m->wi;
  BFAM_ASSUME_ALIGNED(wi ,32);

  for(bfam_locidx_t le = 0; le < sub_g->K; le++)
  {
    bfam_locidx_t e = sub_g->EToEm[le];
    int8_t face = sub_g->EToFm[le];
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

      /* First remove what we already did */
      // bfam_real_t Zsp = Zs[iM];
      // bfam_real_t Zpp = Zp[iM];
      bfam_real_t Zsp = Zsm;
      bfam_real_t Zpp = Zpm;

      /*
      bfam_real_t np[] = {-nm[0],-nm[1],-nm[2]};

      bfam_real_t Tpp[] = {
        np[0]*S11[iM]+np[1]*S12[iM],
        np[0]*S12[iM]+np[1]*S22[iM],
        np[0]*S13[iM]+np[1]*S23[iM],
      };
      bfam_real_t Tnp = Tpp[0]*np[0]+Tpp[1]*np[1];
      Tpp[0] = Tpp[0]-Tnp*np[0];
      Tpp[1] = Tpp[1]-Tnp*np[1];

      bfam_real_t vpp[] = {v1[iM],v2[iM],v3[iM]};
      bfam_real_t vnp = np[0]*vpp[0]+np[1]*vpp[1];
      vpp[0] = vpp[0]-vnp*np[0];
      vpp[1] = vpp[1]-vnp*np[1];
      */


      bfam_real_t Tpp[] = {-Tpm[0],-Tpm[1],-Tpm[2]};
      bfam_real_t Tnp   = Tnm;

      bfam_real_t vpp[] = { vpm[0], vpm[1], vpm[2]};
      bfam_real_t vnp   = -vnm;

      bfam_real_t TnS;
      bfam_real_t TpS[3];
      bfam_real_t vnS;
      bfam_real_t vpS[3];

      bfam_elasticity_dgx_quad_upwind_state_m(&TnS,TpS,&vnS,vpS,
          Tnm, Tnp, Tpm, Tpp, vnm, vnp, vpm, vpp, Zpm, Zpp, Zsm, Zsp);

      /* velocities for the flux */
      bfam_real_t vS[] = {vpS[0] + nm[0]*vnS,vpS[1] + nm[1]*vnS,vpS[2]};

      /* add the flux back in */
      bfam_real_t JI_wi_sJ = wi[0] * sJ[f] * JI[iM];

      bfam_real_t  lamM =  lam[iM];
      bfam_real_t   muM =   mu[iM];
      bfam_real_t rhoiM = rhoi[iM];

      dv1[iM]  -= rhoiM*JI_wi_sJ*(TpS[0]-Tpm[0]+(TnS-Tnm)*nm[0]);
      dv2[iM]  -= rhoiM*JI_wi_sJ*(TpS[1]-Tpm[1]+(TnS-Tnm)*nm[1]);
      dv3[iM]  -= rhoiM*JI_wi_sJ*(TpS[2]-Tpm[2]);

      dS11[iM] -= JI_wi_sJ*((lamM+2*muM)*vS[0]*nm[0] +  lamM       *vS[1]*nm[1]);
      dS22[iM] -= JI_wi_sJ*( lamM       *vS[0]*nm[0] + (lamM+2*muM)*vS[1]*nm[1]);
      dS33[iM] -= JI_wi_sJ*( lamM       *vS[0]*nm[0] +  lamM       *vS[1]*nm[1]);
      dS12[iM] -= JI_wi_sJ*muM*(vS[0]*nm[1] + vS[1]*nm[0]);
      dS13[iM] -= JI_wi_sJ*muM *vS[2]*nm[0];
      dS23[iM] -= JI_wi_sJ*muM *vS[2]*nm[1];

      /* now add the real flux */
      /* Setup stuff for the plus side */
      Tpp[0] = R*Tpm[0];
      Tpp[1] = R*Tpm[1];
      Tpp[2] = R*Tpm[2];
      Tnp    =-R*Tnm;

      vpp[0] = R*vpm[0];
      vpp[1] = R*vpm[1];
      vpp[2] = R*vpm[2];
      vnp    =-R*vnm;

      bfam_elasticity_dgx_quad_upwind_state_m(&TnS,TpS,&vnS,vpS,
          Tnm, Tnp, Tpm, Tpp, vnm, vnp, vpm, vpp, Zpm, Zpp, Zsm, Zsp);

      vS[0] = vpS[0] + nm[0]*vnS;
      vS[1] = vpS[1] + nm[1]*vnS;
      vS[2] = vpS[2];

      dv1[iM] += rhoiM*JI_wi_sJ*(TpS[0]-Tpm[0]+(TnS-Tnm)*nm[0]);
      dv2[iM] += rhoiM*JI_wi_sJ*(TpS[1]-Tpm[1]+(TnS-Tnm)*nm[1]);
      dv3[iM] += rhoiM*JI_wi_sJ*(TpS[2]-Tpm[2]);

      dS11[iM] += JI_wi_sJ*((lamM+2*muM)*vS[0]*nm[0] +  lamM       *vS[1]*nm[1]);
      dS22[iM] += JI_wi_sJ*( lamM       *vS[0]*nm[0] + (lamM+2*muM)*vS[1]*nm[1]);
      dS33[iM] += JI_wi_sJ*( lamM       *vS[0]*nm[0] +  lamM       *vS[1]*nm[1]);
      dS12[iM] += JI_wi_sJ*muM*(vS[0]*nm[1] + vS[1]*nm[0]);
      dS13[iM] += JI_wi_sJ*muM *vS[2]*nm[0];
      dS23[iM] += JI_wi_sJ*muM *vS[2]*nm[1];
    }
  }
}

void BFAM_APPEND_EXPAND(bfam_elasticity_dgx_quad_inter_rhs_interface_,NORDER)(
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

  bfam_real_t *wi  = sub_m->wi;
  BFAM_ASSUME_ALIGNED(wi ,32);

  for(bfam_locidx_t le = 0; le < sub_g->K; le++)
  {
    bfam_locidx_t e = sub_g->EToEm[le];
    int8_t face = sub_g->EToFm[le];
    for(bfam_locidx_t pnt = 0; pnt < Nfp; pnt++)
    {
      bfam_locidx_t f = pnt + Nfp*(face + 4*e);
      bfam_locidx_t iM = sub_m->vmapM[f];
      bfam_locidx_t iG = pnt+le*Nfp;

      /* Setup stuff for the minus side */
      bfam_real_t Zsm = Zs[iM];
      bfam_real_t Zpm = Zp[iM];

      bfam_real_t nm[] = {n1[f],n2[f],0};

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

      /* First remove what we already did */
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

      bfam_elasticity_dgx_quad_upwind_state_m(&TnS,TpS,&vnS,vpS,
          Tnm, Tnp, Tpm, Tpp, vnm, vnp, vpm, vpp, Zpm, Zpp, Zsm, Zsp);

      /* velocities for the flux */
      bfam_real_t vS[] = {vpS[0] + nm[0]*vnS,vpS[1] + nm[1]*vnS,vpS[2]};

      /* add the flux back in */
      bfam_real_t JI_wi_sJ = wi[0] * sJ[f] * JI[iM];

      bfam_real_t  lamM =  lam[iM];
      bfam_real_t   muM =   mu[iM];
      bfam_real_t rhoiM = rhoi[iM];

      dv1[iM]  -= rhoiM*JI_wi_sJ*(TpS[0]-Tpm[0]+(TnS-Tnm)*nm[0]);
      dv2[iM]  -= rhoiM*JI_wi_sJ*(TpS[1]-Tpm[1]+(TnS-Tnm)*nm[1]);
      dv3[iM]  -= rhoiM*JI_wi_sJ*(TpS[2]-Tpm[2]);

      dS11[iM] -= JI_wi_sJ*((lamM+2*muM)*vS[0]*nm[0] +  lamM       *vS[1]*nm[1]);
      dS22[iM] -= JI_wi_sJ*( lamM       *vS[0]*nm[0] + (lamM+2*muM)*vS[1]*nm[1]);
      dS33[iM] -= JI_wi_sJ*( lamM       *vS[0]*nm[0] +  lamM       *vS[1]*nm[1]);
      dS12[iM] -= JI_wi_sJ*muM*(vS[0]*nm[1] + vS[1]*nm[0]);
      dS13[iM] -= JI_wi_sJ*muM *vS[2]*nm[0];
      dS23[iM] -= JI_wi_sJ*muM *vS[2]*nm[1];

      /* now add the real flux */
      /* Setup stuff for the plus side */
      Zsp = Zs[iM];
      Zpp = Zp[iM];

      bfam_real_t np[] = {-nm[0],-nm[1],0};

      Tpp[0] = np[0]*S11_p[iG]+np[1]*S12_p[iG];
      Tpp[1] = np[0]*S12_p[iG]+np[1]*S22_p[iG];
      Tpp[2] = np[0]*S13_p[iG]+np[1]*S23_p[iG];
      Tnp = Tpp[0]*np[0]+Tpp[1]*np[1];
      Tpp[0] = Tpp[0]-Tnp*np[0];
      Tpp[1] = Tpp[1]-Tnp*np[1];

      vpp[0] = v1_p[iG];
      vpp[1] = v2_p[iG];
      vpp[2] = v3_p[iG];
      vnp = np[0]*vpp[0]+np[1]*vpp[1];
      vpp[0] = vpp[0]-vnp*np[0];
      vpp[1] = vpp[1]-vnp*np[1];

      bfam_elasticity_dgx_quad_upwind_state_m(&TnS,TpS,&vnS,vpS,
          Tnm, Tnp, Tpm, Tpp, vnm, vnp, vpm, vpp, Zpm, Zpp, Zsm, Zsp);

      vS[0] = vpS[0] + nm[0]*vnS;
      vS[1] = vpS[1] + nm[1]*vnS;
      vS[2] = vpS[2];

      dv1[iM] += rhoiM*JI_wi_sJ*(TpS[0]-Tpm[0]+(TnS-Tnm)*nm[0]);
      dv2[iM] += rhoiM*JI_wi_sJ*(TpS[1]-Tpm[1]+(TnS-Tnm)*nm[1]);
      dv3[iM] += rhoiM*JI_wi_sJ*(TpS[2]-Tpm[2]);

      dS11[iM] += JI_wi_sJ*((lamM+2*muM)*vS[0]*nm[0] +  lamM       *vS[1]*nm[1]);
      dS22[iM] += JI_wi_sJ*( lamM       *vS[0]*nm[0] + (lamM+2*muM)*vS[1]*nm[1]);
      dS33[iM] += JI_wi_sJ*( lamM       *vS[0]*nm[0] +  lamM       *vS[1]*nm[1]);
      dS12[iM] += JI_wi_sJ*muM*(vS[0]*nm[1] + vS[1]*nm[0]);
      dS13[iM] += JI_wi_sJ*muM *vS[2]*nm[0];
      dS23[iM] += JI_wi_sJ*muM *vS[2]*nm[1];
    }
  }
}

void BFAM_APPEND_EXPAND(bfam_elasticity_dgx_quad_energy_,NORDER)(
    int inN, bfam_real_t *energy_sq, 
    bfam_subdomain_dgx_quad_t *sub, const char *field_prefix)
{
#ifdef USE_GENERIC
  const int N   = inN;
  const int Np  = (inN+1)*(inN+1);
  BFAM_WARNING("Using generic inter rhs function");
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
