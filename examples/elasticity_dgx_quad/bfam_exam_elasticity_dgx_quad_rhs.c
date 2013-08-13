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
  for(int i = 0; i < 3; i++)
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

  char tmp_name[BFAM_BUFSIZ];

  /* get the fields we will need */
  bfam_dictionary_t *fields = &sub->base.fields;
  bfam_dictionary_t *fields_face = &sub->base.fields_face;

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

  /*
  snprintf(tmp_name,BFAM_BUFSIZ,"%sS33",field_prefix);
  bfam_real_t *restrict S33 = bfam_dictionary_get_value_ptr(fields, tmp_name);
  */

  BFAM_ASSUME_ALIGNED(v1,32);
  BFAM_ASSUME_ALIGNED(v2,32);
  BFAM_ASSUME_ALIGNED(v3,32);
  BFAM_ASSUME_ALIGNED(S11,32);
  BFAM_ASSUME_ALIGNED(S12,32);
  BFAM_ASSUME_ALIGNED(S13,32);
  BFAM_ASSUME_ALIGNED(S22,32);
  BFAM_ASSUME_ALIGNED(S23,32);
  /* BFAM_ASSUME_ALIGNED(S33,32); */

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

  BFAM_ASSUME_ALIGNED(dv1,32);
  BFAM_ASSUME_ALIGNED(dv2,32);
  BFAM_ASSUME_ALIGNED(dv3,32);
  BFAM_ASSUME_ALIGNED(dS11,32);
  BFAM_ASSUME_ALIGNED(dS12,32);
  BFAM_ASSUME_ALIGNED(dS13,32);
  BFAM_ASSUME_ALIGNED(dS22,32);
  BFAM_ASSUME_ALIGNED(dS23,32);
  BFAM_ASSUME_ALIGNED(dS33,32);

  /* get the material properties and metric terms */
  bfam_real_t *restrict rhoi = bfam_dictionary_get_value_ptr(fields,"rho_inv");
  bfam_real_t *restrict lam = bfam_dictionary_get_value_ptr(fields,"lam");
  bfam_real_t *restrict mu  = bfam_dictionary_get_value_ptr(fields,"mu" );
  bfam_real_t *restrict Zs  = bfam_dictionary_get_value_ptr(fields,"Zs" );
  bfam_real_t *restrict Zp  = bfam_dictionary_get_value_ptr(fields,"Zp" );
  bfam_real_t *restrict J   = bfam_dictionary_get_value_ptr(fields,"_grid_J");
  bfam_real_t *restrict JI  = bfam_dictionary_get_value_ptr(fields,"_grid_JI");
  bfam_real_t *restrict Jrx = bfam_dictionary_get_value_ptr(fields,"_grid_Jrx");
  bfam_real_t *restrict Jry = bfam_dictionary_get_value_ptr(fields,"_grid_Jry");
  bfam_real_t *restrict Jsx = bfam_dictionary_get_value_ptr(fields,"_grid_Jsx");
  bfam_real_t *restrict Jsy = bfam_dictionary_get_value_ptr(fields,"_grid_Jsy");
  BFAM_ASSUME_ALIGNED(rhoi,32);
  BFAM_ASSUME_ALIGNED(lam,32);
  BFAM_ASSUME_ALIGNED(mu ,32);
  BFAM_ASSUME_ALIGNED(J  ,32);
  BFAM_ASSUME_ALIGNED(Jrx,32);
  BFAM_ASSUME_ALIGNED(Jry,32);
  BFAM_ASSUME_ALIGNED(Jsx,32);
  BFAM_ASSUME_ALIGNED(Jsy,32);

  bfam_real_t *restrict n1 =
    bfam_dictionary_get_value_ptr(fields_face,"_grid_nx");
  bfam_real_t *restrict n2 =
    bfam_dictionary_get_value_ptr(fields_face,"_grid_ny");
  bfam_real_t *restrict sJ =
    bfam_dictionary_get_value_ptr(fields_face,"_grid_sJ");
  BFAM_ASSUME_ALIGNED(n1,32);
  BFAM_ASSUME_ALIGNED(n2,32);
  BFAM_ASSUME_ALIGNED(sJ,32);

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
    for(int i = 0; i < Np; i++)
    {
      aux1[i] = Jrx[off+i] * S11[off+i] + Jry[off+i] * S12[off+i];
      aux2[i] = Jsx[off+i] * S11[off+i] + Jsy[off+i] * S12[off+i];
    }
    BFAM_KRON_IXA   (N+1, Dr     , aux1   , rate); /* Dr */
    BFAM_KRON_AXI_PE(N+1, Dr     , aux2   , rate); /* Ds */

    for(int i = 0; i < Np; i++)
      dv1[off+i] += JI[off+i]*rhoi[off+i]*rate[i];

    /* v2 += JI*rhoi*(Dr*(Jrx*S12+Jry*S22)) + JI*rhoi*(Ds*(Jsx*S12+Jsy*S22));*/
    for(int i = 0; i < Np; i++)
    {
      aux1[i] = Jrx[off+i] * S12[off+i] + Jry[off+i] * S22[off+i];
      aux2[i] = Jsx[off+i] * S12[off+i] + Jsy[off+i] * S22[off+i];
    }
    BFAM_KRON_IXA   (N+1, Dr     , aux1   , rate); /* Dr */
    BFAM_KRON_AXI_PE(N+1, Dr     , aux2   , rate); /* Ds */

    for(int i = 0; i < Np; i++)
      dv2[off+i] += JI[off+i]*rhoi[off+i]*rate[i];

    /* rho*v3 += JI*rhoi*(Dr*(Jrx*S13+Jry*S23)) +  JI*rhoi*(Ds*(Jsx*S13+Jsy*S23));*/
    for(int i = 0; i < Np; i++)
    {
      aux1[i] = Jrx[off+i] * S13[off+i] + Jry[off+i] * S23[off+i];
      aux2[i] = Jsx[off+i] * S13[off+i] + Jsy[off+i] * S23[off+i];
    }
    BFAM_KRON_IXA   (N+1, Dr     , aux1   , rate); /* Dr */
    BFAM_KRON_AXI_PE(N+1, Dr     , aux2   , rate); /* Ds */

    for(int i = 0; i < Np; i++)
      dv3[off+i] += JI[off+i]*rhoi[off+i]*rate[i];

    /* we use these a lot so store them */
    for(int i = 0; i < N+1; i++)
      for(int j = 0; j < N+1; j++)
      {
        int k = i*(N+1)+j;
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
    for(int i = 0; i < N+1; i++)
      for(int j = 0; j < N+1; j++)
      {
        int k = i*(N+1)+j;

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
  char name[BFAM_BUFSIZ];
  for(int f = 0; f_names[f]!=NULL; ++f)
  {
    snprintf(name,BFAM_BUFSIZ,"%s%s",rate_prefix,f_names[f]);
    bfam_real_t *restrict rate = bfam_dictionary_get_value_ptr(fields, name);
    BFAM_ASSUME_ALIGNED(rate,32);
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

  const char *field_names[] = {"v1", "v2", "v3",
    "S11", "S22", "S33", "S12", "S13", "S23", NULL};
  char tmp_name[BFAM_BUFSIZ];

  /* get the fields we will need */
  bfam_dictionary_t *fields = &sub->base.fields;

  for(int f = 0; field_names[f] != NULL; f++)
  {
    snprintf(tmp_name,BFAM_BUFSIZ,"%s%s",field_prefix_lhs,field_names[f]);
    bfam_real_t *lhs = bfam_dictionary_get_value_ptr(fields, tmp_name);
    BFAM_ASSUME_ALIGNED(lhs,32);

    snprintf(tmp_name,BFAM_BUFSIZ,"%s%s",field_prefix_rhs,field_names[f]);
    bfam_real_t *rhs = bfam_dictionary_get_value_ptr(fields, tmp_name);
    BFAM_ASSUME_ALIGNED(rhs,32);

    snprintf(tmp_name,BFAM_BUFSIZ,"%s%s",rate_prefix,field_names[f]);
    bfam_real_t *rate = bfam_dictionary_get_value_ptr(fields, tmp_name);
    BFAM_ASSUME_ALIGNED(rate,32);

    for(int p = 0; p < sub->K*Np; p++) lhs[p] = rhs[p] + a*rate[p];
  }
}

void BFAM_APPEND_EXPAND(bfam_elasticity_dgx_quad_inter_rhs_boundary_,NORDER)(
    int inN, bfam_subdomain_dgx_quad_glue_t *sub_f, const char *rate_prefix,
    const char *field_prefix, const bfam_long_real_t t, const bfam_real_t R)
{
#ifdef USE_GENERIC
  /*const int N   = inN;*/
  /*const int Np  = (inN+1)*(inN+1);*/
  const int Nfp = inN+1;
  BFAM_WARNING("Using generic inter rhs function");
#endif

  bfam_subdomain_dgx_quad_t* sub_m = sub_f->sub_m;

  char tmp_name[BFAM_BUFSIZ];

  /* get the fields we will need */
  bfam_dictionary_t *fields = &sub_m->base.fields;
  bfam_dictionary_t *fields_face = &sub_m->base.fields_face;

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

  /*
  snprintf(tmp_name,BFAM_BUFSIZ,"%sS33",field_prefix);
  bfam_real_t *restrict S33 = bfam_dictionary_get_value_ptr(fields, tmp_name);
  */

  BFAM_ASSUME_ALIGNED(v1,32);
  BFAM_ASSUME_ALIGNED(v2,32);
  BFAM_ASSUME_ALIGNED(v3,32);
  BFAM_ASSUME_ALIGNED(S11,32);
  BFAM_ASSUME_ALIGNED(S12,32);
  BFAM_ASSUME_ALIGNED(S13,32);
  BFAM_ASSUME_ALIGNED(S22,32);
  BFAM_ASSUME_ALIGNED(S23,32);
  /* BFAM_ASSUME_ALIGNED(S33,32); */

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

  BFAM_ASSUME_ALIGNED(dv1,32);
  BFAM_ASSUME_ALIGNED(dv2,32);
  BFAM_ASSUME_ALIGNED(dv3,32);
  BFAM_ASSUME_ALIGNED(dS11,32);
  BFAM_ASSUME_ALIGNED(dS12,32);
  BFAM_ASSUME_ALIGNED(dS13,32);
  BFAM_ASSUME_ALIGNED(dS22,32);
  BFAM_ASSUME_ALIGNED(dS23,32);
  BFAM_ASSUME_ALIGNED(dS33,32);

  /* get the material properties and metric terms */
  bfam_real_t *restrict rhoi = bfam_dictionary_get_value_ptr(fields,"rho_inv");
  bfam_real_t *restrict lam = bfam_dictionary_get_value_ptr(fields,"lam");
  bfam_real_t *restrict mu  = bfam_dictionary_get_value_ptr(fields,"mu" );
  bfam_real_t *restrict Zs  = bfam_dictionary_get_value_ptr(fields,"Zs" );
  bfam_real_t *restrict Zp  = bfam_dictionary_get_value_ptr(fields,"Zp" );
  bfam_real_t *restrict J   = bfam_dictionary_get_value_ptr(fields,"_grid_J");
  bfam_real_t *restrict JI  = bfam_dictionary_get_value_ptr(fields,"_grid_JI");
  BFAM_ASSUME_ALIGNED(rhoi,32);
  BFAM_ASSUME_ALIGNED(lam,32);
  BFAM_ASSUME_ALIGNED(mu ,32);
  BFAM_ASSUME_ALIGNED(J  ,32);

  bfam_real_t *restrict n1 =
    bfam_dictionary_get_value_ptr(fields_face,"_grid_nx");
  bfam_real_t *restrict n2 =
    bfam_dictionary_get_value_ptr(fields_face,"_grid_ny");
  bfam_real_t *restrict sJ =
    bfam_dictionary_get_value_ptr(fields_face,"_grid_sJ");
  BFAM_ASSUME_ALIGNED(n1,32);
  BFAM_ASSUME_ALIGNED(n2,32);
  BFAM_ASSUME_ALIGNED(sJ,32);

  bfam_real_t *wi  = sub_m->wi;
  BFAM_ASSUME_ALIGNED(wi ,32);

  for(bfam_locidx_t le = 0; le < sub_f->K; le++)
  {
    bfam_locidx_t e = sub_f->EToEm[le];
    int8_t face = sub_f->EToFm[le];
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
