#ifdef BEARD_DGX_DIMENSION

/* START DIM 2 */
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

#define BEARD_APPEND_4(a,b,c,d) a ## b ## c ##d
#define BEARD_APPEND_EXPAND_4(a,b,c,d) BEARD_APPEND_4(a,b,c,d)
#define beard_dgx_intra_rhs_elastic \
  BEARD_APPEND_EXPAND_4(beard_dgx_intra_rhs_elastic_,DIM,_,NORDER)

#define beard_dgx_scale_rates_elastic \
  BEARD_APPEND_EXPAND_4(beard_dgx_scale_rates_elastic_,DIM,_,NORDER)

#define beard_dgx_scale_rates_interface \
  BEARD_APPEND_EXPAND_4(beard_dgx_scale_rates_interface_,DIM,_,NORDER)

#define beard_dgx_scale_rates \
  BEARD_APPEND_EXPAND_4(beard_dgx_scale_rates_,DIM,_,NORDER)

#define beard_dgx_add_rates_elastic \
  BEARD_APPEND_EXPAND_4(beard_dgx_add_rates_elastic_,DIM,_,NORDER)

#define beard_dgx_add_rates_interface \
  BEARD_APPEND_EXPAND_4(beard_dgx_add_rates_interface_,DIM,_,NORDER)

#define beard_dgx_add_rates_glue_p \
  BEARD_APPEND_EXPAND_4(beard_dgx_add_rates_glue_p_,DIM,_,NORDER)

#define beard_dgx_add_rates \
  BEARD_APPEND_EXPAND_4(beard_dgx_add_rates_,DIM,_,NORDER)

#define beard_dgx_inter_rhs_boundary \
  BEARD_APPEND_EXPAND_4(beard_dgx_inter_rhs_boundary_,DIM,_,NORDER)

#define beard_dgx_inter_rhs_interface \
  BEARD_APPEND_EXPAND_4(beard_dgx_inter_rhs_interface_,DIM,_,NORDER)

#define beard_dgx_inter_rhs_slip_weakening_interface \
  BEARD_APPEND_EXPAND_4(beard_dgx_inter_rhs_slip_weakening_interface_, \
                        DIM,_,NORDER)

#define beard_dgx_inter_rhs_ageing_law_interface \
  BEARD_APPEND_EXPAND_4(beard_dgx_inter_rhs_ageing_law_interface_, \
                        DIM,_,NORDER)

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


/* Computational choice macros */
#define MAX_ITER         (100)
#define FTOL             (BFAM_REAL_EPS*100)
#define ATOL             (0)
#define RTOL             (0)
#define BEARD_STATE      beard_dgx_upwind_state_m
#define MASSPROJECTION   (1)
#define PROJECTION       (0)
#define NO_OPENING       (0 == 1)

static inline void
beard_dgx_upwind_state_friction_m(
          bfam_real_t *TpS,       bfam_real_t *vpSm,        bfam_real_t *VpS,
    const bfam_real_t    T,
    const bfam_real_t *Tpm, const bfam_real_t *Tpp  , const bfam_real_t *Tp0,
    const bfam_real_t *vpm, const bfam_real_t *vpp  ,
    const bfam_real_t  Zsm, const bfam_real_t  Zsp  )
{
  /*
   * upwind perpendicular velocities and tractions
   *
   * LHS are the volume values and RHS are the upwind states and we find them by
   * forcing these characteristics to remain unchanged which satisfying the
   * constraint that
   *
   *    T = | TpS[i] |
   *
   * and that TpS points in the same direciton as VpS
   *
   * The following characteristics must remain unchanges
   *
   *   wpm = Tpm - Zsm*vpm = TpS - Zsm*vpSm
   *   wpp = Tpp - Zsp*vpp =-TpS - Zsp*vpSp
   *
   * scaling and adding gives the relationship
   *
   *   phi[i] = (wpm[i]*Zsp-wpp[i]*Zsm)/(Zsp+Zsm)
   *          = TpS[i] + eta*VpS[i]
   *
   *   Vps[i] = vpSp[i] - vpSm[i]
   *
   *   eta    = ((Zsm*Zsp)/(Zsm+Zsp))
   *
   *
   * which since TpS and VpS must be parallel gives the directions of TpS and
   * VpS. Thus we have that TpS[i] = A phi[i] / |phi|, but since the magnitude
   * of TpS[i] is T, we have that
   *
   *    TpS[i] = T phi[i]/|phi|
   *
   * Note: that this TpS[i] includes the load / background stress, so that must
   * be subtracted off to get the actual upwind state
   */
  bfam_real_t phi[3];
  bfam_real_t phi_m = 0;
  const bfam_real_t eta = ((Zsp*Zsm)/(Zsm+Zsp));
  for(bfam_locidx_t i = 0; i < 3; i++)
  {
    bfam_real_t wpm = Tpm[i] + Tp0[i] - Zsm*vpm[i];
    bfam_real_t wpp = Tpp[i] - Tp0[i] - Zsp*vpp[i];

    /* phi[i] = TpS[i] + eta*VpS[i] */
    phi[i] = (wpm*Zsp-wpp*Zsm)/(Zsp+Zsm);

    /* phi_m = | phi[i] | */
    phi_m += phi[i]*phi[i];
  }

  phi_m = BFAM_REAL_SQRT(phi_m);
  for(bfam_locidx_t i = 0; i < 3; i++)
  {
    TpS[i]  = T*phi[i]/phi_m - Tp0[i];
    vpSm[i] = (TpS[i]-Tpm[i])/Zsm + vpm[i];
    VpS[i]  = (phi[i] - TpS[i] - Tp0[i])/eta;
  }
}

static inline bfam_locidx_t
beard_dgx_rs_arcsinh_newton(bfam_real_t *V_ptr, bfam_real_t *T_ptr,
    const bfam_real_t eta, const bfam_real_t phi,
    const bfam_real_t a  , const bfam_real_t psi,
    const bfam_real_t Tn , const bfam_real_t V0,
    const bfam_locidx_t max_iterations,
    const bfam_real_t ftol,
    const bfam_real_t atol,
    const bfam_real_t rtol)
{

  /*
   * We solve the relationship
   *
   *    phi - eta V = -Tn f(V; psi, a, V0) = tau
   *    f(V; psi, a, V0) = a arcsinh((V/(2 V0)) exp(psi/a))
   *
   * with a bracketed Newton's method. The initial brackets are:
   *
   *    V = [0, phi/eta]
   *
   * since this corresponds to tau = phi and tau = 0 (respectively)
   */

  bfam_real_t V = *V_ptr;
  bfam_real_t Vmin = 0;
  bfam_real_t Vmax = phi/eta;

  if( V < Vmin )
    V = Vmin;
  else if(V > Vmax)
    V = Vmax;

  const bfam_real_t A  = BFAM_REAL_EXP(psi/a)/(2*V0);

  /* Compute the initial function value */
  bfam_real_t G = phi - eta*V + Tn*a*BFAM_REAL_ASINH(A*V);
  if(G < 0)
    Vmax = V;
  else if(G > 0)
    Vmin = V;
  else
  {
    *V_ptr = V;
    *T_ptr = -Tn*a*BFAM_REAL_ASINH(A*V);
    return 0;
  }

  for(bfam_locidx_t iter = 0; iter<max_iterations; iter++)
  {
    const bfam_real_t dG = - eta + Tn*a*A/BFAM_REAL_SQRT(A*A*V*V+1);
    bfam_real_t dV = -G/dG;
    V += dV;

    /* If outside the bounds use the midpoint */
    if( V < Vmin || V > Vmax)
    {
      V  = 0.5*(Vmin + Vmax);
      dV =      Vmax - Vmin ;
    }

    /* Update the objective function */
    const bfam_real_t G2 = -Tn*a*BFAM_REAL_ASINH(A*V);
    G = phi - eta*V - G2;

    /* Stop if converged */
    if(  BFAM_REAL_ABS(G) < ftol*(1+G2)    ||
        (BFAM_REAL_ABS(dV)<=atol+rtol*(BFAM_REAL_ABS(V)+BFAM_REAL_ABS(dV))))
    {
      *V_ptr = V;
      *T_ptr = -Tn*a*BFAM_REAL_ASINH(A*V);
      return 0;
    }

    /* Update the bounds function */
    if(G < 0)
      Vmax = V;
    else if(G > 0)
      Vmin = V;
  }

  BFAM_INFO("A:     %"BFAM_REAL_FMTe,A);
  BFAM_INFO("phi:   %"BFAM_REAL_FMTe,phi);
  BFAM_INFO("psi:   %"BFAM_REAL_FMTe,psi);
  BFAM_INFO("psi/a: %"BFAM_REAL_FMTe,psi/a);

  G = phi - eta*V + Tn*a*BFAM_REAL_ASINH(A*Vmin);
  BFAM_INFO("Gmin:  %"BFAM_REAL_FMTe,G);
  BFAM_INFO("Vmin:  %"BFAM_REAL_FMTe,Vmin);

  G = phi - eta*V + Tn*a*BFAM_REAL_ASINH(A*V);
  BFAM_INFO("G:     %"BFAM_REAL_FMTe,G);
  BFAM_INFO("V:     %"BFAM_REAL_FMTe,Vmin);

  G = phi - eta*V + Tn*a*BFAM_REAL_ASINH(A*Vmax);
  BFAM_INFO("Gmax:  %"BFAM_REAL_FMTe,G);
  BFAM_INFO("Vmax:  %"BFAM_REAL_FMTe,Vmax);

  return 1;
}

static inline void
beard_dgx_upwind_state_rate_and_state_friction_m(
          bfam_real_t *TpS,       bfam_real_t *vpSm ,       bfam_real_t *VpS,
    const bfam_real_t   Tn, const bfam_real_t  a    ,
    const bfam_real_t   V0, const bfam_real_t  psi  ,
    const bfam_real_t *Tpm, const bfam_real_t *Tpp  , const bfam_real_t *Tp0,
    const bfam_real_t *vpm, const bfam_real_t *vpp  ,
    const bfam_real_t  Zsm, const bfam_real_t  Zsp  )
{
  /*
   * Upwind perpendicular velocities and tractions
   *
   * To solve the rate and state friction law we need to find T and V that
   * satisfy
   *
   *    T = -Tn f(V, psi; a, V0)
   *
   * (note that with our definition Tn is negative in compression)  such that
   * the characteristics propagating into the interface remain unchanged, i.e.,
   * the following characteristics must remain unchanges
   *
   *   wpm = Tpm - Zsm*vpm = TpS - Zsm*vpSm
   *   wpp = Tpp - Zsp*vpp =-TpS - Zsp*vpSp
   *
   * Scaling and adding these relationships gives
   *
   *   phi[i] = (wpm[i]*Zsp-wpp[i]*Zsm)/(Zsp+Zsm)
   *          = TpS[i] + eta*VpS[i]
   *
   *   Vps[i] = vpSp[i] - vpSm[i]
   *
   *   eta    = ((Zsm*Zsp)/(Zsm+Zsp))
   *
   * which since TpS and VpS must be parallel gives the directions of TpS and
   * VpS. That is, we have that
   *
   *    |phi| u[i] = |TpS| u[i] + eta |VpS| u[i]
   *
   * which gives the radiation damping relationship
   *
   *    |phi| = |TpS| + eta |VpS|
   *
   * We thus must solve the following non-linear relationship:
   *
   *    |TpS] = |phi| - eta |VpS| = -Tn f(V; psi, a, V0)
   *
   * Note: that this TpS[i] includes the load / background stress, so that must
   * be subtracted off to get the actual upwind state
   */
  bfam_real_t phi[3];
  bfam_real_t phi_m = 0;
  const bfam_real_t eta = ((Zsp*Zsm)/(Zsm+Zsp));
  bfam_real_t V = 0;
  for(bfam_locidx_t i = 0; i < 3; i++)
  {
    bfam_real_t wpm = Tpm[i] + Tp0[i] - Zsm*vpm[i];
    bfam_real_t wpp = Tpp[i] - Tp0[i] - Zsp*vpp[i];

    /* phi[i] = TpS[i] + eta*VpS[i] */
    phi[i] = (wpm*Zsp-wpp*Zsm)/(Zsp+Zsm);

    /* phi_m = | phi[i] | */
    phi_m += phi[i]*phi[i];

    /* initialize V to the difference of nodal values */
    V = V + (vpp[i]-vpm[i])*(vpp[i]-vpm[i]);
  }
  phi_m = BFAM_REAL_SQRT(phi_m);
  V     = BFAM_REAL_SQRT(V);

  bfam_real_t T = 0;
  const bfam_locidx_t ret_val = beard_dgx_rs_arcsinh_newton(&V, &T, eta, phi_m,
                                                            a, psi, Tn, V0,
                                                            MAX_ITER, FTOL,
                                                            ATOL, RTOL);
  BFAM_ABORT_IF(ret_val,"newton solver returned %"BFAM_LOCIDX_PRId,ret_val);

  bfam_real_t tmp = 0;
  bfam_real_t vmp = 0;
  for(bfam_locidx_t i = 0; i < 3; i++)
  {
    TpS[i]  = T*phi[i]/phi_m - Tp0[i];
    vpSm[i] = (TpS[i]-Tpm[i])/Zsm + vpm[i];
    VpS[i]  = (phi[i] - (TpS[i] + Tp0[i]))/eta;
  }
}

static inline void
beard_dgx_upwind_state_m(
          bfam_real_t *Tns,       bfam_real_t *TpS,
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

  /* upwind perpendicular velocities and tractions */
  /* wpm = Tpm - Zsm*vpm = TpS - Zsm*vpS */
  /* wpp = Tpp - Zsp*vpp =-TpS - Zsp*vpS */
  for(bfam_locidx_t i = 0; i < 3; i++)
  {
    bfam_real_t wpm = Tpm[i] - Zsm*vpm[i];
    bfam_real_t wpp = Tpp[i] - Zsp*vpp[i];
    vps[i] =-(wpm    +wpp    )/(Zsp+Zsm);
    TpS[i] = (wpm*Zsp-wpp*Zsm)/(Zsp+Zsm);
  }
}

static inline void
beard_dgx_central_state_m(
          bfam_real_t *Tns,       bfam_real_t *TpS,
          bfam_real_t *vns,       bfam_real_t *vps,
    const bfam_real_t  Tnm, const bfam_real_t  Tnp,
    const bfam_real_t *Tpm, const bfam_real_t *Tpp,
    const bfam_real_t  vnm, const bfam_real_t  vnp,
    const bfam_real_t *vpm, const bfam_real_t *vpp,
    const bfam_real_t  Zpm, const bfam_real_t  Zpp,
    const bfam_real_t  Zsm, const bfam_real_t  Zsp)
{
  vns[0] = BFAM_REAL(0.5)*(vnm-vnp);
  Tns[0] = BFAM_REAL(0.5)*(Tnm+Tnp);

  for(bfam_locidx_t i = 0; i < 3; i++)
  {
    vps[i] = BFAM_REAL(0.5)*(vpm[i]+vpp[i]);
    TpS[i] = BFAM_REAL(0.5)*(Tpm[i]-Tpp[i]);
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

static inline void
beard_dgx_remove_flux( const int inN,
    bfam_locidx_t face, bfam_locidx_t e, const bfam_locidx_t *vmapM,
    const bfam_real_t *n1, const bfam_real_t *n2,
#if DIM==3
    const bfam_real_t *n3,
#endif
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
  GENERIC_INIT(inN,beard_dgx_remove_flux);

  for(bfam_locidx_t pnt = 0; pnt < Nfp; pnt++)
  {
    bfam_locidx_t f = pnt + Nfp*(face + Nfaces*e);
    bfam_locidx_t iM = vmapM[f];

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

    /* First remove what we already did */
    const bfam_real_t ZsP = ZsM;
    const bfam_real_t ZpP = ZpM;

    const bfam_real_t TpP[] = {-TpM[0],-TpM[1],-TpM[2]};
    const bfam_real_t TnP   = TnM;

    const bfam_real_t vpP[] = { vpM[0], vpM[1], vpM[2]};
    const bfam_real_t vnP   = -vnM;

    bfam_real_t TnS;
    bfam_real_t TpS[3];
    bfam_real_t vnS;
    bfam_real_t vpS[3];

    BEARD_STATE(&TnS,TpS,&vnS,vpS,
        TnM, TnP, TpM, TpP, vnM, vnP, vpM, vpP, ZpM, ZpP, ZsM, ZsP);

    TnS -= TnM;
    TpS[0] -= TpM[0];
    TpS[1] -= TpM[1];
    TpS[2] -= TpM[2];

    beard_dgx_add_flux(-1, TnS,TpS,vnS,vpS,iM,
        dv1,dv2,dv3, dS11,dS22,dS33,dS12,dS13,dS23,
        lam[iM],mu[iM],rhoi[iM],nM,sJ[f],JI[iM],wi[0]);
  }
}

static inline void
beard_dgx_add_boundary_flux( const int inN,
    bfam_locidx_t face, bfam_locidx_t e, const bfam_locidx_t *vmapM,
    const bfam_real_t *n1, const bfam_real_t *n2,
#if DIM==3
    const bfam_real_t *n3,
#endif
    const bfam_real_t *Zs, const bfam_real_t *Zp,
    const bfam_real_t *mu, const bfam_real_t *rhoi, const bfam_real_t *lam,
    const bfam_real_t *sJ, const bfam_real_t *JI,   const bfam_real_t *wi,
    const bfam_real_t *v1,  const bfam_real_t *v2,  const bfam_real_t *v3,
    const bfam_real_t *S11, const bfam_real_t *S22, const bfam_real_t *S33,
    const bfam_real_t *S12, const bfam_real_t *S13, const bfam_real_t *S23,
          bfam_real_t *dv1,  bfam_real_t *dv2,  bfam_real_t *dv3,
          bfam_real_t *dS11, bfam_real_t *dS22, bfam_real_t *dS33,
          bfam_real_t *dS12, bfam_real_t *dS13, bfam_real_t *dS23,
          bfam_real_t R
    )
{
  GENERIC_INIT(inN,beard_dgx_add_boundary_flux);

  for(bfam_locidx_t pnt = 0; pnt < Nfp; pnt++)
  {
    bfam_locidx_t f = pnt + Nfp*(face + Nfaces*e);
    bfam_locidx_t iM = vmapM[f];

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

    /* First remove what we already did */
    const bfam_real_t ZsP = ZsM;
    const bfam_real_t ZpP = ZpM;

    const bfam_real_t TpP[] = {R*TpM[0],R*TpM[1],R*TpM[2]};
    const bfam_real_t TnP   = -R*TnM;

    const bfam_real_t vpP[] = {R*vpM[0],R*vpM[1],R*vpM[2]};
    const bfam_real_t vnP   = -R*vnM;

    bfam_real_t TnS;
    bfam_real_t TpS[3];
    bfam_real_t vnS;
    bfam_real_t vpS[3];

    BEARD_STATE(&TnS,TpS,&vnS,vpS,
        TnM, TnP, TpM, TpP, vnM, vnP, vpM, vpP, ZpM, ZpP, ZsM, ZsP);

    TnS -= TnM;
    TpS[0] -= TpM[0];
    TpS[1] -= TpM[1];
    TpS[2] -= TpM[2];

    beard_dgx_add_flux(1, TnS,TpS,vnS,vpS,iM,
        dv1,dv2,dv3, dS11,dS22,dS33,dS12,dS13,dS23,
        lam[iM],mu[iM],rhoi[iM],nM,sJ[f],JI[iM],wi[0]);
  }
}

static inline void
beard_project_flux(bfam_real_t *Tns,       bfam_real_t *TpS,
                 bfam_real_t *vns,       bfam_real_t *vps,
               bfam_locidx_t  inN,     bfam_locidx_t  Nrpg,
           const bfam_real_t *Tng, const bfam_real_t *Tpg,
           const bfam_real_t *vng, const bfam_real_t *vpg,
           const bfam_real_t * P1
#if DIM==3
          ,const bfam_real_t * P2
#endif
           )
{
  GENERIC_INIT(inN,beard_project_flux);

  /* zero out the components */
  for(int i = 0; i < Nfp; i++)
  {
    for(int n = 0; n < 3; n++)
    {
      TpS[3*i+n] = 0;
      vps[3*i+n] = 0;
    }
    Tns[i] = 0;
    vns[i] = 0;
  }

  for(int m = 0, p = 0; m < Nrpg; m++)
#if DIM==3
    for(int l = 0; l < Nrpg; l++)
#endif
    {
      for(int j = 0, n = 0; j < N+1; j++)
#if DIM==3
        for(int i = 0; i < N+1; i++)
#endif
        {
          const bfam_real_t PR = BEARD_D3_AP(P1[j+m*(N+1)],
                                            *P2[i+l*(N+1)]);
          for(int k = 0; k < 3; k++)
          {
            TpS[3*n+k] += PR*Tpg[3*p+k];
            vps[3*n+k] += PR*vpg[3*p+k];
          }
          Tns[n] += PR*Tng[p];
          vns[n] += PR*vng[p];
          n++;
        }
      p++;
    }
}

static inline void
beard_massproject_flux(bfam_real_t *Tns,       bfam_real_t *TpS,
                 bfam_real_t *vns,       bfam_real_t *vps,
               bfam_locidx_t  inN,     bfam_locidx_t  Nrpg,
           const bfam_real_t *Tng, const bfam_real_t *Tpg,
           const bfam_real_t *vng, const bfam_real_t *vpg,
           const bfam_real_t * MP1,
#if DIM==3
           const bfam_real_t * MP2,
#endif
           const bfam_real_t  *wi)
{
  GENERIC_INIT(inN,beard_massproject_flux);

  /* zero out the components */
  for(int i = 0; i < Nfp; i++)
  {
    for(int n = 0; n < 3; n++)
    {
      TpS[3*i+n] = 0;
      vps[3*i+n] = 0;
    }
    Tns[i] = 0;
    vns[i] = 0;
  }

  for(int m = 0, p = 0; m < Nrpg; m++)
#if DIM==3
    for(int l = 0; l < Nrpg; l++)
#endif
    {
      for(int j = 0, n = 0; j < N+1; j++)
#if DIM==3
        for(int i = 0; i < N+1; i++)
#endif
        {
          const bfam_real_t wi_MP = BEARD_D3_AP(wi[j]*MP1[j+m*(N+1)],
                                               *wi[i]*MP2[i+l*(N+1)]);
          for(int k = 0; k < 3; k++)
          {
            TpS[3*n+k] += wi_MP*Tpg[3*p+k];
            vps[3*n+k] += wi_MP*vpg[3*p+k];
          }
          Tns[n] += wi_MP*Tng[p];
          vns[n] += wi_MP*vng[p];
          n++;
        }
      p++;
    }
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
          const bfam_real_t wi_ijk = BEARD_D3_AP(wi[i]*wi[j],*wi[k]);
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

void beard_dgx_scale_rates_interface(
    int inN, bfam_subdomain_dgx_t *sub, const char *rate_prefix,
    const bfam_long_real_t in_a, const char **f_names)
{
  GENERIC_INIT(inN,beard_dgx_scale_rates_interface);
  const bfam_locidx_t num_pts = sub->K * Nfp;
  bfam_dictionary_t *fields = &sub->base.fields;
  beard_dgx_scale_rates(inN, sub, rate_prefix, in_a, fields, f_names,
                        NULL, NULL, num_pts);
}

void beard_dgx_scale_rates_elastic(
    int inN, bfam_subdomain_dgx_t *sub, const char *rate_prefix,
    const bfam_long_real_t in_a)
{
  GENERIC_INIT(inN,beard_dgx_scale_rates_elastic);
  const bfam_locidx_t num_pts = sub->K * Np;
  bfam_dictionary_t *fields = &sub->base.fields;
  const char *f_names[] =
    {"v1","v2","v3","S11","S22","S33","S12","S13","S23",NULL};
  beard_dgx_scale_rates(inN, sub, rate_prefix, in_a, fields, f_names,
                        NULL, NULL, num_pts);

}

void beard_dgx_scale_rates(
    int inN, bfam_subdomain_dgx_t *sub, const char *rate_prefix,
    const bfam_long_real_t in_a, bfam_dictionary_t *fields,
    const char** scalars, const char** vectors, const char** tensors,
    const bfam_locidx_t num_pts)
{
  GENERIC_INIT(inN,beard_dgx_glue_p);
  const bfam_real_t a = (bfam_real_t) in_a;

  if(scalars)
  {
    for(int s = 0; scalars[s] != NULL;s++)
    {
      BFAM_LOAD_FIELD_RESTRICT_ALIGNED(rate,rate_prefix,scalars[s],fields);
      for(bfam_locidx_t n = 0; n < num_pts; n++) rate[n] *= a;
    }
  }

  if(vectors)
  {
    const char* postfix[] = {"n","p1","p2","p3",NULL};
    for(int s = 0; vectors[s] != NULL;s++)
    {
      char name[BFAM_BUFSIZ];
      for(int k = 0; postfix[k] != NULL;k++)
      {
        snprintf(name,BFAM_BUFSIZ,"%s%s",vectors[s],postfix[k]);
        BFAM_LOAD_FIELD_RESTRICT_ALIGNED(rate,rate_prefix,name,fields);
        for(bfam_locidx_t n = 0; n < num_pts; n++) rate[n] *= a;
      }
    }
  }

  if(tensors)
  {
    const char* postfix[] = {"n","p1","p2","p3",NULL};
    for(int s = 0; tensors[s] != NULL;s++)
    {
      char name[BFAM_BUFSIZ];
      for(int k = 0; postfix[k] != NULL;k++)
      {
        snprintf(name,BFAM_BUFSIZ,"%s%s",tensors[s],postfix[k]);
        BFAM_LOAD_FIELD_RESTRICT_ALIGNED(rate,rate_prefix,name,fields);
        for(bfam_locidx_t n = 0; n < num_pts; n++) rate[n] *= a;
      }
    }
  }
}


void beard_dgx_add_rates_interface(
    int inN, bfam_subdomain_dgx_t *sub, const char *field_prefix_lhs,
    const char *field_prefix_rhs, const char *rate_prefix,
    const bfam_long_real_t in_a, const char **f_names)
{
  GENERIC_INIT(inN,beard_dgx_add_rates_interface);
  bfam_dictionary_t *fields = &sub->base.fields;
  const bfam_locidx_t num_pts = sub->K * Nfp;
  beard_dgx_add_rates(inN, sub, field_prefix_lhs, field_prefix_rhs, rate_prefix,
                      in_a, fields, f_names, NULL, NULL, num_pts);
}

void beard_dgx_add_rates_elastic(
    int inN, bfam_subdomain_dgx_t *sub, const char *field_prefix_lhs,
    const char *field_prefix_rhs, const char *rate_prefix,
    const bfam_long_real_t in_a)
{
  GENERIC_INIT(inN,beard_dgx_add_rates_elastic);
  const char *f_names[] = {"v1", "v2", "v3",
    "S11", "S22", "S33", "S12", "S13", "S23", NULL};
  const bfam_locidx_t num_pts = sub->K * Np;
  bfam_dictionary_t *fields = &sub->base.fields;
  beard_dgx_add_rates(inN, sub, field_prefix_lhs, field_prefix_rhs, rate_prefix,
                      in_a, fields, f_names, NULL, NULL, num_pts);
}

void beard_dgx_add_rates_glue_p(
    int inN, bfam_subdomain_dgx_t *sub, const char *field_prefix_lhs,
    const char *field_prefix_rhs, const char *rate_prefix,
    const bfam_long_real_t in_a,
    const char** scalars, const char** vectors, const char** tensors)
{
  GENERIC_INIT(inN,beard_dgx_glue_p);
  const bfam_locidx_t num_pts = sub->K * Nfp;

  /* get the fields we will need */
  BFAM_ASSERT(sub->base.glue_p);
  bfam_dictionary_t *fields = &sub->base.glue_p->fields;

  beard_dgx_add_rates(inN, sub, field_prefix_lhs, field_prefix_rhs, rate_prefix,
                      in_a, fields, scalars, vectors, tensors, num_pts);
}

void beard_dgx_add_rates(
    int inN, bfam_subdomain_dgx_t *sub, const char *field_prefix_lhs,
    const char *field_prefix_rhs, const char *rate_prefix,
    const bfam_long_real_t in_a, bfam_dictionary_t *fields,
    const char** scalars, const char** vectors, const char** tensors,
    const bfam_locidx_t num_pts)
{
  GENERIC_INIT(inN,beard_dgx_glue_p);
  const bfam_real_t a = (bfam_real_t) in_a;

  if(scalars)
  {
    for(int s = 0; scalars[s] != NULL;s++)
    {
      BFAM_LOAD_FIELD_ALIGNED(         lhs ,field_prefix_lhs,scalars[s],fields);
      BFAM_LOAD_FIELD_ALIGNED(         rhs ,field_prefix_rhs,scalars[s],fields);
      BFAM_LOAD_FIELD_RESTRICT_ALIGNED(rate,rate_prefix     ,scalars[s],fields);
      for(bfam_locidx_t n = 0; n < num_pts; n++) lhs[n] = rhs[n] + a*rate[n];
    }
  }

  if(vectors)
  {
    const char* postfix[] = {"n","p1","p2","p3",NULL};
    for(int s = 0; vectors[s] != NULL;s++)
    {
      char name[BFAM_BUFSIZ];
      for(int k = 0; postfix[k] != NULL;k++)
      {
        snprintf(name,BFAM_BUFSIZ,"%s%s",vectors[s],postfix[k]);
        BFAM_LOAD_FIELD_ALIGNED(         lhs ,field_prefix_lhs,name,fields);
        BFAM_LOAD_FIELD_ALIGNED(         rhs ,field_prefix_rhs,name,fields);
        BFAM_LOAD_FIELD_RESTRICT_ALIGNED(rate,rate_prefix     ,name,fields);
        for(bfam_locidx_t n = 0; n < num_pts; n++) lhs[n] = rhs[n] + a*rate[n];
      }
    }
  }

  if(tensors)
  {
    const char* postfix[] = {"n","p1","p2","p3",NULL};
    for(int s = 0; tensors[s] != NULL;s++)
    {
      char name[BFAM_BUFSIZ];
      for(int k = 0; postfix[k] != NULL;k++)
      {
        snprintf(name,BFAM_BUFSIZ,"%s%s",tensors[s],postfix[k]);
        BFAM_LOAD_FIELD_ALIGNED(         lhs ,field_prefix_lhs,name,fields);
        BFAM_LOAD_FIELD_ALIGNED(         rhs ,field_prefix_rhs,name,fields);
        BFAM_LOAD_FIELD_RESTRICT_ALIGNED(rate,rate_prefix     ,name,fields);
        for(bfam_locidx_t n = 0; n < num_pts; n++) lhs[n] = rhs[n] + a*rate[n];
      }
    }
  }
}

void beard_dgx_inter_rhs_boundary(
    int inN, bfam_subdomain_dgx_t *sub_g, const char *rate_prefix,
    const char *field_prefix, const bfam_long_real_t t, const bfam_real_t R)
{
  GENERIC_INIT(inN,beard_dgx_inter_rhs_boundary);

  bfam_subdomain_dgx_t* sub_m =
    (bfam_subdomain_dgx_t*) sub_g->base.glue_m->sub_m;

  /* get the fields we will need */
  bfam_subdomain_dgx_glue_data_t* glue_p =
    (bfam_subdomain_dgx_glue_data_t*) sub_g->base.glue_p;
  BFAM_ASSERT(glue_p != NULL);
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

  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(n1,"","_grid_nx0",fields_face);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(n2,"","_grid_nx1",fields_face);
  BEARD_D3_OP(BFAM_LOAD_FIELD_RESTRICT_ALIGNED(n3,"","_grid_nx2",fields_face));
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(sJ  ,"","_grid_sJ",fields_face);

  bfam_real_t *wi  = sub_m->wi;
  BFAM_ASSUME_ALIGNED(wi ,32);

  BFAM_ASSERT(glue_p->EToEm);
  BFAM_ASSERT(glue_p->EToFm);

  for(bfam_locidx_t le = 0; le < sub_g->K; le++)
  {
    bfam_locidx_t e = glue_p->EToEm[le];
    int8_t face = glue_p->EToFm[le];


    beard_dgx_remove_flux(N,face,e,sub_m->vmapM,
        n1,n2,
#if DIM==3
        n3,
#endif
        Zs,Zp,
        mu,rhoi,lam,sJ,JI,wi,
        v1,v2,v3,S11,S22,S33,S12,S13,S23,
        dv1,dv2,dv3,dS11,dS22,dS33,dS12,dS13,dS23);

    beard_dgx_add_boundary_flux(N,face,e,sub_m->vmapM,
        n1,n2,
#if DIM==3
        n3,
#endif
        Zs,Zp,
        mu,rhoi,lam,sJ,JI,wi,
        v1,v2,v3,S11,S22,S33,S12,S13,S23,
        dv1,dv2,dv3,dS11,dS22,dS33,dS12,dS13,dS23,R);
  }

}

void beard_dgx_inter_rhs_interface(
    int inN, bfam_subdomain_dgx_t *sub_g, const char *rate_prefix,
    const char *field_prefix, const bfam_long_real_t t)
{
  GENERIC_INIT(inN,beard_dgx_inter_rhs_interface);

  bfam_subdomain_dgx_t* sub_m =
    (bfam_subdomain_dgx_t*) sub_g->base.glue_m->sub_m;

  /* get the fields we will need */
  bfam_subdomain_dgx_glue_data_t* glue_m =
    (bfam_subdomain_dgx_glue_data_t*) sub_g->base.glue_m;
  bfam_subdomain_dgx_glue_data_t* glue_p =
    (bfam_subdomain_dgx_glue_data_t*) sub_g->base.glue_p;
  BFAM_ASSERT(glue_m != NULL);
  BFAM_ASSERT(glue_p != NULL);
  bfam_dictionary_t *fields_m    = &glue_m->base.fields;
  bfam_dictionary_t *fields_p    = &glue_p->base.fields;
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

  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(vn_M ,field_prefix,"vn" ,fields_m);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(vp1_M,field_prefix,"vp1",fields_m);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(vp2_M,field_prefix,"vp2",fields_m);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(vp3_M,field_prefix,"vp3",fields_m);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Tn_M ,field_prefix,"Tn" ,fields_m);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Tp1_M,field_prefix,"Tp1",fields_m);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Tp2_M,field_prefix,"Tp2",fields_m);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Tp3_M,field_prefix,"Tp3",fields_m);

  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(vn_P ,field_prefix,"vn" ,fields_p);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(vp1_P,field_prefix,"vp1",fields_p);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(vp2_P,field_prefix,"vp2",fields_p);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(vp3_P,field_prefix,"vp3",fields_p);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Tn_P ,field_prefix,"Tn" ,fields_p);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Tp1_P,field_prefix,"Tp1",fields_p);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Tp2_P,field_prefix,"Tp2",fields_p);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Tp3_P,field_prefix,"Tp3",fields_p);

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

  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(n1,"","_grid_nx0",fields_face);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(n2,"","_grid_nx1",fields_face);
  BEARD_D3_OP(BFAM_LOAD_FIELD_RESTRICT_ALIGNED(n3,"","_grid_nx2",fields_face));
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(sJ  ,"","_grid_sJ",fields_face);

  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Zs_M  ,"","Zs"       ,fields_m);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Zp_M  ,"","Zp"       ,fields_m);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Zs_P  ,"","Zs"       ,fields_p);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Zp_P  ,"","Zp"       ,fields_p);

  bfam_real_t *wi  = sub_m->wi;
  BFAM_ASSUME_ALIGNED(wi ,32);

  BFAM_ASSERT(glue_p->EToEm);
  BFAM_ASSERT(glue_p->EToFm);
  BFAM_ASSERT(glue_p->EToHm);

  for(bfam_locidx_t le = 0; le < sub_g->K; le++)
  {
    bfam_locidx_t e = glue_p->EToEm[le];
    int8_t face = glue_p->EToFm[le];


    if(glue_p->EToHm[le] < 2)
      beard_dgx_remove_flux(N,face,e,sub_m->vmapM,
          n1,n2,
#if DIM==3
          n3,
#endif
          Zs,Zp,
          mu,rhoi,lam,sJ,JI,wi,
          v1,v2,v3,S11,S22,S33,S12,S13,S23,
          dv1,dv2,dv3,dS11,dS22,dS33,dS12,dS13,dS23);


#ifndef USE_GENERIC
#undef Np
#undef N
#endif
    bfam_locidx_t Np_g  = sub_g->Np;
    bfam_locidx_t Nrp_g = sub_g->N+1;
#ifndef USE_GENERIC
#define Np Np_BACK
#define N  NORDER
#endif

    bfam_real_t TpS_g[3*Np_g];
    bfam_real_t TnS_g[  Np_g];
    bfam_real_t vpS_g[3*Np_g];
    bfam_real_t vnS_g[  Np_g];

    for(bfam_locidx_t pnt = 0; pnt < Np_g; pnt++)
    {
      bfam_locidx_t iG = pnt+le*Np_g;

      /* Setup stuff for the minus side */
      bfam_real_t ZsM = Zs_M[iG];
      bfam_real_t ZpM = Zp_M[iG];

      bfam_real_t TpM[] = {Tp1_M[iG], Tp2_M[iG], Tp3_M[iG]};
      bfam_real_t TnM = Tn_M[iG];

      bfam_real_t vpM[] = {vp1_M[iG],vp2_M[iG],vp3_M[iG]};
      bfam_real_t vnM = vn_M[iG];

      /* now add the real flux */
      /* Setup stuff for the plus side */
      bfam_real_t ZsP = Zs_P[iG];
      bfam_real_t ZpP = Zp_P[iG];

      bfam_real_t TpP[] = {Tp1_P[iG], Tp2_P[iG], Tp3_P[iG]};
      bfam_real_t TnP = Tn_P[iG];

      bfam_real_t vpP[] = {vp1_P[iG],vp2_P[iG],vp3_P[iG]};
      bfam_real_t vnP = vn_P[iG];

      BEARD_STATE(
          &TnS_g[pnt],&TpS_g[3*pnt],&vnS_g[pnt],&vpS_g[3*pnt],
          TnM, TnP, TpM, TpP, vnM, vnP, vpM, vpP, ZpM, ZpP, ZsM, ZsP);

      /*
      BFAM_INFO("TP: %f %f %f %f",TnP,TpP[0],TpP[1],TpP[2]);
      BFAM_INFO("TM: %f %f %f %f",TnM,TpM[0],TpM[1],TpM[2]);
      BFAM_INFO("vP: %f %f %f %f",vnP,vpP[0],vpP[1],vpP[2]);
      BFAM_INFO("vM: %f %f %f %f",vnM,vpM[0],vpM[1],vpM[2]);
      */

      /* substract off the grid values */
      TpS_g[3*pnt+0] -= TpM[0];
      TpS_g[3*pnt+1] -= TpM[1];
      TpS_g[3*pnt+2] -= TpM[2];
      TnS_g[pnt]     -= TnM;
    }

    bfam_real_t *restrict TpS_M;
    bfam_real_t *restrict TnS_M;
    bfam_real_t *restrict vpS_M;
    bfam_real_t *restrict vnS_M;

    /* these will be used to the store the projected values if we need them */
    BFAM_ALIGN(32) bfam_real_t TpS_m_STORAGE[3*Nfp];
    BFAM_ALIGN(32) bfam_real_t TnS_m_STORAGE[  Nfp];
    BFAM_ALIGN(32) bfam_real_t vpS_m_STORAGE[3*Nfp];
    BFAM_ALIGN(32) bfam_real_t vnS_m_STORAGE[  Nfp];

    /* check to tee if projection */
    /* locked */
    {
      /* set to the correct Mass times projection */
      TpS_M = TpS_m_STORAGE;
      TnS_M = TnS_m_STORAGE;
      vpS_M = vpS_m_STORAGE;
      vnS_M = vnS_m_STORAGE;

      if(MASSPROJECTION && (!glue_p->same_order || glue_p->EToHm[le] || glue_p->EToHp[le]))
      {
        BFAM_ASSERT(glue_m->massprojection);
#if   DIM == 2
        const bfam_real_t *MP1 = glue_m->massprojection[glue_p->EToHm[le]];
#elif DIM == 3
        const int I1 = (glue_p->EToHm[le] == 0) ? 0 : (glue_p->EToHm[le]-1)/2+1;
        const int I2 = (glue_p->EToHm[le] == 0) ? 0 : (glue_p->EToHm[le]-1)%2+1;
        const bfam_real_t *MP1 = glue_m->massprojection[I1];
        const bfam_real_t *MP2 = glue_m->massprojection[I2];
#else
#error "Bad Dimension"
#endif
        beard_massproject_flux(TnS_M, TpS_M, vnS_M, vpS_M, N, Nrp_g, TnS_g, TpS_g,
                               vnS_g, vpS_g,
                               MP1,
#if DIM==3
                               MP2,
#endif
                               wi);
      }
      else if(PROJECTION  && (!glue_p->same_order || glue_p->EToHm[le] || glue_p->EToHp[le]))
      {
        BFAM_ASSERT(glue_m->projection);
#if   DIM == 2
        const bfam_real_t *P1 = glue_m->projection[glue_p->EToHm[le]];
#elif DIM == 3
        const int I1 = (glue_p->EToHm[le] == 0) ? 0 : (glue_p->EToHm[le]-1)/2+1;
        const int I2 = (glue_p->EToHm[le] == 0) ? 0 : (glue_p->EToHm[le]-1)%2+1;
        const bfam_real_t *P1 = glue_m->projection[I1];
        const bfam_real_t *P2 = glue_m->projection[I2];
#else
#error "Bad Dimension"
#endif
        beard_project_flux(TnS_M, TpS_M, vnS_M, vpS_M, N, Nrp_g, TnS_g, TpS_g,
                           vnS_g, vpS_g,
                           P1
#if DIM==3
                          ,P2
#endif
                           );
      }
      else
      {
        TpS_M = TpS_g;
        TnS_M = TnS_g;
        vpS_M = vpS_g;
        vnS_M = vnS_g;
      }
    }

    for(bfam_locidx_t pnt = 0; pnt < Nfp; pnt++)
    {
      bfam_locidx_t f = pnt + Nfp*(face + Nfaces*e);
      bfam_locidx_t iM = sub_m->vmapM[f];
      const bfam_real_t nM[] = {n1[f],n2[f],BEARD_D3_AP(0,+n3[f])};

      beard_dgx_add_flux(1,
          TnS_M[pnt],&TpS_M[3*pnt],vnS_M[pnt],&vpS_M[3*pnt],iM,
          dv1,dv2,dv3, dS11,dS22,dS33,dS12,dS13,dS23,
          lam[iM],mu[iM],rhoi[iM],nM,sJ[f],JI[iM],wi[0]);
    }
  }
}

void beard_dgx_inter_rhs_slip_weakening_interface(
    int inN, bfam_subdomain_dgx_t *sub_g, const char *rate_prefix,
    const char *minus_rate_prefix, const char *field_prefix,
    const bfam_long_real_t t)
{
  GENERIC_INIT(inN,beard_dgx_inter_rhs_slip_weakening_interface);

  bfam_subdomain_dgx_t* sub_m =
    (bfam_subdomain_dgx_t*) sub_g->base.glue_m->sub_m;

  /* get the fields we will need */
  bfam_subdomain_dgx_glue_data_t* glue_m =
    (bfam_subdomain_dgx_glue_data_t*) sub_g->base.glue_m;
  bfam_subdomain_dgx_glue_data_t* glue_p =
    (bfam_subdomain_dgx_glue_data_t*) sub_g->base.glue_p;
  BFAM_ASSERT(glue_m != NULL);
  BFAM_ASSERT(glue_p != NULL);
  bfam_dictionary_t *fields_m    = &glue_m->base.fields;
  bfam_dictionary_t *fields_p    = &glue_p->base.fields;
  bfam_dictionary_t *fields_g    = &sub_g ->base.fields;
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

  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(vn_M ,field_prefix,"vn" ,fields_m);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(vp1_M,field_prefix,"vp1",fields_m);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(vp2_M,field_prefix,"vp2",fields_m);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(vp3_M,field_prefix,"vp3",fields_m);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Tn_M ,field_prefix,"Tn" ,fields_m);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Tp1_M,field_prefix,"Tp1",fields_m);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Tp2_M,field_prefix,"Tp2",fields_m);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Tp3_M,field_prefix,"Tp3",fields_m);

  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(vn_P ,field_prefix,"vn" ,fields_p);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(vp1_P,field_prefix,"vp1",fields_p);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(vp2_P,field_prefix,"vp2",fields_p);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(vp3_P,field_prefix,"vp3",fields_p);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Tn_P ,field_prefix,"Tn" ,fields_p);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Tp1_P,field_prefix,"Tp1",fields_p);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Tp2_P,field_prefix,"Tp2",fields_p);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Tp3_P,field_prefix,"Tp3",fields_p);

  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(dv1 ,minus_rate_prefix,"v1" ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(dv2 ,minus_rate_prefix,"v2" ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(dv3 ,minus_rate_prefix,"v3" ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(dS11,minus_rate_prefix,"S11",fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(dS22,minus_rate_prefix,"S22",fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(dS33,minus_rate_prefix,"S33",fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(dS12,minus_rate_prefix,"S12",fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(dS13,minus_rate_prefix,"S13",fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(dS23,minus_rate_prefix,"S23",fields);

  /* get the material properties and metric terms */
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(rhoi,"","rho_inv"  ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(lam ,"","lam"      ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(mu  ,"","mu"       ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Zs  ,"","Zs"       ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Zp  ,"","Zp"       ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(J   ,"","_grid_J"  ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(JI  ,"","_grid_JI" ,fields);

  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(n1,"","_grid_nx0",fields_face);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(n2,"","_grid_nx1",fields_face);
  BEARD_D3_OP(BFAM_LOAD_FIELD_RESTRICT_ALIGNED(n3,"","_grid_nx2",fields_face));
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(sJ  ,"","_grid_sJ",fields_face);

  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Zs_M  ,"","Zs"       ,fields_m);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Zp_M  ,"","Zp"       ,fields_m);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Zs_P  ,"","Zs"       ,fields_p);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Zp_P  ,"","Zp"       ,fields_p);

  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Tp1_0   ,"","Tp1_0"   ,fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Tp2_0   ,"","Tp2_0"   ,fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Tp3_0   ,"","Tp3_0"   ,fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Tn_0    ,"","Tn_0"    ,fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Tp1     ,"","Tp1"     ,fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Tp2     ,"","Tp2"     ,fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Tp3     ,"","Tp3"     ,fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Tn      ,"","Tn"      ,fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(V       ,"","V"       ,fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Vp1     ,"","Vp1"     ,fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Vp2     ,"","Vp2"     ,fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Vp3     ,"","Vp3"     ,fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Dc      ,"","Dc"      ,fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Dp      ,"","Dp"      ,fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(fs      ,"","fs"      ,fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(fc      ,"","fc"      ,fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(fd      ,"","fd"      ,fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(c0      ,"","c0"      ,fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(T       ,"","Tforce"  ,fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(T0      ,"","Tforce_0",fields_g);

  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(dDp ,rate_prefix,"Dp", fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(dDn ,rate_prefix,"Dn", fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(dDp1,rate_prefix,"Dp1",fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(dDp2,rate_prefix,"Dp2",fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(dDp3,rate_prefix,"Dp3",fields_g);

  bfam_real_t *wi  = sub_m->wi;
  BFAM_ASSUME_ALIGNED(wi ,32);

  BFAM_ASSERT(glue_p->EToEm);
  BFAM_ASSERT(glue_p->EToFm);
  BFAM_ASSERT(glue_p->EToHm);

  for(bfam_locidx_t le = 0; le < sub_g->K; le++)
  {
    bfam_locidx_t e = glue_p->EToEm[le];
    int8_t face = glue_p->EToFm[le];


    if(glue_p->EToHm[le] < 2)
      beard_dgx_remove_flux(N,face,e,sub_m->vmapM,
          n1,n2,
#if DIM==3
          n3,
#endif
          Zs,Zp,
          mu,rhoi,lam,sJ,JI,wi,
          v1,v2,v3,S11,S22,S33,S12,S13,S23,
          dv1,dv2,dv3,dS11,dS22,dS33,dS12,dS13,dS23);


#ifndef USE_GENERIC
#undef Np
#undef N
#endif
    bfam_locidx_t Np_g  = sub_g->Np;
    bfam_locidx_t Nrp_g = sub_g->N+1;
#ifndef USE_GENERIC
#define Np Np_BACK
#define N  NORDER
#endif

    bfam_real_t TpS_g[3*Np_g];
    bfam_real_t TnS_g[  Np_g];
    bfam_real_t vpS_g[3*Np_g];
    bfam_real_t vnS_g[  Np_g];

    for(bfam_locidx_t pnt = 0; pnt < Np_g; pnt++)
    {
      bfam_locidx_t iG = pnt+le*Np_g;

      /* Setup stuff for the minus side */
      bfam_real_t ZsM = Zs_M[iG];
      bfam_real_t ZpM = Zp_M[iG];

      bfam_real_t TpM[] = {Tp1_M[iG], Tp2_M[iG], Tp3_M[iG]};
      bfam_real_t TnM = Tn_M[iG];

      bfam_real_t vpM[] = {vp1_M[iG],vp2_M[iG],vp3_M[iG]};
      bfam_real_t vnM = vn_M[iG];

      /* now add the real flux */
      /* Setup stuff for the plus side */
      bfam_real_t ZsP = Zs_P[iG];
      bfam_real_t ZpP = Zp_P[iG];

      bfam_real_t TpP[] = {Tp1_P[iG], Tp2_P[iG], Tp3_P[iG]};
      bfam_real_t TnP = Tn_P[iG];

      bfam_real_t vpP[] = {vp1_P[iG],vp2_P[iG],vp3_P[iG]};
      bfam_real_t vnP = vn_P[iG];

      /* compute the flux assume a locked fault */
      beard_dgx_upwind_state_m(
          &TnS_g[pnt],&TpS_g[3*pnt],&vnS_g[pnt],&vpS_g[3*pnt],
          TnM, TnP, TpM, TpP, vnM, vnP, vpM, vpP, ZpM, ZpP, ZsM, ZsP);

      Tn[iG] = TnS_g[pnt]+Tn_0[iG];
      if(NO_OPENING && Tn[iG] > 0)
      {
        BFAM_LOAD_FIELD_RESTRICT_ALIGNED(x  ,"","_grid_x0" ,fields_g);
        BFAM_LOAD_FIELD_RESTRICT_ALIGNED(y  ,"","_grid_x1" ,fields_g);
#if DIM==2
        BFAM_ABORT("fault opening not implemented: point"
            " %"BFAM_REAL_FMTe
            " %"BFAM_REAL_FMTe,
            x[iG], y[iG]);
#elif DIM==3
        BFAM_LOAD_FIELD_RESTRICT_ALIGNED(z  ,"","_grid_x2" ,fields_g);
        BFAM_ABORT("fault opening not implemented: point"
            " %"BFAM_REAL_FMTe
            " %"BFAM_REAL_FMTe
            " %"BFAM_REAL_FMTe,
            x[iG], y[iG], z[iG]);
#else
#error "invalid DIM"
#endif
      }

      const bfam_real_t Slock2 =
        + (TpS_g[3*pnt+0] + Tp1_0[iG])*(TpS_g[3*pnt+0] + Tp1_0[iG])
        + (TpS_g[3*pnt+1] + Tp2_0[iG])*(TpS_g[3*pnt+1] + Tp2_0[iG])
        + (TpS_g[3*pnt+2] + Tp3_0[iG])*(TpS_g[3*pnt+2] + Tp3_0[iG]);

      bfam_real_t f_val = BFAM_MIN(Dp[iG],Dc[iG])/Dc[iG];

      /* START SCEC
       * This forces the rupture for SCEC test problems :: eventually should
       * move to the input file
       */
      if(T0[iG] > 0)
        f_val = BFAM_MAX(f_val,
                         BFAM_MIN(BFAM_MAX((bfam_real_t)(t-T[iG])/T0[iG],0),1));
      /* END SCEC */

      fc[iG] = fs[iG]-(fs[iG]-fd[iG])*f_val;
      const bfam_real_t Sfric = c0[iG]-Tn[iG]*fc[iG];

      if(Sfric*Sfric < Slock2)
      {
        bfam_real_t VpS[3];
        const bfam_real_t Tp0[] =
            {Tp1_0[iG], Tp2_0[iG], Tp3_0[iG]};

        beard_dgx_upwind_state_friction_m(&TpS_g[3*pnt], &vpS_g[3*pnt], VpS,
            Sfric, TpM, TpP, Tp0, vpM, vpP, ZsM, ZsP);

        Vp1[iG] = VpS[0];
        Vp2[iG] = VpS[1];
        Vp3[iG] = VpS[2];
        V[iG]   = BFAM_REAL_SQRT(VpS[0]*VpS[0] + VpS[1]*VpS[1] + VpS[2]*VpS[2]);
      }
      else
      {
        Vp1[iG] = 0;
        Vp2[iG] = 0;
        Vp3[iG] = 0;
        V[iG]   = 0;
      }

      dDn[iG]  += 0;
      dDp[iG]  += V[iG];
      dDp1[iG] += Vp1[iG];
      dDp2[iG] += Vp2[iG];
      dDp3[iG] += Vp3[iG];
      Tp1[iG]   = TpS_g[3*pnt+0] + Tp1_0[iG];
      Tp2[iG]   = TpS_g[3*pnt+1] + Tp2_0[iG];
      Tp3[iG]   = TpS_g[3*pnt+2] + Tp3_0[iG];

      /* substract off the grid values */
      TpS_g[3*pnt+0] -= TpM[0];
      TpS_g[3*pnt+1] -= TpM[1];
      TpS_g[3*pnt+2] -= TpM[2];
      TnS_g[pnt]     -= TnM;
    }

    bfam_real_t *restrict TpS_M;
    bfam_real_t *restrict TnS_M;
    bfam_real_t *restrict vpS_M;
    bfam_real_t *restrict vnS_M;

    /* these will be used to the store the projected values if we need them */
    BFAM_ALIGN(32) bfam_real_t TpS_m_STORAGE[3*Nfp];
    BFAM_ALIGN(32) bfam_real_t TnS_m_STORAGE[  Nfp];
    BFAM_ALIGN(32) bfam_real_t vpS_m_STORAGE[3*Nfp];
    BFAM_ALIGN(32) bfam_real_t vnS_m_STORAGE[  Nfp];

    /* check to tee if projection */
    /* locked */
    {
      /* set to the correct Mass times projection */
      TpS_M = TpS_m_STORAGE;
      TnS_M = TnS_m_STORAGE;
      vpS_M = vpS_m_STORAGE;
      vnS_M = vnS_m_STORAGE;

      if(MASSPROJECTION && (!glue_p->same_order || glue_p->EToHm[le] || glue_p->EToHp[le]))
      {
        BFAM_ASSERT(glue_m->massprojection);
#if   DIM == 2
        const bfam_real_t *MP1 = glue_m->massprojection[glue_p->EToHm[le]];
#elif DIM == 3
        const int I1 = (glue_p->EToHm[le] == 0) ? 0 : (glue_p->EToHm[le]-1)/2+1;
        const int I2 = (glue_p->EToHm[le] == 0) ? 0 : (glue_p->EToHm[le]-1)%2+1;
        const bfam_real_t *MP1 = glue_m->massprojection[I1];
        const bfam_real_t *MP2 = glue_m->massprojection[I2];
#else
#error "Bad Dimension"
#endif
        beard_massproject_flux(TnS_M, TpS_M, vnS_M, vpS_M, N, Nrp_g, TnS_g, TpS_g,
                               vnS_g, vpS_g,
                               MP1,
#if DIM==3
                               MP2,
#endif
                               wi);
      }
      else if(PROJECTION  && (!glue_p->same_order || glue_p->EToHm[le] || glue_p->EToHp[le]))
      {
        BFAM_ASSERT(glue_m->projection);
#if   DIM == 2
        const bfam_real_t *P1 = glue_m->projection[glue_p->EToHm[le]];
#elif DIM == 3
        const int I1 = (glue_p->EToHm[le] == 0) ? 0 : (glue_p->EToHm[le]-1)/2+1;
        const int I2 = (glue_p->EToHm[le] == 0) ? 0 : (glue_p->EToHm[le]-1)%2+1;
        const bfam_real_t *P1 = glue_m->projection[I1];
        const bfam_real_t *P2 = glue_m->projection[I2];
#else
#error "Bad Dimension"
#endif
        beard_project_flux(TnS_M, TpS_M, vnS_M, vpS_M, N, Nrp_g, TnS_g, TpS_g,
                           vnS_g, vpS_g,
                           P1
#if DIM==3
                          ,P2
#endif
                           );
      }
      else
      {
        TpS_M = TpS_g;
        TnS_M = TnS_g;
        vpS_M = vpS_g;
        vnS_M = vnS_g;
      }
    }

    for(bfam_locidx_t pnt = 0; pnt < Nfp; pnt++)
    {
      bfam_locidx_t f = pnt + Nfp*(face + Nfaces*e);
      bfam_locidx_t iM = sub_m->vmapM[f];
      const bfam_real_t nM[] = {n1[f],n2[f],BEARD_D3_AP(0,+n3[f])};

      beard_dgx_add_flux(1,
          TnS_M[pnt],&TpS_M[3*pnt],vnS_M[pnt],&vpS_M[3*pnt],iM,
          dv1,dv2,dv3, dS11,dS22,dS33,dS12,dS13,dS23,
          lam[iM],mu[iM],rhoi[iM],nM,sJ[f],JI[iM],wi[0]);
    }
  }
}

void beard_dgx_inter_rhs_ageing_law_interface(
    int inN, bfam_subdomain_dgx_t *sub_g, const char *rate_prefix,
    const char *minus_rate_prefix, const char *field_prefix,
    const bfam_long_real_t t)
{
  GENERIC_INIT(inN,beard_dgx_inter_rhs_ageing_law_interface);

  bfam_subdomain_dgx_t* sub_m =
    (bfam_subdomain_dgx_t*) sub_g->base.glue_m->sub_m;

  /* get the fields we will need */
  bfam_subdomain_dgx_glue_data_t* glue_m =
    (bfam_subdomain_dgx_glue_data_t*) sub_g->base.glue_m;
  bfam_subdomain_dgx_glue_data_t* glue_p =
    (bfam_subdomain_dgx_glue_data_t*) sub_g->base.glue_p;
  BFAM_ASSERT(glue_m != NULL);
  BFAM_ASSERT(glue_p != NULL);
  bfam_dictionary_t *fields_m    = &glue_m->base.fields;
  bfam_dictionary_t *fields_p    = &glue_p->base.fields;
  bfam_dictionary_t *fields_g    = &sub_g ->base.fields;
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

  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(vn_M ,field_prefix,"vn" ,fields_m);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(vp1_M,field_prefix,"vp1",fields_m);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(vp2_M,field_prefix,"vp2",fields_m);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(vp3_M,field_prefix,"vp3",fields_m);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Tn_M ,field_prefix,"Tn" ,fields_m);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Tp1_M,field_prefix,"Tp1",fields_m);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Tp2_M,field_prefix,"Tp2",fields_m);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Tp3_M,field_prefix,"Tp3",fields_m);

  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(vn_P ,field_prefix,"vn" ,fields_p);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(vp1_P,field_prefix,"vp1",fields_p);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(vp2_P,field_prefix,"vp2",fields_p);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(vp3_P,field_prefix,"vp3",fields_p);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Tn_P ,field_prefix,"Tn" ,fields_p);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Tp1_P,field_prefix,"Tp1",fields_p);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Tp2_P,field_prefix,"Tp2",fields_p);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Tp3_P,field_prefix,"Tp3",fields_p);

  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(dv1 ,minus_rate_prefix,"v1" ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(dv2 ,minus_rate_prefix,"v2" ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(dv3 ,minus_rate_prefix,"v3" ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(dS11,minus_rate_prefix,"S11",fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(dS22,minus_rate_prefix,"S22",fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(dS33,minus_rate_prefix,"S33",fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(dS12,minus_rate_prefix,"S12",fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(dS13,minus_rate_prefix,"S13",fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(dS23,minus_rate_prefix,"S23",fields);

  /* get the material properties and metric terms */
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(rhoi,"","rho_inv"  ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(lam ,"","lam"      ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(mu  ,"","mu"       ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Zs  ,"","Zs"       ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Zp  ,"","Zp"       ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(J   ,"","_grid_J"  ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(JI  ,"","_grid_JI" ,fields);

  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(n1,"","_grid_nx0",fields_face);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(n2,"","_grid_nx1",fields_face);
  BEARD_D3_OP(BFAM_LOAD_FIELD_RESTRICT_ALIGNED(n3,"","_grid_nx2",fields_face));
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(sJ  ,"","_grid_sJ",fields_face);

  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Zs_M  ,"","Zs"       ,fields_m);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Zp_M  ,"","Zp"       ,fields_m);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Zs_P  ,"","Zs"       ,fields_p);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Zp_P  ,"","Zp"       ,fields_p);

  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Tp1_0 ,"","Tp1_0" ,fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Tp2_0 ,"","Tp2_0" ,fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Tp3_0 ,"","Tp3_0" ,fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Tn_0  ,"","Tn_0"  ,fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Tp1   ,"","Tp1"   ,fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Tp2   ,"","Tp2"   ,fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Tp3   ,"","Tp3"   ,fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Tn    ,"","Tn"    ,fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(V     ,"","V"     ,fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Vp1   ,"","Vp1"   ,fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Vp2   ,"","Vp2"   ,fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Vp3   ,"","Vp3"   ,fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(a     ,"","a"     ,fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(psi   ,"","psi"   ,fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(V0    ,"","V0"    ,fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(f0    ,"","f0"    ,fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(L     ,"","L"     ,fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(b     ,"","b"     ,fields_g);

  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(dpsi,rate_prefix,"psi",fields_g);

  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(dDp ,rate_prefix,"Dp", fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(dDn ,rate_prefix,"Dn", fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(dDp1,rate_prefix,"Dp1",fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(dDp2,rate_prefix,"Dp2",fields_g);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(dDp3,rate_prefix,"Dp3",fields_g);

  bfam_real_t *wi  = sub_m->wi;
  BFAM_ASSUME_ALIGNED(wi ,32);

  BFAM_ASSERT(glue_p->EToEm);
  BFAM_ASSERT(glue_p->EToFm);
  BFAM_ASSERT(glue_p->EToHm);

  for(bfam_locidx_t le = 0; le < sub_g->K; le++)
  {
    bfam_locidx_t e = glue_p->EToEm[le];
    int8_t face = glue_p->EToFm[le];


    if(glue_p->EToHm[le] < 2)
      beard_dgx_remove_flux(N,face,e,sub_m->vmapM,
          n1,n2,
#if DIM==3
          n3,
#endif
          Zs,Zp,
          mu,rhoi,lam,sJ,JI,wi,
          v1,v2,v3,S11,S22,S33,S12,S13,S23,
          dv1,dv2,dv3,dS11,dS22,dS33,dS12,dS13,dS23);


#ifndef USE_GENERIC
#undef Np
#undef N
#endif
    bfam_locidx_t Np_g  = sub_g->Np;
    bfam_locidx_t Nrp_g = sub_g->N+1;
#ifndef USE_GENERIC
#define Np Np_BACK
#define N  NORDER
#endif

    bfam_real_t TpS_g[3*Np_g];
    bfam_real_t TnS_g[  Np_g];
    bfam_real_t vpS_g[3*Np_g];
    bfam_real_t vnS_g[  Np_g];

    for(bfam_locidx_t pnt = 0; pnt < Np_g; pnt++)
    {
      bfam_locidx_t iG = pnt+le*Np_g;

      /* Setup stuff for the minus side */
      bfam_real_t ZsM = Zs_M[iG];
      bfam_real_t ZpM = Zp_M[iG];

      bfam_real_t TpM[] = {Tp1_M[iG], Tp2_M[iG], Tp3_M[iG]};
      bfam_real_t TnM = Tn_M[iG];

      bfam_real_t vpM[] = {vp1_M[iG],vp2_M[iG],vp3_M[iG]};
      bfam_real_t vnM = vn_M[iG];

      /* now add the real flux */
      /* Setup stuff for the plus side */
      bfam_real_t ZsP = Zs_P[iG];
      bfam_real_t ZpP = Zp_P[iG];

      bfam_real_t TpP[] = {Tp1_P[iG], Tp2_P[iG], Tp3_P[iG]};
      bfam_real_t TnP = Tn_P[iG];

      bfam_real_t vpP[] = {vp1_P[iG],vp2_P[iG],vp3_P[iG]};
      bfam_real_t vnP = vn_P[iG];

      /* compute the flux assume a locked fault */
      beard_dgx_upwind_state_m(
          &TnS_g[pnt],&TpS_g[3*pnt],&vnS_g[pnt],&vpS_g[3*pnt],
          TnM, TnP, TpM, TpP, vnM, vnP, vpM, vpP, ZpM, ZpP, ZsM, ZsP);

      Tn[iG] = TnS_g[pnt]+Tn_0[iG];
      if(NO_OPENING && Tn[iG] > 0)
      {
        BFAM_LOAD_FIELD_RESTRICT_ALIGNED(x  ,"","_grid_x0" ,fields_g);
        BFAM_LOAD_FIELD_RESTRICT_ALIGNED(y  ,"","_grid_x1" ,fields_g);
#if DIM==2
        BFAM_ABORT("fault opening not implemented: point"
            " %"BFAM_REAL_FMTe
            " %"BFAM_REAL_FMTe,
            x[iG], y[iG]);
#elif DIM==3
        BFAM_LOAD_FIELD_RESTRICT_ALIGNED(z  ,"","_grid_x2" ,fields_g);
        BFAM_ABORT("fault opening not implemented: point"
            " %"BFAM_REAL_FMTe
            " %"BFAM_REAL_FMTe
            " %"BFAM_REAL_FMTe,
            x[iG], y[iG], z[iG]);
#else
#error "invalid DIM"
#endif
      }

      /*
       * Call bracketed Newton solver to solve:
       *    Slock - eta*V = N*f(V,psi)
       * where eta = 2/(1/ZsP+1/ZsM)
       */
      const bfam_real_t Slock2 =
        + (TpS_g[3*pnt+0] + Tp1_0[iG])*(TpS_g[3*pnt+0] + Tp1_0[iG])
        + (TpS_g[3*pnt+1] + Tp2_0[iG])*(TpS_g[3*pnt+1] + Tp2_0[iG])
        + (TpS_g[3*pnt+2] + Tp3_0[iG])*(TpS_g[3*pnt+2] + Tp3_0[iG]);

      //JK fc[iG] = fs[iG]-(fs[iG]-fd[iG])*BFAM_MIN(Dp[iG],Dc[iG])/Dc[iG];
      //JK const bfam_real_t Sfric = c0[iG]-Tn[iG]*fc[iG];
      const bfam_real_t Sfric = 0;

      bfam_real_t VpS[3];
      const bfam_real_t Tp0[] =
      {Tp1_0[iG], Tp2_0[iG], Tp3_0[iG]};

      beard_dgx_upwind_state_rate_and_state_friction_m(
          &TpS_g[3*pnt], &vpS_g[3*pnt], VpS,
          Tn[pnt], a[iG],
          V0[iG] , psi[iG],
          TpM, TpP, Tp0,
          vpM, vpP,
          ZsM, ZsP);

      Vp1[iG] = VpS[0];
      Vp2[iG] = VpS[1];
      Vp3[iG] = VpS[2];
      V[iG]   = BFAM_REAL_SQRT(VpS[0]*VpS[0] + VpS[1]*VpS[1] + VpS[2]*VpS[2]);

      dDn[iG]  += 0;
      dDp[iG]  += V[iG];
      dDp1[iG] += Vp1[iG];
      dDp2[iG] += Vp2[iG];
      dDp3[iG] += Vp3[iG];
      Tp1[iG]   = TpS_g[3*pnt+0] + Tp1_0[iG];
      Tp2[iG]   = TpS_g[3*pnt+1] + Tp2_0[iG];
      Tp3[iG]   = TpS_g[3*pnt+2] + Tp3_0[iG];

      const bfam_real_t theta =
        (L[iG]/V0[iG])*BFAM_REAL_EXP((psi[iG]-f0[iG])/b[iG]);
      dpsi[iG] += b[iG]*(1.0/theta-V[iG]/L[iG]);
      /* dpsi[iG] += (b[iG]/theta)*(1-V[iG]*theta/L[iG]); */

      /* substract off the grid values */
      TpS_g[3*pnt+0] -= TpM[0];
      TpS_g[3*pnt+1] -= TpM[1];
      TpS_g[3*pnt+2] -= TpM[2];
      TnS_g[pnt]     -= TnM;
    }

    bfam_real_t *restrict TpS_M;
    bfam_real_t *restrict TnS_M;
    bfam_real_t *restrict vpS_M;
    bfam_real_t *restrict vnS_M;

    /* these will be used to the store the projected values if we need them */
    BFAM_ALIGN(32) bfam_real_t TpS_m_STORAGE[3*Nfp];
    BFAM_ALIGN(32) bfam_real_t TnS_m_STORAGE[  Nfp];
    BFAM_ALIGN(32) bfam_real_t vpS_m_STORAGE[3*Nfp];
    BFAM_ALIGN(32) bfam_real_t vnS_m_STORAGE[  Nfp];

    /* check to tee if projection */
    /* locked */
    {
      /* set to the correct Mass times projection */
      TpS_M = TpS_m_STORAGE;
      TnS_M = TnS_m_STORAGE;
      vpS_M = vpS_m_STORAGE;
      vnS_M = vnS_m_STORAGE;

      if(MASSPROJECTION && (!glue_p->same_order || glue_p->EToHm[le] ||
                             glue_p->EToHp[le]))
      {
        BFAM_ASSERT(glue_m->massprojection);
#if   DIM == 2
        const bfam_real_t *MP1 = glue_m->massprojection[glue_p->EToHm[le]];
#elif DIM == 3
        const int I1 = (glue_p->EToHm[le] == 0) ? 0 : (glue_p->EToHm[le]-1)/2+1;
        const int I2 = (glue_p->EToHm[le] == 0) ? 0 : (glue_p->EToHm[le]-1)%2+1;
        const bfam_real_t *MP1 = glue_m->massprojection[I1];
        const bfam_real_t *MP2 = glue_m->massprojection[I2];
#else
#error "Bad Dimension"
#endif
        beard_massproject_flux(TnS_M, TpS_M, vnS_M, vpS_M, N, Nrp_g,
                               TnS_g, TpS_g,
                               vnS_g, vpS_g,
                               MP1,
#if DIM==3
                               MP2,
#endif
                               wi);
      }
      else if(PROJECTION  && (!glue_p->same_order || glue_p->EToHm[le] ||
                               glue_p->EToHp[le]))
      {
        BFAM_ASSERT(glue_m->projection);
#if   DIM == 2
        const bfam_real_t *P1 = glue_m->projection[glue_p->EToHm[le]];
#elif DIM == 3
        const int I1 = (glue_p->EToHm[le] == 0) ? 0 : (glue_p->EToHm[le]-1)/2+1;
        const int I2 = (glue_p->EToHm[le] == 0) ? 0 : (glue_p->EToHm[le]-1)%2+1;
        const bfam_real_t *P1 = glue_m->projection[I1];
        const bfam_real_t *P2 = glue_m->projection[I2];
#else
#error "Bad Dimension"
#endif
        beard_project_flux(TnS_M, TpS_M, vnS_M, vpS_M, N, Nrp_g, TnS_g, TpS_g,
                           vnS_g, vpS_g,
                           P1
#if DIM==3
                          ,P2
#endif
                           );
      }
      else
      {
        TpS_M = TpS_g;
        TnS_M = TnS_g;
        vpS_M = vpS_g;
        vnS_M = vnS_g;
      }
    }

    for(bfam_locidx_t pnt = 0; pnt < Nfp; pnt++)
    {
      bfam_locidx_t f = pnt + Nfp*(face + Nfaces*e);
      bfam_locidx_t iM = sub_m->vmapM[f];
      const bfam_real_t nM[] = {n1[f],n2[f],BEARD_D3_AP(0,+n3[f])};

      beard_dgx_add_flux(1,
          TnS_M[pnt],&TpS_M[3*pnt],vnS_M[pnt],&vpS_M[3*pnt],iM,
          dv1,dv2,dv3, dS11,dS22,dS33,dS12,dS13,dS23,
          lam[iM],mu[iM],rhoi[iM],nM,sJ[f],JI[iM],wi[0]);
    }
  }
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
          bfam_real_t mean_stress = (S11[n]+S22[n]+S33[n])/3;
          bfam_real_t s11 = S11[n] - mean_stress;
          bfam_real_t s22 = S22[n] - mean_stress;
          bfam_real_t s33 = S33[n] - mean_stress;
          bfam_real_t s12 = S12[n];
          bfam_real_t s13 = S13[n];
          bfam_real_t s23 = S23[n];

          /* bulk modulus */
          bfam_real_t K = lam[n] + 2*mu[n]/3;

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
