#include "beard_dgx_rhs.h"

#ifdef BEARD_DGX_DIMENSION

#if    BEARD_DGX_DIMENSION==2

#define DIM 2

#ifndef NORDER
#define NORDER
#define USE_GENERIC
#define GENERIC_NORDER(INN,FNAME)\
  int N   = INN; \
  int Np  = (INN+1)*(INN+1); \
  int Nfp = (INN+1); \
  BFAM_WARNING("using generic FNAME for 2d of order %d",(int)INN);
#else
#define N    NORDER
#define Np  (NORDER+1)*(NORDER+1)
#define Nfp (NORDER+1)
#endif

#elif  BEARD_DGX_DIMENSION==3

#define DIM 3

#ifndef NORDER
#define NORDER
#define USE_GENERIC
#define GENERIC_NORDER(INN,FNAME)\
  int N   = INN; \
  int Np  = (INN+1)*(INN+1)*(INN+1); \
  int Nfp = (INN+1)*(INN+1); \
  BFAM_WARNING("using generic FNAME for 3d of order %d",(int)INN);
#else
#define N    NORDER
#define Np  (NORDER+1)*(NORDER+1)*(NORDER+1)
#define Nfp (NORDER+1)*(NORDER+1)
#endif

#else
#error "bad dimension"
#endif

#define BEARD_APPEND(a,b,c,d) a ## b ## c ##d
#define BEARD_APPEND_EXPAND(a,b,c,d) BEARD_APPEND(a,b,c,d)
#define beard_dgx_intra_rhs_elastic \
  BEARD_APPEND_EXPAND(beard_dgx_intra_rhs_elastic_,DIM,_,NORDER)

#define beard_dgx_scale_rates_elastic \
  BEARD_APPEND_EXPAND(beard_dgx_scale_rates_elastic_,DIM,_,NORDER)

#define beard_dgx_add_rates_elastic \
  BEARD_APPEND_EXPAND(beard_dgx_add_rates_elastic_,DIM,_,NORDER)

#define beard_dgx_inter_rhs_boundary \
  BEARD_APPEND_EXPAND(beard_dgx_inter_rhs_boundary_,DIM,_,NORDER)

#define beard_dgx_inter_rhs_interface \
  BEARD_APPEND_EXPAND(beard_dgx_inter_rhs_interface_,DIM,_,NORDER)

#define beard_dgx_energy \
  BEARD_APPEND_EXPAND(beard_dgx_energy_,DIM,_,NORDER)


void beard_dgx_intra_rhs_elastic(
    int inN, bfam_subdomain_dgx_quad_t *sub, const char *rate_prefix,
    const char *field_prefix, const bfam_long_real_t t)
{
#ifdef USE_GENERIC
  GENERIC_INIT(inN,beard_dgx_intra_rhs_elastic);
#endif
}

void beard_dgx_scale_rates_elastic(
    int inN, bfam_subdomain_dgx_quad_t *sub, const char *rate_prefix,
    const bfam_long_real_t a)
{
#ifdef USE_GENERIC
  GENERIC_INIT(inN,beard_dgx_scale_rates_elastic);
#endif
}

void beard_dgx_add_rates_elastic(
    int inN, bfam_subdomain_dgx_quad_t *sub, const char *field_prefix_lhs,
    const char *field_prefix_rhs, const char *rate_prefix,
    const bfam_long_real_t a)
{
#ifdef USE_GENERIC
  GENERIC_INIT(inN,beard_dgx_add_rates_elastic);
#endif
}

void beard_dgx_inter_rhs_boundary(
    int inN, bfam_subdomain_dgx_quad_glue_t *sub_g, const char *rate_prefix,
    const char *field_prefix, const bfam_long_real_t t, const bfam_real_t R)
{
#ifdef USE_GENERIC
  GENERIC_INIT(inN,beard_dgx_inter_rhs_boundary);
#endif
}

void beard_dgx_inter_rhs_interface(
    int inN, bfam_subdomain_dgx_quad_glue_t *sub_g, const char *rate_prefix,
    const char *field_prefix, const bfam_long_real_t t)
{
#ifdef USE_GENERIC
  GENERIC_INIT(inN,beard_dgx_inter_rhs_interface);
#endif
}

void beard_dgx_energy(
    int inN, bfam_real_t *energy_sq,
    bfam_subdomain_dgx_quad_t *sub, const char *field_prefix)
{
#ifdef USE_GENERIC
  GENERIC_INIT(inN,beard_dgx_energy);
#endif
}

#endif
