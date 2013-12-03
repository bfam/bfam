#ifdef BEARD_DGX_DIMENSION

#if    BEARD_DGX_DIMENSION==2
#include "beard_dgx_rhs_2.h"

#define DIM 2

#ifndef NORDER
#define NORDER
#define USE_GENERIC
#define GENERIC_INIT(INN,NFM)\
  int N   = INN; \
  int Np  = (INN+1)*(INN+1); \
  int Nfp = (INN+1); \
  BFAM_WARNING(\
      "using generic NFM for 2d of order %d with Np = %d and Nfp = %d",\
      N, Np, Nfp);
#else
#define GENERIC_INIT(INN,NFM) BFAM_NOOP()
#define N    NORDER
#define Np  (NORDER+1)*(NORDER+1)
#define Nfp (NORDER+1)
#endif

#elif  BEARD_DGX_DIMENSION==3
#include "beard_dgx_rhs_3.h"

#define DIM 3

#ifndef NORDER
#define NORDER
#define USE_GENERIC
#define GENERIC_INIT(INN,NFM)\
  int N   = INN; \
  int Np  = (INN+1)*(INN+1)*(INN+1); \
  int Nfp = (INN+1)*(INN+1); \
  BFAM_WARNING(\
      "using generic NFM for 3d of order %d with Np = %d and Nfp = %d",\
      N, Np, Nfp);
#else
#define GENERIC_INIT(INN,NFM) BFAM_NOOP()
#define N    NORDER
#define Np  (NORDER+1)*(NORDER+1)*(NORDER+1)
#define Nfp (NORDER+1)*(NORDER+1)
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



void beard_dgx_intra_rhs_elastic(
    int inN, bfam_subdomain_dgx_t *sub, const char *rate_prefix,
    const char *field_prefix, const bfam_long_real_t t)
{
  GENERIC_INIT(inN,beard_dgx_intra_rhs_elastic);

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

}

#endif
