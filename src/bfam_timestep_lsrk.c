#include <bfam_timestep_lsrk.h>

typedef struct bfam_ts_lsrk_sd
{
  bfam_subdomain_t *subdomain; /**< subdomain this guy updates */
  void *fields; /**< pointer to the fields I update */
  void *rates;  /**< pointer to the rates I use */

  /**< scale the rates */
  void (*scale_rates) (bfam_subdomain_t *thisSubdomain, void *rates,
      const bfam_real_t a);

  /**< Do the intra work for this subdomain: dq += RHS(q,t) */
  void (*intra_rhs) (bfam_subdomain_t *thisSubdomain, void *rates,
      const void *fields, const bfam_real_t t);

  /**< Do the inter work for this subdomain: dq += RHS(q,t) */
  void (*inter_rhs) (bfam_subdomain_t *thisSubdomain, void *rates,
      const void *fields, const bfam_real_t t);

  /**< add rates to fields: q += b*dq */
  void (*add_rates) (bfam_subdomain_t *thisSubdomain, void *fields,
      const void *rates, const bfam_real_t a);

} bfam_ts_lsrk_sd_t;

bfam_ts_lsrk_t*
bfam_ts_lsrk_new(bfam_domain_t* dom, bfam_communicator_t *comm,
    bfam_ts_lsrk_method_t method)
{
  bfam_ts_lsrk_t* newTS = bfam_malloc(sizeof(bfam_ts_lsrk_t));
  bfam_ts_lsrk_init(newTS, dom, comm, method);
  return newTS;
}


void
bfam_ts_lsrk_init(bfam_ts_lsrk_t* ts, bfam_domain_t* dom,
    bfam_communicator_t *comm, bfam_ts_lsrk_method_t method)
{
  bfam_ts_init(&ts->base, dom);
  bfam_dictionary_init(&ts->elems);
  ts->t  = 0.0;
  ts->comm = comm;
  switch(method)
  {
    case BFAM_TS_LSRK_KC54:
      ts->nStages = 5;
      ts->A = bfam_malloc_aligned(ts->nStages*sizeof(bfam_long_real_t));
      ts->B = bfam_malloc_aligned(ts->nStages*sizeof(bfam_long_real_t));
      ts->C = bfam_malloc_aligned((ts->nStages+1)*sizeof(bfam_long_real_t));

      ts->A[0] = 0;
      ts->A[1] = -567301805773.0/1357537059087.0;
      ts->A[2] = -2404267990393.0/2016746695238.0;
      ts->A[3] = -3550918686646.0/2091501179385.0;
      ts->A[4] = -1275806237668.0/842570457699.0;

      ts->B[0] = 1432997174477.0/9575080441755.0;
      ts->B[1] = 5161836677717.0/13612068292357.0;
      ts->B[2] = 1720146321549.0/2090206949498.0;
      ts->B[3] = 3134564353537.0/4481467310338.0;
      ts->B[4] = 2277821191437.0/14882151754819.0;

      ts->C[0] = 0.0;
      ts->C[1] = 1432997174477.0/9575080441755.0;
      ts->C[2] = 2526269341429.0/6820363962896.0;
      ts->C[3] = 2006345519317.0/3224310063776.0;
      ts->C[4] = 2802321613138.0/2924317926251.0;
      ts->C[5] = 1.0;
      break;
    case BFAM_TS_LSRK_W33:
      ts->nStages = 3;
      ts->A = bfam_malloc_aligned(ts->nStages*sizeof(bfam_long_real_t));
      ts->B = bfam_malloc_aligned(ts->nStages*sizeof(bfam_long_real_t));
      ts->C = bfam_malloc_aligned((ts->nStages+1)*sizeof(bfam_long_real_t));

      ts->A[0] = 0.0;
      ts->A[1] = -5.0/9.0;
      ts->A[2] = -153.0/128.0;

      ts->B[0] = 1.0/3.0;
      ts->B[1] = 15.0/16.0;
      ts->B[2] = 8.0/15.0;

      ts->C[0] = 0.0;
      ts->C[1] = 1.0/3.0;
      ts->C[2] = 3.0/4.0;
      ts->C[3] = 1.0;
      break;
    case BFAM_TS_LSRK_HEUN:
      ts->nStages = 2;
      ts->A = bfam_malloc_aligned(ts->nStages*sizeof(bfam_long_real_t));
      ts->B = bfam_malloc_aligned(ts->nStages*sizeof(bfam_long_real_t));
      ts->C = bfam_malloc_aligned((ts->nStages+1)*sizeof(bfam_long_real_t));

      ts->A[0] = 0.0;
      ts->A[1] = -1.0;

      ts->B[0] = 1.0;
      ts->B[1] = 1.0/2.0;

      ts->C[0] = 0.0;
      ts->C[1] = 1.0;
      ts->C[2] = 1.0;
      break;
    case BFAM_TS_LSRK_FE:
      ts->nStages = 1;
      ts->A = bfam_malloc_aligned(sizeof(bfam_long_real_t));
      ts->B = bfam_malloc_aligned(sizeof(bfam_long_real_t));
      ts->C = bfam_malloc_aligned((ts->nStages+1)*sizeof(bfam_long_real_t));

      ts->A[0] = 0.0;

      ts->B[0] = 1.0;

      ts->C[0] = 0.0;
      ts->C[1] = 1.0;
      break;
  }
}

void
bfam_ts_lsrk_free(bfam_ts_lsrk_t* ts)
{
  bfam_dictionary_clear(&ts->elems);
  bfam_free_aligned(ts->A);
  bfam_free_aligned(ts->B);
  bfam_free_aligned(ts->C);
  ts->nStages = 0;
  ts->t  = NAN;
  bfam_ts_free(&ts->base);
}

void
bfam_ts_lsrk_set_time(bfam_ts_lsrk_t* ts,bfam_long_real_t time)
{
  ts->t = time;
}

bfam_long_real_t
bfam_ts_lsrk_get_time(bfam_ts_lsrk_t* ts)
{
  return ts->t;
}

void
bfam_ts_lsrk_add_subdomains(bfam_ts_lsrk_t *ts, bfam_domain_match_t match,
                            const char *tags[], const char *fields[],
  void (*scale_rates) (bfam_subdomain_t *thisSubdomain, void *rates,
      const bfam_real_t a),
  void (*intra_rhs) (bfam_subdomain_t *thisSubdomain, void *rates,
      const void *fields, const bfam_real_t t),
  void (*inter_rhs) (bfam_subdomain_t *thisSubdomain, void *rates,
      const void *fields, const bfam_real_t t),
  void (*add_rates) (bfam_subdomain_t *thisSubdomain, void *fields,
      const void *rates, const bfam_real_t a))
{
}
