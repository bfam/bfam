#include <bfam_timestep_lsrk.h>
#include <bfam_log.h>

#define BFAM_LSKR_PREFIX ("_lsrk_rate_")

typedef struct bfam_ts_allprefix
{
  bfam_ts_lsrk_t *ts;
  bfam_long_real_t arg;
} bfam_ts_lsrk_allprefix_t;

static int
bfam_ts_lsrk_scale_rates(const char * key, void *val, void *arg)
{
  bfam_ts_lsrk_allprefix_t *data = (bfam_ts_lsrk_allprefix_t *) arg;
  bfam_subdomain_t* sub = (bfam_subdomain_t*) val;
  data->ts->scale_rates(sub, BFAM_LSKR_PREFIX, data->arg);
  return 1;
}

static int
bfam_ts_lsrk_intra_rhs(const char * key, void *val, void *arg)
{
  bfam_ts_lsrk_allprefix_t *data = (bfam_ts_lsrk_allprefix_t *) arg;
  bfam_subdomain_t* sub = (bfam_subdomain_t*) val;
  data->ts->intra_rhs(sub, BFAM_LSKR_PREFIX, "", data->arg);
  return 1;
}

static int
bfam_ts_lsrk_inter_rhs(const char * key, void *val, void *arg)
{
  bfam_ts_lsrk_allprefix_t *data = (bfam_ts_lsrk_allprefix_t *) arg;
  bfam_subdomain_t* sub = (bfam_subdomain_t*) val;
  data->ts->inter_rhs(sub, BFAM_LSKR_PREFIX, "", data->arg);
  return 1;
}

static int
bfam_ts_lsrk_add_rates(const char * key, void *val, void *arg)
{
  bfam_ts_lsrk_allprefix_t *data = (bfam_ts_lsrk_allprefix_t *) arg;
  bfam_subdomain_t* sub = (bfam_subdomain_t*) val;
  data->ts->add_rates(sub, "", "", BFAM_LSKR_PREFIX, data->arg);
  return 1;
}



void
bfam_ts_lsrk_step(bfam_ts_t *a_ts, bfam_long_real_t dt)
{
  bfam_ts_lsrk_t *ts = (bfam_ts_lsrk_t*) a_ts;
  bfam_ts_lsrk_allprefix_t data;
  data.ts = ts;

  for(int s = 0; s < ts->nStages;s++)
  {
    /* set the stage time */
    bfam_long_real_t t = ts->t + ts->C[s]*dt;

    /*
     * start the communication
     */
    bfam_communicator_start(ts->comm);

    /*
     * scale the rate
     */
    data.arg = ts->A[s];
    bfam_dictionary_allprefixed_ptr(&ts->elems,
        "",&bfam_ts_lsrk_scale_rates,&data);

    /*
     * do the intra work
     */
    data.arg = t;
    bfam_dictionary_allprefixed_ptr(&ts->elems,
        "",&bfam_ts_lsrk_intra_rhs,&data);

    /*
     * finish the communication
     */
    bfam_communicator_finish(ts->comm);

    /*
     * do the inter work
     */
    data.arg = t;
    bfam_dictionary_allprefixed_ptr(&ts->elems,
        "",&bfam_ts_lsrk_inter_rhs,&data);

    /*
     * add the rate
     */
    data.arg = ts->B[s]*dt;
    bfam_dictionary_allprefixed_ptr(&ts->elems,
        "",&bfam_ts_lsrk_add_rates,&data);
  }
  ts->t += ts->C[ts->nStages]*dt;
}

bfam_ts_lsrk_t*
bfam_ts_lsrk_new(bfam_domain_t* dom, bfam_ts_lsrk_method_t method,
    bfam_domain_match_t subdom_match, const char** subdom_tags,
    bfam_domain_match_t comm_match, const char** comm_tags,
    MPI_Comm mpicomm, int mpitag,
    void (*aux_rates) (bfam_subdomain_t *thisSubdomain, const char *prefix),
    void (*scale_rates) (bfam_subdomain_t *thisSubdomain,
      const char *rate_prefix, const bfam_long_real_t a),
    void (*intra_rhs) (bfam_subdomain_t *thisSubdomain,
      const char *rate_prefix, const char *field_prefix,
      const bfam_long_real_t t),
    void (*inter_rhs) (bfam_subdomain_t *thisSubdomain,
      const char *rate_prefix, const char *field_prefix,
      const bfam_long_real_t t),
    void (*add_rates) (bfam_subdomain_t *thisSubdomain,
      const char *field_prefix_lhs, const char *field_prefix_rhs,
      const char *rate_prefix, const bfam_long_real_t a))
{
  bfam_ts_lsrk_t* newTS = bfam_malloc(sizeof(bfam_ts_lsrk_t));
  bfam_ts_lsrk_init(newTS, dom, method, subdom_match, subdom_tags,
      comm_match, comm_tags, mpicomm, mpitag, aux_rates,
      scale_rates,intra_rhs,inter_rhs,add_rates);
  return newTS;
}


void
bfam_ts_lsrk_init(bfam_ts_lsrk_t* ts,
    bfam_domain_t* dom, bfam_ts_lsrk_method_t method,
    bfam_domain_match_t subdom_match, const char** subdom_tags,
    bfam_domain_match_t comm_match, const char** comm_tags,
    MPI_Comm mpicomm, int mpitag,
    void (*aux_rates) (bfam_subdomain_t *thisSubdomain, const char *prefix),
    void (*scale_rates) (bfam_subdomain_t *thisSubdomain,
      const char *rate_prefix, const bfam_long_real_t a),
    void (*intra_rhs) (bfam_subdomain_t *thisSubdomain,
      const char *rate_prefix, const char *field_prefix,
      const bfam_long_real_t t),
    void (*inter_rhs) (bfam_subdomain_t *thisSubdomain,
      const char *rate_prefix, const char *field_prefix,
      const bfam_long_real_t t),
    void (*add_rates) (bfam_subdomain_t *thisSubdomain,
      const char *field_prefix_lhs, const char *field_prefix_rhs,
      const char *rate_prefix, const bfam_long_real_t a))
{
  BFAM_LDEBUG("LSRK INIT");

  /*
   * set up some preliminaries
   */
  bfam_ts_init(&ts->base, dom);
  bfam_dictionary_init(&ts->elems);
  ts->t  = 0.0;
  ts->base.step = &bfam_ts_lsrk_step;

  /*
   * store the function calls
   */
  ts->scale_rates = scale_rates;
  ts->intra_rhs   = intra_rhs;
  ts->inter_rhs   = inter_rhs;
  ts->add_rates   = add_rates;

  /*
   * get the subdomains and create rates we will need
   */
   bfam_subdomain_t *subs[dom->numSubdomains+1];
   bfam_locidx_t numSubs = 0;
   bfam_domain_get_subdomains(dom,subdom_match,subdom_tags,
       dom->numSubdomains,subs,&numSubs);
   for(int s = 0; s < numSubs;s++)
   {
     int rval = bfam_dictionary_insert_ptr(&ts->elems,subs[s]->name,subs[s]);
     BFAM_ASSERT(rval != 1);

     aux_rates(subs[s],BFAM_LSKR_PREFIX);
   }

  /*
   * Set up the communicator we will use
   */
   ts->comm = bfam_communicator_new(dom,comm_match,comm_tags,mpicomm,mpitag);

  switch(method)
  {
    default:
      BFAM_WARNING("Invalid LSRK scheme, using KC54");
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
  BFAM_LDEBUG("LSRK FREE");
  bfam_communicator_free(ts->comm);
  bfam_free(ts->comm);
  ts->comm = NULL;
  bfam_dictionary_clear(&ts->elems);
  bfam_free_aligned(ts->A);
  ts->A = NULL;
  bfam_free_aligned(ts->B);
  ts->B = NULL;
  bfam_free_aligned(ts->C);
  ts->C = NULL;
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
