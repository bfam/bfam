#include <bfam_timestep_adams.h>
#include <bfam_log.h>

#define BFAM_ADAMS_PREFIX ("_adams_rate_")

typedef struct bfam_ts_allprefix
{
  bfam_ts_adams_t *ts;
  bfam_long_real_t dt;
} bfam_ts_adams_allprefix_t;

bfam_ts_adams_t*
bfam_ts_adams_new(bfam_domain_t* dom, bfam_ts_adams_method_t method,
    bfam_domain_match_t subdom_match, const char** subdom_tags,
    bfam_domain_match_t comm_match, const char** comm_tags,
    MPI_Comm mpicomm, int mpitag, void *comm_data,
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
      const char *rate_prefix, const bfam_long_real_t a),
    const int RK_init)
{
  bfam_ts_adams_t* newTS = bfam_malloc(sizeof(bfam_ts_adams_t));
  bfam_ts_adams_init(newTS, dom, method, subdom_match, subdom_tags,
      comm_match, comm_tags, mpicomm, mpitag, comm_data, aux_rates,
      scale_rates,intra_rhs,inter_rhs,add_rates,RK_init);
  return newTS;
}

static inline int
bfam_ts_adams_do_update(bfam_subdomain_t* sub, const bfam_long_real_t* A,
    const bfam_ts_adams_t* ts, const bfam_long_real_t dt, const int nStages)
{
  BFAM_LDEBUG("BFAM_TS_ADAMS_DO_UPDATE");

  /* Loop through the stages to scale rates and add in */
  /*
   * nStages is the computing number of stages whereas ts->nStages is the
   * storage number of stages
   */
  for(int k = 0; k < nStages;k++)
  {
    char prefix[BFAM_BUFSIZ];
    snprintf(prefix,BFAM_BUFSIZ,"%s%d_",BFAM_ADAMS_PREFIX,
        (ts->currentStage+ts->nStages-k)%ts->nStages);
    BFAM_LDEBUG("Adams step: stage %d of %d using prefix %s",k,nStages,prefix);
    ts->add_rates(sub, "", "", prefix, dt*A[k]);
  }
  return 1;
}

static int
bfam_ts_adams_update(const char * key, void *val, void *arg)
{
  bfam_ts_adams_allprefix_t *data = (bfam_ts_adams_allprefix_t *) arg;
  bfam_subdomain_t* sub = (bfam_subdomain_t*) val;
  bfam_ts_adams_t* ts = data->ts;

  switch(BFAM_MIN(ts->numSteps+1,ts->nStages))
  {
    case 1:
      {
        bfam_long_real_t A[1] = {BFAM_LONG_REAL(1.0)};
        bfam_ts_adams_do_update(sub, A, ts, data->dt, 1);
      }
      break;
    case 2:
      {
        bfam_long_real_t A[2] = {
          BFAM_LONG_REAL( 3.0) / BFAM_LONG_REAL( 2.0),
          BFAM_LONG_REAL(-1.0) / BFAM_LONG_REAL( 2.0),
        };
        bfam_ts_adams_do_update(sub, A, ts, data->dt, 2);
      }
      break;
    case 3:
      {
        bfam_long_real_t A[3] = {
          BFAM_LONG_REAL(23.0)/ BFAM_LONG_REAL(12.0),
          BFAM_LONG_REAL(-4.0)/ BFAM_LONG_REAL( 3.0),
          BFAM_LONG_REAL( 5.0)/ BFAM_LONG_REAL(12.0),
        };
        bfam_ts_adams_do_update(sub, A, ts, data->dt, 3);
      }
      break;
    case 4:
      {
        bfam_long_real_t A[4] = {
          BFAM_LONG_REAL( 55.0)/ BFAM_LONG_REAL( 24.0),
          BFAM_LONG_REAL(-59.0)/ BFAM_LONG_REAL( 24.0),
          BFAM_LONG_REAL( 37.0)/ BFAM_LONG_REAL( 24.0),
          BFAM_LONG_REAL(  3.0)/ BFAM_LONG_REAL(  8.0),
        };
        bfam_ts_adams_do_update(sub, A, ts, data->dt, 4);
      }
      break;
    default:
      BFAM_ABORT("Adams-Bashforth order %d not implemented",ts->nStages);
  }
  return 1;
}

static int
bfam_ts_adams_intra_rhs(const char * key, void *val, void *arg)
{
  bfam_ts_adams_allprefix_t *data = (bfam_ts_adams_allprefix_t *) arg;
  bfam_subdomain_t* sub = (bfam_subdomain_t*) val;
  char prefix[BFAM_BUFSIZ];
  snprintf(prefix,BFAM_BUFSIZ,"%s%d_",BFAM_ADAMS_PREFIX,
      data->ts->currentStage%data->ts->nStages);
  BFAM_LDEBUG("Adams intra: using prefix %s",prefix);
  data->ts->scale_rates(sub, prefix, 0);
  data->ts->intra_rhs(sub, prefix, "", data->ts->t);
  return 1;
}

static int
bfam_ts_adams_inter_rhs(const char * key, void *val, void *arg)
{
  bfam_ts_adams_allprefix_t *data = (bfam_ts_adams_allprefix_t *) arg;
  bfam_subdomain_t* sub = (bfam_subdomain_t*) val;
  char prefix[BFAM_BUFSIZ];
  snprintf(prefix,BFAM_BUFSIZ,"%s%d_",BFAM_ADAMS_PREFIX,
      data->ts->currentStage%data->ts->nStages);
  BFAM_LDEBUG("Adams inter: using prefix %s",prefix);
  data->ts->inter_rhs(sub, prefix, "", data->ts->t);
  return 1;
}


void
bfam_ts_adams_step(bfam_ts_t *a_ts, bfam_long_real_t dt)
{
  /* cast to the proper type */
  bfam_ts_adams_t *ts = (bfam_ts_adams_t*) a_ts;

  bfam_ts_adams_allprefix_t data;
  data.ts = ts;
  data.dt = dt;

  /*
   * start the communication
   */
  bfam_communicator_start(ts->comm);

  /*
   * do the intra work
   */
  bfam_dictionary_allprefixed_ptr(&ts->elems,
      "",&bfam_ts_adams_intra_rhs,&data);

  /*
   * finish the communication
   */
  bfam_communicator_finish(ts->comm);

  /*
   * do the inter work
   */
  bfam_dictionary_allprefixed_ptr(&ts->elems,
      "",&bfam_ts_adams_inter_rhs,&data);


  if(ts->lsrk)
  {
    /*
     * If we are using RK, we want to tell the RK scheme to use the next rate
     * for storage (not the current rate which is valid
     */
    ts->lsrk->t = ts->t;
    char rate_prefix[BFAM_BUFSIZ];
    snprintf(rate_prefix,BFAM_BUFSIZ,"%s%d_",BFAM_ADAMS_PREFIX,
        (ts->currentStage+1)%ts->nStages);
    BFAM_LDEBUG("Adams step: RK rate rate_prefix %s",rate_prefix);
    BFAM_ASSERT(ts->lsrk->step_extended);
    ts->lsrk->step_extended((bfam_ts_t*)ts->lsrk,dt,rate_prefix,"","");
    if(ts->numSteps+2 >= ts->nStages)
    {
      bfam_ts_lsrk_free(ts->lsrk);
      ts->lsrk = NULL;
    }
  }
  else
    /* q_{n+1}  := q_{n} + dt \sum_{k=0}^{m} a_{k} dq_{n-k} */
    bfam_dictionary_allprefixed_ptr(&ts->elems,
        "",&bfam_ts_adams_update,&data);

  /* shift the stage counter */
  ts->currentStage = (ts->currentStage+1)%ts->nStages;
  ts->numSteps++;

  /* update the stage time */
  ts->t += dt;

}

void
bfam_ts_adams_init(
    bfam_ts_adams_t*       ts,
    bfam_domain_t*         dom,
    bfam_ts_adams_method_t method,
    bfam_domain_match_t    subdom_match,
    const char**           subdom_tags,
    bfam_domain_match_t    comm_match,
    const char**           comm_tags,
    MPI_Comm               mpicomm,
    int                    mpitag,
    void*                  comm_data,
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
      const char *rate_prefix, const bfam_long_real_t a),
    const int RK_init)
{
  BFAM_LDEBUG("ADAMS INIT");

  /*
   * set up some preliminaries
   */
  bfam_ts_init(&ts->base, dom);
  bfam_dictionary_init(&ts->elems);
  ts->t  = BFAM_LONG_REAL(0.0);
  ts->base.step = &bfam_ts_adams_step;

  /*
   * store the function calls
   */
  ts->scale_rates = scale_rates;
  ts->intra_rhs   = intra_rhs;
  ts->inter_rhs   = inter_rhs;
  ts->add_rates   = add_rates;

  ts->currentStage = 0;
  ts->numSteps     = 0;

  ts->lsrk         = NULL;


  switch(method)
  {
    default:
      BFAM_WARNING("Invalid Adams scheme, using ADAMS_3");
    case BFAM_TS_ADAMS_3:
      ts->nStages = 3;
      /*
      ts->A = bfam_malloc_aligned(ts->nStages*sizeof(bfam_long_real_t));

      ts->A[0] = BFAM_LONG_REAL(23.0)/
                 BFAM_LONG_REAL(12.0);
      ts->A[1] = BFAM_LONG_REAL(-4.0)/
                 BFAM_LONG_REAL( 3.0);
      ts->A[2] = BFAM_LONG_REAL( 5.0)/
                 BFAM_LONG_REAL(12.0);
      */

      /* if necessary initialize the RK scheme */
      if(RK_init)
      {
        ts->lsrk = bfam_ts_lsrk_new_extended(dom, BFAM_TS_LSRK_KC54,
            subdom_match, subdom_tags, comm_match, comm_tags, mpicomm, mpitag,
            comm_data, aux_rates, scale_rates, intra_rhs, inter_rhs, add_rates,
            0);
      }

      break;
    case BFAM_TS_ADAMS_1:
      ts->nStages = 1;
      /*
      ts->A = bfam_malloc_aligned(ts->nStages*sizeof(bfam_long_real_t));

      ts->A[0] = BFAM_LONG_REAL(1.0);
      */

      /* if necessary initialize the RK scheme */
      if(RK_init)
      {
        ts->lsrk = bfam_ts_lsrk_new_extended(dom, BFAM_TS_LSRK_FE,
            subdom_match, subdom_tags, comm_match, comm_tags, mpicomm, mpitag,
            comm_data, aux_rates, scale_rates, intra_rhs, inter_rhs, add_rates,
            0);
      }

      break;
    case BFAM_TS_ADAMS_2:
      ts->nStages = 2;
      /*
      ts->A = bfam_malloc_aligned(ts->nStages*sizeof(bfam_long_real_t));

      ts->A[0] = BFAM_LONG_REAL( 3.0)/
                 BFAM_LONG_REAL( 2.0);
      ts->A[1] = BFAM_LONG_REAL(-1.0)/
                 BFAM_LONG_REAL( 2.0);
      */

      /* if necessary initialize the RK scheme */
      if(RK_init)
      {
        ts->lsrk = bfam_ts_lsrk_new_extended(dom, BFAM_TS_LSRK_HEUN,
            subdom_match, subdom_tags, comm_match, comm_tags, mpicomm, mpitag,
            comm_data, aux_rates, scale_rates, intra_rhs, inter_rhs, add_rates,
            0);
      }

      break;
    case BFAM_TS_ADAMS_4:
      ts->nStages = 4;
      /*
      ts->A = bfam_malloc_aligned(ts->nStages*sizeof(bfam_long_real_t));

      ts->A[0] = BFAM_LONG_REAL( 55.0)/
                 BFAM_LONG_REAL( 24.0);
      ts->A[1] = BFAM_LONG_REAL(-59.0)/
                 BFAM_LONG_REAL( 24.0);
      ts->A[2] = BFAM_LONG_REAL( 37.0)/
                 BFAM_LONG_REAL( 24.0);
      ts->A[3] = BFAM_LONG_REAL(  3.0)/
                 BFAM_LONG_REAL(  8.0);
      */

      /* if necessary initialize the RK scheme */
      if(RK_init)
      {
        ts->lsrk = bfam_ts_lsrk_new_extended(dom, BFAM_TS_LSRK_KC54,
            subdom_match, subdom_tags, comm_match, comm_tags, mpicomm, mpitag,
            comm_data, aux_rates, scale_rates, intra_rhs, inter_rhs, add_rates,
            0);
      }

      break;
  }

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
     BFAM_ABORT_IF_NOT(rval != 1, "Issue adding subdomain %s", subs[s]->name);

     for(int n = 0; n < ts->nStages; n++)
     {
       char aux_rates_name[BFAM_BUFSIZ];
       snprintf(aux_rates_name,BFAM_BUFSIZ,"%s%d_",BFAM_ADAMS_PREFIX,n);
       aux_rates(subs[s],aux_rates_name);
     }
   }

  /*
   * Set up the communicator we will use
   */
   ts->comm = bfam_communicator_new(dom,comm_match,comm_tags,mpicomm,mpitag,
       comm_data);
}

void
bfam_ts_adams_free(bfam_ts_adams_t* ts)
{
  BFAM_LDEBUG("ADAMS FREE");
  if(ts->lsrk != NULL)
    bfam_ts_lsrk_free(ts->lsrk);
  ts->lsrk = NULL;
  bfam_communicator_free(ts->comm);
  bfam_free(ts->comm);
  ts->comm = NULL;
  bfam_dictionary_clear(&ts->elems);
  /*
  bfam_free_aligned(ts->A);
  ts->A = NULL;
  */
  ts->nStages = 0;
  ts->t  = NAN;
  bfam_ts_free(&ts->base);
}

void
bfam_ts_adams_set_time(bfam_ts_adams_t* ts,bfam_long_real_t time)
{
  ts->t = time;
}

bfam_long_real_t
bfam_ts_adams_get_time(bfam_ts_adams_t* ts)
{
  return ts->t;
}
