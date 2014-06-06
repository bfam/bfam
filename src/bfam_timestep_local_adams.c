#include <bfam_timestep_local_adams.h>
#include <bfam_log.h>

#define BFAM_LOCAL_ADAMS_PREFIX ("_local_adams_rate_")
#define BFAM_LOCAL_ADAMS_LVL_PREFIX "_local_adams_lvl_"
#define BFAM_LOCAL_ADAMS_COMM_LVL_PREFIX ("_local_adams_comm_lvl_")

typedef struct bfam_ts_local_allprefix
{
  bfam_ts_local_adams_t *ts;
  bfam_locidx_t    lvl; /* this marks the level to the updated */
  bfam_long_real_t dt;  /* this is the fastest time step */
} bfam_ts_local_adams_allprefix_t;

void
bfam_ts_local_adams_fill_level_tag(char* tag, size_t buf_sz, int level)
{
  snprintf(tag,buf_sz,"%s%d",BFAM_LOCAL_ADAMS_LVL_PREFIX,level);
}

static int
get_tag_level_number(const char * key, void *val)
{
  BFAM_ASSERT(*(bfam_locidx_t*)val < 0);
  return sscanf(key, BFAM_LOCAL_ADAMS_LVL_PREFIX"%"BFAM_LOCIDX_PRId,
      (bfam_locidx_t*) val);
}

void
bfam_ts_local_adams_fill_comm_level_tag(char* tag, size_t buf_sz, int level)
{
  snprintf(tag,buf_sz,"%s%d",BFAM_LOCAL_ADAMS_COMM_LVL_PREFIX,level);
}

bfam_ts_local_adams_t*
bfam_ts_local_adams_new(bfam_domain_t* dom, bfam_ts_local_adams_method_t method,
    bfam_locidx_t num_lvl, bfam_domain_match_t subdom_match,
    const char** subdom_tags, bfam_domain_match_t comm_match, const char**
    comm_tags, MPI_Comm mpicomm, int mpitag, void *comm_data,
    void (*aux_rates) (bfam_subdomain_t *thisSubdomain, const char *prefix),
    void (*glue_rates) (bfam_subdomain_t *thisSubdomain, const char *prefix),
    void (*scale_rates) (bfam_subdomain_t *thisSubdomain,
      const char *rate_prefix, const bfam_long_real_t a),
    void (*intra_rhs) (bfam_subdomain_t *thisSubdomain,
      const char *rate_prefix, const char *minus_rate_prefix,
      const char *field_prefix, const bfam_long_real_t t),
    void (*inter_rhs) (bfam_subdomain_t *thisSubdomain,
      const char *rate_prefix, const char *minus_rate_prefix,
      const char *field_prefix, const bfam_long_real_t t),
    void (*add_rates) (bfam_subdomain_t *thisSubdomain,
      const char *field_prefix_lhs, const char *field_prefix_rhs,
      const char *rate_prefix, const bfam_long_real_t a),
    const int RK_init)
{
  bfam_ts_local_adams_t* newTS = bfam_malloc(sizeof(bfam_ts_local_adams_t));
  bfam_ts_local_adams_init(newTS, dom, method, num_lvl, subdom_match,
      subdom_tags, comm_match, comm_tags, mpicomm, mpitag, comm_data, aux_rates,
      glue_rates, scale_rates,intra_rhs,inter_rhs,add_rates,RK_init);
  return newTS;
}

static int
bfam_ts_local_adams_intra_rhs(const char * key, void *val, void *arg)
{
  bfam_ts_local_adams_allprefix_t *data =
    (bfam_ts_local_adams_allprefix_t *) arg;
  bfam_subdomain_t* sub = (bfam_subdomain_t*) val;

  /* We have to determine if this domain is updated or not */
  bfam_locidx_t lvl = -1;
  bfam_critbit0_allprefixed(&sub->tags, BFAM_LOCAL_ADAMS_LVL_PREFIX,
      get_tag_level_number,&lvl);
  char *rate_prefix = NULL;
  char rate_prefix_storage[BFAM_BUFSIZ];
  if(lvl > -1 && lvl <= data->lvl)
  {
    rate_prefix = rate_prefix_storage;
    snprintf(rate_prefix_storage,BFAM_BUFSIZ,"%s%d_",BFAM_LOCAL_ADAMS_PREFIX,
        data->ts->currentStageArray[lvl]%data->ts->nStages);
    BFAM_LDEBUG("Local Adams intra: level %"BFAM_LOCIDX_PRId
        " using rate prefix %s",lvl,rate_prefix_storage);
  }

  bfam_locidx_t m_lvl = -1;
  if(sub->glue_m && sub->glue_m->sub_m)
    bfam_critbit0_allprefixed(&sub->glue_m->sub_m->tags,
        BFAM_LOCAL_ADAMS_LVL_PREFIX, get_tag_level_number,&m_lvl);
  char *minus_rate_prefix = NULL;
  char minus_rate_prefix_storage[BFAM_BUFSIZ];
  if(m_lvl > -1 && m_lvl <= data->lvl)
  {
    minus_rate_prefix = minus_rate_prefix_storage;
    snprintf(minus_rate_prefix_storage,BFAM_BUFSIZ,"%s%d_",BFAM_LOCAL_ADAMS_PREFIX,
        data->ts->currentStageArray[m_lvl]%data->ts->nStages);
    BFAM_LDEBUG("Local Adams intra: level %"BFAM_LOCIDX_PRId
        " using minus rate prefix %s",m_lvl,minus_rate_prefix_storage);
  }

  if(rate_prefix)
    data->ts->scale_rates(sub, rate_prefix, 0);
  if(rate_prefix || minus_rate_prefix)
    data->ts->intra_rhs(sub, rate_prefix, minus_rate_prefix, "", data->ts->t);

  return 1;
}

static int
bfam_ts_local_adams_inter_rhs(const char * key, void *val, void *arg)
{
  bfam_ts_local_adams_allprefix_t *data =
    (bfam_ts_local_adams_allprefix_t *) arg;
  bfam_subdomain_t* sub = (bfam_subdomain_t*) val;

  /* We have to determine if this domain is updated or not */
  bfam_locidx_t lvl = -1;
  bfam_critbit0_allprefixed(&sub->tags, BFAM_LOCAL_ADAMS_LVL_PREFIX,
      get_tag_level_number,&lvl);
  char *rate_prefix = NULL;
  char rate_prefix_storage[BFAM_BUFSIZ];
  if(lvl > -1 && lvl <= data->lvl)
  {
    rate_prefix = rate_prefix_storage;
    snprintf(rate_prefix_storage,BFAM_BUFSIZ,"%s%d_",BFAM_LOCAL_ADAMS_PREFIX,
        data->ts->currentStageArray[lvl]%data->ts->nStages);
    BFAM_LDEBUG("Local Adams inter: level %"BFAM_LOCIDX_PRId
        " using rate prefix %s",lvl,rate_prefix_storage);
  }

  bfam_locidx_t m_lvl = -1;
  if(sub->glue_m && sub->glue_m->sub_m)
    bfam_critbit0_allprefixed(&sub->glue_m->sub_m->tags,
        BFAM_LOCAL_ADAMS_LVL_PREFIX, get_tag_level_number,&m_lvl);
  char *minus_rate_prefix = NULL;
  char minus_rate_prefix_storage[BFAM_BUFSIZ];
  if(m_lvl > -1 && m_lvl <= data->lvl)
  {
    minus_rate_prefix = minus_rate_prefix_storage;
    snprintf(minus_rate_prefix_storage,BFAM_BUFSIZ,"%s%d_",BFAM_LOCAL_ADAMS_PREFIX,
        data->ts->currentStageArray[m_lvl]%data->ts->nStages);
    BFAM_LDEBUG("Local Adams inter: level %"BFAM_LOCIDX_PRId
        " using minus rate prefix %s",m_lvl,minus_rate_prefix_storage);
  }

  if(rate_prefix || minus_rate_prefix)
    data->ts->inter_rhs(sub, rate_prefix, minus_rate_prefix, "", data->ts->t);

  return 1;
}

static void
bfam_ts_local_adams_step(bfam_ts_t *a_ts, bfam_long_real_t dt)
{
  bfam_ts_local_adams_t* ts = (bfam_ts_local_adams_t*) a_ts;
  bfam_locidx_t num_steps = 1<<(ts->numLevels-1);
  BFAM_LDEBUG("Number of steps for the local time stepper %"BFAM_LOCIDX_PRId,
      num_steps);

  dt /= ts->numLevels;
  bfam_ts_local_adams_allprefix_t data;
  data.ts = ts;
  data.dt = dt;
  data.lvl = -1;

  /* determine the level of comm to do: max of this update and last update */
  bfam_locidx_t last_lvl = 0;
  for(bfam_locidx_t s = 0; s < num_steps; s++)
  {
    BFAM_LDEBUG("local time step number %"BFAM_LOCIDX_PRId, s);
    data.lvl = 0;

    /* loop through the levels */
    for(bfam_locidx_t lvl = 0; lvl < ts->numLevels; lvl++)
    {
      bfam_locidx_t chk = 1 << lvl;
      if(!(s%chk))
      {
        BFAM_LDEBUG("level %"BFAM_LOCIDX_PRId" to be updated",lvl);
        data.lvl = BFAM_MAX(data.lvl, lvl);
      }
    }
    bfam_locidx_t comm_lvl = BFAM_MAX(data.lvl, last_lvl);
    BFAM_LDEBUG("step %"BFAM_LOCIDX_PRId
        ": update level %"BFAM_LOCIDX_PRId
        ": and communication level %"BFAM_LOCIDX_PRId,
        s, data.lvl, comm_lvl);

    /* start the communication */
    bfam_communicator_start(ts->comm_array[comm_lvl-1]);

    /* Do the intra work for the levels to be udpated */
    bfam_dictionary_allprefixed_ptr(&ts->elems,
        "",&bfam_ts_local_adams_intra_rhs,&data);

    /* finish the communication */
    bfam_communicator_finish(ts->comm_array[comm_lvl-1]);

    /* Do the inter work for the levels to be udpated */
    bfam_dictionary_allprefixed_ptr(&ts->elems,
        "",&bfam_ts_local_adams_inter_rhs,&data);

    /* set last update to next update */
    last_lvl = data.lvl;

    /* update the stage counters */
    for(bfam_locidx_t s=0; s <= data.lvl;s++)
      ts->currentStageArray[s] = (ts->currentStageArray[s]+1)%ts->nStages;
  }
}

void
bfam_ts_local_adams_init(
    bfam_ts_local_adams_t*       ts,
    bfam_domain_t*               dom,
    bfam_ts_local_adams_method_t method,
    bfam_locidx_t                num_lvl,
    bfam_domain_match_t          subdom_match,
    const char**                 subdom_tags,
    bfam_domain_match_t          comm_match,
    const char**                 comm_tags,
    MPI_Comm                     mpicomm,
    int                          mpitag,
    void*                        comm_data,
    void (*aux_rates) (bfam_subdomain_t *thisSubdomain, const char *prefix),
    void (*glue_rates) (bfam_subdomain_t *thisSubdomain, const char *prefix),
    void (*scale_rates) (bfam_subdomain_t *thisSubdomain,
      const char *rate_prefix, const bfam_long_real_t a),
    void (*intra_rhs) (bfam_subdomain_t *thisSubdomain,
      const char *rate_prefix, const char *minus_rate_prefix,
      const char *field_prefix, const bfam_long_real_t t),
    void (*inter_rhs) (bfam_subdomain_t *thisSubdomain,
      const char *rate_prefix, const char *minus_rate_prefix,
      const char *field_prefix, const bfam_long_real_t t),
    void (*add_rates) (bfam_subdomain_t *thisSubdomain,
      const char *field_prefix_lhs, const char *field_prefix_rhs,
      const char *rate_prefix, const bfam_long_real_t a),
    const int RK_init)
{
  BFAM_LDEBUG("LOCAL ADAMS INIT");

  /*
   * set up some preliminaries
   */
  bfam_ts_init(&ts->base, dom);
  bfam_dictionary_init(&ts->elems);
  ts->t  = BFAM_LONG_REAL(0.0);
  ts->base.step = &bfam_ts_local_adams_step;

  /*
   * store the function calls
   */
  ts->scale_rates = scale_rates;
  ts->intra_rhs   = intra_rhs;
  ts->inter_rhs   = inter_rhs;
  ts->add_rates   = add_rates;

  ts->numSteps     = 0;

  /*
   * fast log2 computation for ints using bitwise operations from
   * http://stackoverflow.com/questions/994593/how-to-do-an-integer-log2-in-c
   */
  ts->numLevels = num_lvl;

  ts->currentStageArray = bfam_malloc(ts->numLevels*sizeof(bfam_locidx_t));
  for(bfam_locidx_t k = 0; k < ts->numLevels; k++)
    ts->currentStageArray[k] = 0;

  ts->lsrk         = NULL;
  ts->comm_array   = NULL;

  switch(method)
  {
    default:
      BFAM_WARNING("Invalid Adams scheme, using ADAMS_3");
    case BFAM_TS_LOCAL_ADAMS_3:
      ts->nStages = 3;

      /* if necessary initialize the RK scheme */
      if(RK_init)
      {
        ts->lsrk = bfam_ts_lsrk_new_extended(dom, BFAM_TS_LSRK_KC54,
            subdom_match, subdom_tags, comm_match, comm_tags, mpicomm, mpitag,
            comm_data, aux_rates, scale_rates, intra_rhs, inter_rhs, add_rates,
            0);
      }

      break;
    case BFAM_TS_LOCAL_ADAMS_1:
      ts->nStages = 1;

      /* if necessary initialize the RK scheme */
      if(RK_init)
      {
        ts->lsrk = bfam_ts_lsrk_new_extended(dom, BFAM_TS_LSRK_FE,
            subdom_match, subdom_tags, comm_match, comm_tags, mpicomm, mpitag,
            comm_data, aux_rates, scale_rates, intra_rhs, inter_rhs, add_rates,
            0);
      }

      break;
    case BFAM_TS_LOCAL_ADAMS_2:
      ts->nStages = 2;

      /* if necessary initialize the RK scheme */
      if(RK_init)
      {
        ts->lsrk = bfam_ts_lsrk_new_extended(dom, BFAM_TS_LSRK_HEUN,
            subdom_match, subdom_tags, comm_match, comm_tags, mpicomm, mpitag,
            comm_data, aux_rates, scale_rates, intra_rhs, inter_rhs, add_rates,
            0);
      }

      break;
    case BFAM_TS_LOCAL_ADAMS_4:
      ts->nStages = 4;

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
      snprintf(aux_rates_name,BFAM_BUFSIZ,"%s%d_",BFAM_LOCAL_ADAMS_PREFIX,n);
      aux_rates(subs[s],aux_rates_name);
      glue_rates(subs[s],aux_rates_name);
    }
  }

  numSubs = 0;
  bfam_domain_get_subdomains(dom,comm_match,comm_tags,
      dom->numSubdomains,subs,&numSubs);

  /* this tracks the max level that I know about */
#ifdef BFAM_DEBUG
  int local_max_levels = 0;
#endif

  /* find the plus and minus levels for the glue grids */
  for(int s = 0; s < numSubs;s++)
  {
    bfam_locidx_t m_lvl = -1;
    bfam_locidx_t p_lvl = -1;

    BFAM_ASSERT(subs[s]->glue_m);
    BFAM_ASSERT(subs[s]->glue_p);

    bfam_critbit0_allprefixed(&subs[s]->glue_p->tags,
        BFAM_LOCAL_ADAMS_LVL_PREFIX, get_tag_level_number,&p_lvl);
    bfam_critbit0_allprefixed(&subs[s]->glue_m->tags,
        BFAM_LOCAL_ADAMS_LVL_PREFIX, get_tag_level_number,&m_lvl);

    BFAM_ASSERT(m_lvl >= 0);
    BFAM_ASSERT(p_lvl >= 0);

#ifdef BFAM_DEBUG
    /* just a sanity check since these should match */
    BFAM_ASSERT(subs[s]->glue_m->sub_m);
    bfam_locidx_t s_lvl = -1;
    bfam_critbit0_allprefixed(&subs[s]->glue_m->sub_m->tags,
        BFAM_LOCAL_ADAMS_LVL_PREFIX, get_tag_level_number,&s_lvl);
    BFAM_ASSERT(s_lvl == m_lvl);
#endif

    /* communicate as infrequently as needed, so use the higher level */
    char comm_lvl_tag[BFAM_BUFSIZ];
    bfam_ts_local_adams_fill_comm_level_tag(comm_lvl_tag,BFAM_BUFSIZ,
        BFAM_MAX(m_lvl, p_lvl));
    bfam_subdomain_add_tag(subs[s],comm_lvl_tag);

#ifdef BFAM_DEBUG
    /* update the max level */
    local_max_levels = BFAM_MAX(local_max_levels,BFAM_MAX(m_lvl,p_lvl));
#endif
  }


#ifdef BFAM_DEBUG
  BFAM_ASSERT(local_max_levels);
  BFAM_ASSERT(local_max_levels <= num_lvl);
#endif

  /* loop through all possible commmunication tags */
  char *local_comm_tags[ts->numLevels+1];
  char tag_stor[BFAM_BUFSIZ*ts->numLevels];
  ts->comm_array = bfam_malloc(ts->numLevels*sizeof(bfam_ts_local_adams_t*));

  for(int lvl = 1, k=0; k < ts->numLevels; lvl*=2, k++)
  {
    local_comm_tags[k] = &tag_stor[BFAM_BUFSIZ*k];
    bfam_ts_local_adams_fill_comm_level_tag(&tag_stor[BFAM_BUFSIZ*k],
        BFAM_BUFSIZ, lvl);
    local_comm_tags[k+1] = NULL;

    /*
     * Set up the communicator we will use
     */
    ts->comm_array[k] = bfam_communicator_new(dom, BFAM_DOMAIN_OR,
        (const char**)local_comm_tags, mpicomm, mpitag, comm_data);
  }
}

void
bfam_ts_local_adams_free(bfam_ts_local_adams_t* ts)
{
  BFAM_LDEBUG("LOCAL ADAMS FREE");
  if(ts->lsrk != NULL)
  {
    bfam_ts_lsrk_free(ts->lsrk);
    bfam_free(ts->lsrk);
  }
  ts->lsrk = NULL;
  for(int k = 0; k < ts->numLevels; k++)
  {
    bfam_communicator_free(ts->comm_array[k]);
    bfam_free(ts->comm_array[k]);
    ts->comm_array[k] = NULL;
  }
  bfam_free(ts->comm_array);
  bfam_dictionary_clear(&ts->elems);
  /*
  bfam_free_aligned(ts->A);
  ts->A = NULL;
  */
  ts->nStages = 0;
  bfam_free(ts->currentStageArray);
  ts->t  = NAN;
  bfam_ts_free(&ts->base);
}
