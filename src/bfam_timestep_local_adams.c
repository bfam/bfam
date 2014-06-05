#include <bfam_timestep_local_adams.h>
#include <bfam_log.h>

#define BFAM_LOCAL_ADAMS_PREFIX ("_local_adams_rate_")
#define BFAM_LOCAL_ADAMS_LVL_PREFIX "_local_adams_lvl_"
#define BFAM_LOCAL_ADAMS_COMM_LVL_PREFIX ("_local_adams_comm_lvl_")

void
bfam_ts_local_adams_fill_level_tag(char* tag, size_t buf_sz, int level)
{
  snprintf(tag,buf_sz,"%s%d",BFAM_LOCAL_ADAMS_LVL_PREFIX,level);
}

static int
get_tag_level_number(const char * key, void *val)
{
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
    bfam_locidx_t max_level, bfam_domain_match_t subdom_match,
    const char** subdom_tags, bfam_domain_match_t comm_match, const char**
    comm_tags, MPI_Comm mpicomm, int mpitag, void *comm_data,
    void (*aux_rates) (bfam_subdomain_t *thisSubdomain, const char *prefix),
    void (*glue_rates) (bfam_subdomain_t *thisSubdomain, const char *prefix),
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
  bfam_ts_local_adams_t* newTS = bfam_malloc(sizeof(bfam_ts_local_adams_t));
  bfam_ts_local_adams_init(newTS, dom, method, max_level, subdom_match,
      subdom_tags, comm_match, comm_tags, mpicomm, mpitag, comm_data, aux_rates,
      glue_rates, scale_rates,intra_rhs,inter_rhs,add_rates,RK_init);
  return newTS;
}

static void
do_local_step(bfam_ts_local_adams_t* ts, bfam_long_real_t dt,
    bfam_locidx_t lvl)
{
  /* First we tell the next level down to do two steps */
  if(lvl > 0)
  {
    do_local_step(ts, BFAM_LONG_REAL(0.5)*dt, lvl-1);
    do_local_step(ts, BFAM_LONG_REAL(0.5)*dt, lvl-1);
  }

  /* And now we can do a step */

  /*
   * start the communication
   */

  /*
   * do the intra work
   */

  /*
   * finish the communication
   */

  /*
   * do the inter work
   */
}

static void
bfam_ts_local_adams_step(bfam_ts_t *a_ts, bfam_long_real_t dt)
{
  bfam_ts_local_adams_t* ts = (bfam_ts_local_adams_t*) a_ts;
  do_local_step(ts, dt, ts->numLevels);
}

void
bfam_ts_local_adams_init(
    bfam_ts_local_adams_t*       ts,
    bfam_domain_t*               dom,
    bfam_ts_local_adams_method_t method,
    bfam_locidx_t                max_level,
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
  int tmp = max_level;
  ts->numLevels = 1;
  while (tmp >>= 1) ++ts->numLevels;

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

    BFAM_ASSERT(m_lvl > 0);
    BFAM_ASSERT(p_lvl > 0);

#ifdef BFAM_DEBUG
    /* just a sanity check since these should match */
    BFAM_ASSERT(subs[s]->glue_m->sub_m);
    bfam_locidx_t s_lvl = -1;
    bfam_critbit0_allprefixed(&subs[s]->glue_m->sub_m->tags,
        BFAM_LOCAL_ADAMS_LVL_PREFIX, get_tag_level_number,&s_lvl);
    BFAM_ASSERT(s_lvl == m_lvl);
#endif

    /* communicate as infrequently as needed, so use the higher level */
    char lvl_tag[BFAM_BUFSIZ];
    bfam_ts_local_adams_fill_comm_level_tag(lvl_tag,BFAM_BUFSIZ,
        BFAM_MAX(m_lvl, p_lvl));
    bfam_subdomain_add_tag(subs[s],lvl_tag);

#ifdef BFAM_DEBUG
    /* update the max level */
    local_max_levels = BFAM_MAX(local_max_levels,BFAM_MAX(m_lvl,p_lvl));
#endif
  }


#ifdef BFAM_DEBUG
  BFAM_ASSERT(local_max_levels);
  BFAM_ASSERT(local_max_levels <= max_level);
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
