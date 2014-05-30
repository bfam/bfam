#include <bfam_timestep_local_adams.h>
#include <bfam_log.h>

#define BFAM_LOCAL_ADAMS_PREFIX ("_local_adams_rate_")
#define BFAM_LOCAL_ADAMS_LVL_PREFIX ("_local_adams_lvl_")

void
bfam_ts_local_adams_fill_level_tag(char* tag, size_t buf_sz, int level)
{
  snprintf(tag,buf_sz,"%s%d",BFAM_LOCAL_ADAMS_LVL_PREFIX,level);
}

bfam_ts_local_adams_t*
bfam_ts_local_adams_new(bfam_domain_t* dom, bfam_ts_local_adams_method_t method,
    bfam_domain_match_t subdom_match, const char** subdom_tags,
    bfam_domain_match_t comm_match, const char** comm_tags,
    MPI_Comm mpicomm, int mpitag, void *comm_data,
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
  bfam_ts_local_adams_init(newTS, dom, method, subdom_match, subdom_tags,
      comm_match, comm_tags, mpicomm, mpitag, comm_data, aux_rates, glue_rates,
      scale_rates,intra_rhs,inter_rhs,add_rates,RK_init);
  return newTS;
}

void
bfam_ts_local_adams_init(
    bfam_ts_local_adams_t*       ts,
    bfam_domain_t*               dom,
    bfam_ts_local_adams_method_t method,
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
  // ts->base.step = &bfam_ts_local_adams_step;

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

  /*
   * Set up the communicator we will use
   */
   ts->comm = bfam_communicator_new(dom,comm_match,comm_tags,mpicomm,mpitag,
       comm_data);
}
