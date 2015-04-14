#include <bfam_timestep_local_adams.h>
#include <bfam_log.h>

#define BFAM_LOCAL_ADAMS_PREFIX ("_local_adams_rate_")
#define BFAM_LOCAL_ADAMS_LVL_PREFIX "_local_adams_lvl_"
#define BFAM_LOCAL_ADAMS_COMM_LVL_PREFIX ("_local_adams_comm_lvl_")

/* #define BFAM_LOCAL_ADAMS_ALWAYS_INTERP */

typedef struct bfam_ts_local_allprefix
{
  bfam_ts_local_adams_t *ts;
  bfam_locidx_t step;  /* what number step are we on */
  bfam_locidx_t lvl;   /* this marks the level to the updated */
  bfam_long_real_t dt; /* this is the fastest time step */
} bfam_ts_local_adams_allprefix_t;

void bfam_ts_local_adams_fill_level_tag(char *tag, size_t buf_sz, int level)
{
  snprintf(tag, buf_sz, "%s%d", BFAM_LOCAL_ADAMS_LVL_PREFIX, level);
}

static int get_tag_level_number(const char *key, void *val)
{
  BFAM_ASSERT(*(bfam_locidx_t *)val < 0);
  return sscanf(key, BFAM_LOCAL_ADAMS_LVL_PREFIX "%" BFAM_LOCIDX_PRId,
                (bfam_locidx_t *)val);
}

void bfam_ts_local_adams_fill_comm_level_tag(char *tag, size_t buf_sz,
                                             int level)
{
  snprintf(tag, buf_sz, "%s%d", BFAM_LOCAL_ADAMS_COMM_LVL_PREFIX, level);
}

static void comm_send_prefix(bfam_subdomain_t *sub, char *prefix,
                             size_t buf_siz, void *user_data)
{
  BFAM_ASSERT(sub->glue_m);
  BFAM_ASSERT(sub->glue_p);

  /* Get the plus and minus side levels */
  bfam_locidx_t m_lvl = -1;
  bfam_critbit0_allprefixed(&sub->glue_m->tags, BFAM_LOCAL_ADAMS_LVL_PREFIX,
                            get_tag_level_number, &m_lvl);
  BFAM_ASSERT(m_lvl >= 0);

  bfam_locidx_t p_lvl = -1;
  bfam_critbit0_allprefixed(&sub->glue_p->tags, BFAM_LOCAL_ADAMS_LVL_PREFIX,
                            get_tag_level_number, &p_lvl);
  BFAM_ASSERT(p_lvl >= 0);

  bfam_ts_local_adams_allprefix_t *data =
      (bfam_ts_local_adams_allprefix_t *)user_data;
  BFAM_ASSERT(data);

/*
 * send the rates if my level number is greater than my neighbors
 * and I am not being not being updated (i.e., my level is greater
 * than the updated level number)
 *
 * We have to shift back on b/c the stage number has already been moved
 * forward by 1
 */
#ifdef BFAM_LOCAL_ADAMS_ALWAYS_INTERP
  if (data->ts->numStepsArray[m_lvl] > 0)
#else
  if (p_lvl < m_lvl && data->lvl < m_lvl)
#endif
  {
    snprintf(prefix, buf_siz, "%s%d_", BFAM_LOCAL_ADAMS_PREFIX,
             (data->ts->currentStageArray[m_lvl] + data->ts->nStages - 1) %
                 data->ts->nStages);
    BFAM_LDEBUG("send: %d %d %d %s", data->lvl, m_lvl, p_lvl, prefix);
  }
}

static void comm_recv_prefix(bfam_subdomain_t *sub, char *prefix,
                             size_t buf_siz, void *user_data)
{
  BFAM_ASSERT(sub->glue_m);
  BFAM_ASSERT(sub->glue_p);

  /* Get the plus and minus side levels */
  bfam_locidx_t m_lvl = -1;
  bfam_critbit0_allprefixed(&sub->glue_m->tags, BFAM_LOCAL_ADAMS_LVL_PREFIX,
                            get_tag_level_number, &m_lvl);
  BFAM_ASSERT(m_lvl >= 0);

  bfam_locidx_t p_lvl = -1;
  bfam_critbit0_allprefixed(&sub->glue_p->tags, BFAM_LOCAL_ADAMS_LVL_PREFIX,
                            get_tag_level_number, &p_lvl);
  BFAM_ASSERT(p_lvl >= 0);

  bfam_ts_local_adams_allprefix_t *data =
      (bfam_ts_local_adams_allprefix_t *)user_data;
  BFAM_ASSERT(data);

/*
 * recv the rates if neigh level number is greater than my level
 * and neighbor is not being not being updated (i.e., neighbors level is
 * greater than the updated level number)
 *
 * We have to shift back on b/c the stage number has already been moved
 * forward by 1
 */
#ifdef BFAM_LOCAL_ADAMS_ALWAYS_INTERP
  if (data->ts->numStepsArray[m_lvl] > 0)
#else
  if (m_lvl < p_lvl && data->lvl < p_lvl)
#endif
  {
    snprintf(prefix, buf_siz, "%s%d_", BFAM_LOCAL_ADAMS_PREFIX,
             (data->ts->currentStageArray[p_lvl] + data->ts->nStages - 1) %
                 data->ts->nStages);
    BFAM_LDEBUG("recv: %d %d %d %s", data->lvl, m_lvl, p_lvl, prefix);
  }
}

bfam_ts_local_adams_t *bfam_ts_local_adams_new(
    bfam_domain_t *dom, bfam_ts_local_adams_method_t method,
    bfam_locidx_t num_lvl, bfam_domain_match_t subdom_match,
    const char **subdom_tags, bfam_domain_match_t comm_match,
    const char **comm_tags, MPI_Comm mpicomm, int mpitag,
    bfam_subdomain_comm_args_t *comm_data,
    void (*aux_rates)(bfam_subdomain_t *thisSubdomain, const char *prefix),
    void (*glue_rates)(bfam_subdomain_t *thisSubdomain, const char *prefix),
    void (*scale_rates)(bfam_subdomain_t *thisSubdomain,
                        const char *rate_prefix, const bfam_long_real_t a),
    void (*intra_rhs)(bfam_subdomain_t *thisSubdomain, const char *rate_prefix,
                      const char *minus_rate_prefix, const char *field_prefix,
                      const bfam_long_real_t t),
    void (*inter_rhs)(bfam_subdomain_t *thisSubdomain, const char *rate_prefix,
                      const char *minus_rate_prefix, const char *field_prefix,
                      const bfam_long_real_t t),
    void (*add_rates)(bfam_subdomain_t *thisSubdomain,
                      const char *field_prefix_lhs,
                      const char *field_prefix_rhs, const char *rate_prefix,
                      const bfam_long_real_t a),
    void (*add_rates_glue_p)(bfam_subdomain_t *thisSubdomain,
                             const char *field_prefix_lhs,
                             const char *field_prefix_rhs,
                             const char *rate_prefix, const bfam_long_real_t a),
    const int RK_init)
{
  bfam_ts_local_adams_t *newTS = bfam_malloc(sizeof(bfam_ts_local_adams_t));
  bfam_ts_local_adams_init(
      newTS, dom, method, num_lvl, subdom_match, subdom_tags, comm_match,
      comm_tags, mpicomm, mpitag, comm_data, aux_rates, glue_rates, scale_rates,
      intra_rhs, inter_rhs, add_rates, add_rates_glue_p, RK_init);
  return newTS;
}

static int bfam_ts_local_adams_intra_rhs(const char *key, void *val, void *arg)
{
  bfam_ts_local_adams_allprefix_t *data =
      (bfam_ts_local_adams_allprefix_t *)arg;
  bfam_subdomain_t *sub = (bfam_subdomain_t *)val;

  /* We have to determine if this domain is updated or not */
  bfam_locidx_t lvl = -1;
  bfam_critbit0_allprefixed(&sub->tags, BFAM_LOCAL_ADAMS_LVL_PREFIX,
                            get_tag_level_number, &lvl);
  char *rate_prefix = NULL;
  char rate_prefix_storage[BFAM_BUFSIZ];
  if (lvl > -1 && lvl <= data->lvl)
  {
    rate_prefix = rate_prefix_storage;
    snprintf(rate_prefix_storage, BFAM_BUFSIZ, "%s%d_", BFAM_LOCAL_ADAMS_PREFIX,
             data->ts->currentStageArray[lvl] % data->ts->nStages);
    BFAM_LDEBUG("Local Adams intra: level %" BFAM_LOCIDX_PRId
                " using rate prefix %s",
                lvl, rate_prefix_storage);
  }

  bfam_locidx_t m_lvl = -1;
  if (sub->glue_m && sub->glue_m->sub_m)
    bfam_critbit0_allprefixed(&sub->glue_m->sub_m->tags,
                              BFAM_LOCAL_ADAMS_LVL_PREFIX, get_tag_level_number,
                              &m_lvl);
  char *minus_rate_prefix = NULL;
  char minus_rate_prefix_storage[BFAM_BUFSIZ];
#ifdef BFAM_LOCAL_ADAMS_ALWAYS_INTERP
  if (data->ts->numStepsArray[m_lvl] > 0)
#else
  if (m_lvl > -1 && m_lvl <= data->lvl)
#endif
  {
    minus_rate_prefix = minus_rate_prefix_storage;
    snprintf(minus_rate_prefix_storage, BFAM_BUFSIZ, "%s%d_",
             BFAM_LOCAL_ADAMS_PREFIX,
             data->ts->currentStageArray[m_lvl] % data->ts->nStages);
    BFAM_LDEBUG("Local Adams intra: level %" BFAM_LOCIDX_PRId
                " using minus rate prefix %s",
                m_lvl, minus_rate_prefix_storage);
  }

  if (rate_prefix)
    data->ts->scale_rates(sub, rate_prefix, 0);
  if (rate_prefix || minus_rate_prefix)
    data->ts->intra_rhs(sub, rate_prefix, minus_rate_prefix, "", data->ts->t);

  return 1;
}

static void interp_fields(bfam_ts_local_adams_allprefix_t *data,
                          bfam_subdomain_t *sub, bfam_locidx_t m_lvl,
                          bfam_locidx_t p_lvl)
{
  /*
   * First we we interp, then we don't have our current fields so first we get
   * those
   */
  bfam_subdomain_comm_args_t *sub_comm_data =
      (bfam_subdomain_comm_args_t *)data->ts->comm_array[m_lvl]->user_args;
  BFAM_ASSERT(sub_comm_data);
  BFAM_ASSERT(sub_comm_data->user_data == NULL);
  BFAM_ASSERT(sub_comm_data->user_prefix_function == NULL);
  sub->glue_put_send_buffer(sub, NULL, 0, sub_comm_data);
#ifndef BFAM_LOCAL_ADAMS_ALWAYS_INTERP
  BFAM_ASSERT(m_lvl < p_lvl && p_lvl > data->lvl);
#endif

  /* last time the plus and minu sides were updated */
  bfam_locidx_t m_last = (1 << m_lvl) * (data->step / (1 << m_lvl) - 1);
  bfam_locidx_t p_last = (1 << p_lvl) * (data->step / (1 << p_lvl));

  /* the amount the plus field is ahead of actual field */
  bfam_locidx_t p_delta = m_last - p_last;

  /*
   * When doing the interpolation we assume that the guy whose rates we are
   * updating have steps that are size 1 apart, thus a is the fraction of the
   * step he has already done and b is the fraction of the step we are about to
   * do, that is, we need to add the rate
   *   \int_{a}^{b} f(q,t) dt
   * to the current field values
   */
  bfam_long_real_t a =
      (bfam_long_real_t)p_delta / (bfam_long_real_t)(1 << p_lvl);

  bfam_long_real_t b = (bfam_long_real_t)(p_delta + (1 << m_lvl)) /
                       (bfam_long_real_t)(1 << p_lvl);

  BFAM_LDEBUG("lvl %2d :: ngh lvl %2d"
              " :: last  %2d :: ngh last  %2d :: ngh delta %d"
              " :: a = %" BFAM_LONG_REAL_PRIe " :: b = %" BFAM_LONG_REAL_PRIe,
              m_lvl, p_lvl, m_last, p_last, p_delta, a, b);

  /*
   * we have to look at the p_lvl because these are the guys we are
   * interpolating
   */
  bfam_long_real_t A[4];
  bfam_locidx_t num_stages =
      BFAM_MIN(data->ts->numStepsArray[p_lvl], data->ts->nStages);
  switch (num_stages)
  {
  case 1:
  {
    /* \int_{a}^{b} 1 dt */
    A[0] = b - a;
  }
  break;
  case 2:
  {
    /* \int_{a}^{b} (t+1) dt */
    A[0] = (-a * a - 2 * a + b * b + 2 * b) / 2;

    /* \int_{a}^{b} -t dt */
    A[1] = (a * a - b * b) / 2;
  }
  break;
  case 3:
  {
    /* \int_{a}^{b} (t+1)(t+2)/2 dt */
    A[0] = (-2 * a * a * a - 9 * a * a - 12 * a + 2 * b * b * b + 9 * b * b +
            12 * b) /
           12;

    /* \int_{a}^{b} t(t+2)/(-1) dt */
    A[1] = (a * a * a + 3 * a * a - b * b * b - 3 * b * b) / 3;

    /* \int_{a}^{b} t(t+1)/2 dt */
    A[2] = (-2 * a * a * a - 3 * a * a + 2 * b * b * b + 3 * b * b) / 12;
    break;
  }
  case 4:
  {
    /* \int_{a}^{b} (t+1)(t+2)(t+3)/6 dt */
    A[0] = (-a * a * a * a - 8 * a * a * a - 22 * a * a - 24 * a +
            b * b * b * b + 8 * b * b * b + 22 * b * b + 24 * b) /
           24;

    /* \int_{a}^{b} t(t+2)(t+3)/(-2) dt */
    A[1] = (3 * a * a * a * a + 20 * a * a * a + 36 * a * a -
            3 * b * b * b * b - 20 * b * b * b - 36 * b * b) /
           24;

    /* \int_{a}^{b} t(t+1)(t+3)/(2) dt */
    A[2] = (-3 * a * a * a * a - 16 * a * a * a - 18 * a * a +
            3 * b * b * b * b + 16 * b * b * b + 18 * b * b) /
           24;

    /* \int_{a}^{b} t(t+1)(t+2)/(-6) dt */
    A[3] = (a * a * a * a + 4 * a * a * a + 4 * a * a - b * b * b * b -
            4 * b * b * b - 4 * b * b) /
           24;
  }
  break;
  default:
    BFAM_ABORT("Adams-Bashforth order %d not implemented", num_stages);
  }

  bfam_long_real_t dt = data->dt * (1 << m_lvl);
  for (int k = 0; k < num_stages; k++)
  {
    char prefix[BFAM_BUFSIZ];
    /* minus 1 here b/c we need the last guys (not the current guys */
    snprintf(prefix, BFAM_BUFSIZ, "%s%d_", BFAM_LOCAL_ADAMS_PREFIX,
             (data->ts->currentStageArray[p_lvl] + data->ts->nStages - k - 1) %
                 data->ts->nStages);
    BFAM_LDEBUG("Adams step: stage %d of %d using prefix %s", k, num_stages,
                prefix);
    data->ts->add_rates_glue_p(sub, "", "", prefix, dt * A[k]);
  }
}

static int bfam_ts_local_adams_inter_rhs(const char *key, void *val, void *arg)
{
  bfam_ts_local_adams_allprefix_t *data =
      (bfam_ts_local_adams_allprefix_t *)arg;
  bfam_subdomain_t *sub = (bfam_subdomain_t *)val;

  /* We have to determine if this domain is updated or not */
  bfam_locidx_t lvl = -1;
  bfam_critbit0_allprefixed(&sub->tags, BFAM_LOCAL_ADAMS_LVL_PREFIX,
                            get_tag_level_number, &lvl);
  char *rate_prefix = NULL;
  char rate_prefix_storage[BFAM_BUFSIZ];
  if (lvl > -1 && lvl <= data->lvl)
  {
    rate_prefix = rate_prefix_storage;
    snprintf(rate_prefix_storage, BFAM_BUFSIZ, "%s%d_", BFAM_LOCAL_ADAMS_PREFIX,
             data->ts->currentStageArray[lvl] % data->ts->nStages);
    BFAM_LDEBUG("Local Adams inter: level %" BFAM_LOCIDX_PRId
                " using rate prefix %s",
                lvl, rate_prefix_storage);
  }

  /*
   * Determine if the minus side exists and whether imterpolation need to occure
   * on the plus side
   */
  bfam_locidx_t m_lvl = -1;
  if (sub->glue_m && sub->glue_m->sub_m)
    bfam_critbit0_allprefixed(&sub->glue_m->sub_m->tags,
                              BFAM_LOCAL_ADAMS_LVL_PREFIX, get_tag_level_number,
                              &m_lvl);
  char *minus_rate_prefix = NULL;
  char minus_rate_prefix_storage[BFAM_BUFSIZ];
  if (m_lvl > -1 && m_lvl <= data->lvl)
  {
    minus_rate_prefix = minus_rate_prefix_storage;
    snprintf(minus_rate_prefix_storage, BFAM_BUFSIZ, "%s%d_",
             BFAM_LOCAL_ADAMS_PREFIX,
             data->ts->currentStageArray[m_lvl] % data->ts->nStages);
    BFAM_LDEBUG("Local Adams inter: level %" BFAM_LOCIDX_PRId
                " using minus rate prefix %s",
                m_lvl, minus_rate_prefix_storage);

    bfam_locidx_t p_lvl = -1;
    if (sub->glue_p)
      bfam_critbit0_allprefixed(&sub->glue_p->tags, BFAM_LOCAL_ADAMS_LVL_PREFIX,
                                get_tag_level_number, &p_lvl);

/*
 * If my level is less than my neighbor and my neighbor is not being updated
 * then do the interpolation
 */
#ifdef BFAM_LOCAL_ADAMS_ALWAYS_INTERP
    if (data->ts->numStepsArray[m_lvl] > 0)
#else
    if (m_lvl < p_lvl && p_lvl > data->lvl)
#endif
    {
      interp_fields(data, sub, m_lvl, p_lvl);
    }
  }

  if (rate_prefix || minus_rate_prefix)
    data->ts->inter_rhs(sub, rate_prefix, minus_rate_prefix, "", data->ts->t);

  return 1;
}

static inline int bfam_ts_local_adams_do_update(bfam_subdomain_t *sub,
                                                const bfam_long_real_t *A,
                                                const bfam_ts_local_adams_t *ts,
                                                const bfam_long_real_t dt,
                                                const int nStages,
                                                const bfam_locidx_t lvl)
{
  BFAM_LDEBUG("BFAM_TS_LOCAL_ADAMS_DO_UPDATE");

  /* Loop through the stages to scale rates and add in */
  /*
   * nStages is the computing number of stages whereas ts->nStages is the
   * storage number of stages
   */
  for (int k = 0; k < nStages; k++)
  {
    char prefix[BFAM_BUFSIZ];
    snprintf(prefix, BFAM_BUFSIZ, "%s%d_", BFAM_LOCAL_ADAMS_PREFIX,
             (ts->currentStageArray[lvl] + ts->nStages - k) % ts->nStages);
    BFAM_LDEBUG("Adams step: stage %d of %d using prefix %s", k, nStages,
                prefix);
    ts->add_rates(sub, "", "", prefix, dt * A[k]);
  }
  return 1;
}

static int bfam_ts_local_adams_update(const char *key, void *val, void *arg)
{
  bfam_ts_local_adams_allprefix_t *data =
      (bfam_ts_local_adams_allprefix_t *)arg;
  bfam_subdomain_t *sub = (bfam_subdomain_t *)val;
  bfam_ts_local_adams_t *ts = data->ts;

  bfam_locidx_t lvl = -1;
  bfam_critbit0_allprefixed(&sub->tags, BFAM_LOCAL_ADAMS_LVL_PREFIX,
                            get_tag_level_number, &lvl);

  if (lvl >= 0 && !(data->step % (1 << lvl)))
    switch (BFAM_MIN(ts->numStepsArray[lvl] + 1, ts->nStages))
    {
    case 1:
    {
      bfam_long_real_t A[1] = {BFAM_LONG_REAL(1.0)};
      bfam_ts_local_adams_do_update(sub, A, ts, (1 << lvl) * data->dt, 1, lvl);
    }
    break;
    case 2:
    {
      bfam_long_real_t A[2] = {
          BFAM_LONG_REAL(3.0) / BFAM_LONG_REAL(2.0),
          BFAM_LONG_REAL(-1.0) / BFAM_LONG_REAL(2.0),
      };
      bfam_ts_local_adams_do_update(sub, A, ts, (1 << lvl) * data->dt, 2, lvl);
    }
    break;
    case 3:
    {
      bfam_long_real_t A[3] = {
          BFAM_LONG_REAL(23.0) / BFAM_LONG_REAL(12.0),
          BFAM_LONG_REAL(-4.0) / BFAM_LONG_REAL(3.0),
          BFAM_LONG_REAL(5.0) / BFAM_LONG_REAL(12.0),
      };
      bfam_ts_local_adams_do_update(sub, A, ts, (1 << lvl) * data->dt, 3, lvl);
    }
    break;
    case 4:
    {
      bfam_long_real_t A[4] = {
          BFAM_LONG_REAL(55.0) / BFAM_LONG_REAL(24.0),
          BFAM_LONG_REAL(-59.0) / BFAM_LONG_REAL(24.0),
          BFAM_LONG_REAL(37.0) / BFAM_LONG_REAL(24.0),
          BFAM_LONG_REAL(-3.0) / BFAM_LONG_REAL(8.0),
      };
      bfam_ts_local_adams_do_update(sub, A, ts, (1 << lvl) * data->dt, 4, lvl);
    }
    break;
    default:
      BFAM_ABORT("Adams-Bashforth order %d not implemented",
                 BFAM_MIN(ts->numStepsArray[lvl] + 1, ts->nStages));
    }
  return 1;
}

static void bfam_ts_local_adams_step(bfam_ts_t *a_ts, bfam_long_real_t dt,
                                     void *user_data)
{
  bfam_ts_local_adams_t *ts = (bfam_ts_local_adams_t *)a_ts;
  bfam_locidx_t num_steps = 1 << (ts->numLevels - 1);
  BFAM_LDEBUG("Number of steps for the local time stepper %" BFAM_LOCIDX_PRId,
              num_steps);

  dt /= num_steps;
  bfam_ts_local_adams_allprefix_t data;
  data.ts = ts;
  data.dt = dt;
  data.lvl = -1;

  /* determine the level of comm to do: max of this update and last update */
  bfam_locidx_t last_lvl = 0;
  for (bfam_locidx_t step = 0; step < num_steps; step++)
  {
    data.step = step;

    BFAM_LDEBUG("local time step number %" BFAM_LOCIDX_PRId, step);
    data.lvl = 0;

    /* loop through the levels */
    for (bfam_locidx_t lvl = 0; lvl < ts->numLevels; lvl++)
    {
      bfam_locidx_t chk = 1 << lvl;
      if (!(step % chk))
      {
        BFAM_LDEBUG("level %" BFAM_LOCIDX_PRId " to be updated", lvl);
        data.lvl = BFAM_MAX(data.lvl, lvl);
      }
    }
    bfam_locidx_t comm_lvl = BFAM_MAX(data.lvl, last_lvl);
    BFAM_LDEBUG("step %02" BFAM_LOCIDX_PRId
                ": update level %02" BFAM_LOCIDX_PRId
                ": and communication level %02" BFAM_LOCIDX_PRId,
                step, data.lvl, comm_lvl);

    /* start the communication */
    bfam_subdomain_comm_args_t *sub_comm_data =
        (bfam_subdomain_comm_args_t *)ts->comm_array[comm_lvl]->user_args;
    BFAM_ASSERT(sub_comm_data);
    BFAM_ASSERT(sub_comm_data->user_data == NULL);
    BFAM_ASSERT(sub_comm_data->user_prefix_function == NULL);

    sub_comm_data->user_prefix_function = comm_send_prefix;
    sub_comm_data->user_data = &data;
    bfam_communicator_start(ts->comm_array[comm_lvl]);

    /* Do the intra work for the levels to be udpated */
    bfam_dictionary_allprefixed_ptr(&ts->elems, "",
                                    &bfam_ts_local_adams_intra_rhs, &data);

    /* finish the communication */
    sub_comm_data->user_prefix_function = comm_recv_prefix;
    bfam_communicator_finish(ts->comm_array[comm_lvl]);

    /* reset the pointer to NULL */
    sub_comm_data->user_data = NULL;
    sub_comm_data->user_prefix_function = NULL;

    /* Do the inter work for the levels to be udpated */
    bfam_dictionary_allprefixed_ptr(&ts->elems, "",
                                    &bfam_ts_local_adams_inter_rhs, &data);

    /* Do the inter work for the levels to be udpated */
    bfam_dictionary_allprefixed_ptr(&ts->elems, "", &bfam_ts_local_adams_update,
                                    &data);

    /* set last update to next update */
    last_lvl = data.lvl;

    /* update the stage counters */
    for (bfam_locidx_t k = 0; k <= data.lvl; k++)
    {
      ts->numStepsArray[k]++;
      ts->currentStageArray[k] = (ts->currentStageArray[k] + 1) % ts->nStages;
    }
    ts->t += dt;
  }
}

void bfam_ts_local_adams_init(
    bfam_ts_local_adams_t *ts, bfam_domain_t *dom,
    bfam_ts_local_adams_method_t method, bfam_locidx_t num_lvl,
    bfam_domain_match_t subdom_match, const char **subdom_tags,
    bfam_domain_match_t comm_match, const char **comm_tags, MPI_Comm mpicomm,
    int mpitag, bfam_subdomain_comm_args_t *comm_data,
    void (*aux_rates)(bfam_subdomain_t *thisSubdomain, const char *prefix),
    void (*glue_rates)(bfam_subdomain_t *thisSubdomain, const char *prefix),
    void (*scale_rates)(bfam_subdomain_t *thisSubdomain,
                        const char *rate_prefix, const bfam_long_real_t a),
    void (*intra_rhs)(bfam_subdomain_t *thisSubdomain, const char *rate_prefix,
                      const char *minus_rate_prefix, const char *field_prefix,
                      const bfam_long_real_t t),
    void (*inter_rhs)(bfam_subdomain_t *thisSubdomain, const char *rate_prefix,
                      const char *minus_rate_prefix, const char *field_prefix,
                      const bfam_long_real_t t),
    void (*add_rates)(bfam_subdomain_t *thisSubdomain,
                      const char *field_prefix_lhs,
                      const char *field_prefix_rhs, const char *rate_prefix,
                      const bfam_long_real_t a),
    void (*add_rates_glue_p)(bfam_subdomain_t *thisSubdomain,
                             const char *field_prefix_lhs,
                             const char *field_prefix_rhs,
                             const char *rate_prefix, const bfam_long_real_t a),
    const int RK_init)
{
  BFAM_LDEBUG("LOCAL ADAMS INIT");

  /*
   * set up some preliminaries
   */
  bfam_ts_init(&ts->base, dom);
  bfam_dictionary_init(&ts->elems);
  ts->t = BFAM_LONG_REAL(0.0);
  ts->base.step = &bfam_ts_local_adams_step;

  /*
   * store the function calls
   */
  ts->scale_rates = scale_rates;
  ts->intra_rhs = intra_rhs;
  ts->inter_rhs = inter_rhs;
  ts->add_rates = add_rates;
  ts->add_rates_glue_p = add_rates_glue_p;

  /*
   * fast log2 computation for ints using bitwise operations from
   * http://stackoverflow.com/questions/994593/how-to-do-an-integer-log2-in-c
   */
  ts->numLevels = num_lvl;

  ts->currentStageArray = bfam_malloc(ts->numLevels * sizeof(bfam_locidx_t));
  ts->numStepsArray = bfam_malloc(ts->numLevels * sizeof(bfam_locidx_t));
  for (bfam_locidx_t k = 0; k < ts->numLevels; k++)
  {
    ts->numStepsArray[k] = 0;
    ts->currentStageArray[k] = 0;
  }

  ts->lsrk = NULL;
  ts->comm_array = NULL;

  switch (method)
  {
  default:
    BFAM_WARNING("Invalid Adams scheme, using ADAMS_3");
  case BFAM_TS_LOCAL_ADAMS_3:
    ts->nStages = 3;
    break;
  case BFAM_TS_LOCAL_ADAMS_1:
    ts->nStages = 1;
    break;
  case BFAM_TS_LOCAL_ADAMS_2:
    ts->nStages = 2;
    break;
  case BFAM_TS_LOCAL_ADAMS_4:
    ts->nStages = 4;
    break;
  }

  /*
   * get the subdomains and create rates we will need
   */
  bfam_subdomain_t *subs[dom->numSubdomains + 1];
  bfam_locidx_t numSubs = 0;
  bfam_domain_get_subdomains(dom, subdom_match, subdom_tags, dom->numSubdomains,
                             subs, &numSubs);
  for (int s = 0; s < numSubs; s++)
  {
    int rval = bfam_dictionary_insert_ptr(&ts->elems, subs[s]->name, subs[s]);
    BFAM_ABORT_IF_NOT(rval != 1, "Issue adding subdomain %s", subs[s]->name);

    for (int n = 0; n < ts->nStages; n++)
    {
      char aux_rates_name[BFAM_BUFSIZ];
      snprintf(aux_rates_name, BFAM_BUFSIZ, "%s%d_", BFAM_LOCAL_ADAMS_PREFIX,
               n);
      aux_rates(subs[s], aux_rates_name);
      glue_rates(subs[s], aux_rates_name);
    }
  }

  numSubs = 0;
  bfam_domain_get_subdomains(dom, comm_match, comm_tags, dom->numSubdomains,
                             subs, &numSubs);

/* this tracks the max level that I know about */
#ifdef BFAM_DEBUG
  int local_max_levels = 0;
#endif

  /* find the plus and minus levels for the glue grids */
  for (int s = 0; s < numSubs; s++)
  {
    bfam_locidx_t m_lvl = -1;
    bfam_locidx_t p_lvl = -1;

    BFAM_ASSERT(subs[s]->glue_m);
    BFAM_ASSERT(subs[s]->glue_p);

    bfam_critbit0_allprefixed(&subs[s]->glue_p->tags,
                              BFAM_LOCAL_ADAMS_LVL_PREFIX, get_tag_level_number,
                              &p_lvl);
    bfam_critbit0_allprefixed(&subs[s]->glue_m->tags,
                              BFAM_LOCAL_ADAMS_LVL_PREFIX, get_tag_level_number,
                              &m_lvl);

    BFAM_ASSERT(m_lvl >= 0);
    BFAM_ASSERT(p_lvl >= 0);

#ifdef BFAM_DEBUG
    /* just a sanity check since these should match */
    BFAM_ASSERT(subs[s]->glue_m->sub_m);
    bfam_locidx_t s_lvl = -1;
    bfam_critbit0_allprefixed(&subs[s]->glue_m->sub_m->tags,
                              BFAM_LOCAL_ADAMS_LVL_PREFIX, get_tag_level_number,
                              &s_lvl);
    BFAM_ASSERT(s_lvl == m_lvl);
#endif

    /* communicate as infrequently as needed, so use the higher level */
    char comm_lvl_tag[BFAM_BUFSIZ];
    bfam_ts_local_adams_fill_comm_level_tag(comm_lvl_tag, BFAM_BUFSIZ,
                                            BFAM_MAX(m_lvl, p_lvl));
    bfam_subdomain_add_tag(subs[s], comm_lvl_tag);

#ifdef BFAM_DEBUG
    /* update the max level */
    local_max_levels = BFAM_MAX(local_max_levels, BFAM_MAX(m_lvl, p_lvl));
#endif
  }

#ifdef BFAM_DEBUG
  BFAM_ASSERT(local_max_levels <= num_lvl);
#endif

  /* loop through all possible commmunication tags */
  char *local_comm_tags[ts->numLevels + 1];
  char tag_stor[BFAM_BUFSIZ * ts->numLevels];
  ts->comm_array = bfam_malloc(ts->numLevels * sizeof(bfam_ts_local_adams_t *));

  /* since we will use these make sure they are NULL */
  BFAM_ASSERT(comm_data->user_data == NULL);
  BFAM_ASSERT(comm_data->user_prefix_function == NULL);
  BFAM_ASSERT(comm_data->user_comm_info == NULL);
  BFAM_ASSERT(comm_data->user_put_send_buffer == NULL);
  BFAM_ASSERT(comm_data->user_get_recv_buffer == NULL);

  for (int k = 0; k < ts->numLevels; k++)
  {
    local_comm_tags[k] = &tag_stor[BFAM_BUFSIZ * k];
    bfam_ts_local_adams_fill_comm_level_tag(&tag_stor[BFAM_BUFSIZ * k],
                                            BFAM_BUFSIZ, k);
    local_comm_tags[k + 1] = NULL;

    /*
     * Set up the communicator we will use
     */
    ts->comm_array[k] = bfam_communicator_new(dom, BFAM_DOMAIN_OR,
                                              (const char **)local_comm_tags,
                                              mpicomm, mpitag, comm_data);
  }
}

void bfam_ts_local_adams_free(bfam_ts_local_adams_t *ts)
{
  BFAM_LDEBUG("LOCAL ADAMS FREE");
  if (ts->lsrk != NULL)
  {
    bfam_ts_lsrk_free(ts->lsrk);
    bfam_free(ts->lsrk);
  }
  ts->lsrk = NULL;
  for (int k = 0; k < ts->numLevels; k++)
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
  bfam_free(ts->numStepsArray);
  bfam_free(ts->currentStageArray);
  ts->t = NAN;
  bfam_ts_free(&ts->base);
}
