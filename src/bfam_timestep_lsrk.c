#include <bfam_timestep_lsrk.h>

bfam_ts_lsrk_t*
bfam_ts_lsrk_new(bfam_domain_t* dom, bfam_ts_lsrk_method_t method)
{
  bfam_ts_lsrk_t* newTS = bfam_malloc(sizeof(bfam_ts_lsrk_t));
  bfam_ts_lsrk_init(newTS, dom, method);
  return newTS;
}


void
bfam_ts_lsrk_init(bfam_ts_lsrk_t* ts, bfam_domain_t* dom,
    bfam_ts_lsrk_method_t method)
{
  bfam_ts_init(&ts->p_ts, dom);
  ts->t  = 0.0;
  ts->dt = NAN;
  switch(method)
  {
    case BFAM_TS_LSRK_KC54:
      break;
    case BFAM_TS_LSRK_FE:
      break;
    case BFAM_TS_LSRK_HEUN:
      break;
    case BFAM_TS_LSRK_W33:
      break;
  }
}

void
bfam_ts_lsrk_free(bfam_ts_lsrk_t* ts)
{
  bfam_free(ts->A);
  bfam_free(ts->B);
  bfam_free(ts->C);
  ts->nStages = 0;
  ts->t  = NAN;
  ts->dt = NAN;
  bfam_ts_free(&ts->p_ts);
}
