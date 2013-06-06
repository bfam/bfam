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
}
