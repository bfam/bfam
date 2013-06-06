#include <bfam_timestep.h>

void
bfam_ts_init(bfam_ts_t* ts, bfam_domain_t* dom)
{
  ts->domain = dom;
}
