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
      ts->A = bfam_malloc_aligned(ts->nStages*sizeof(bfam_long_real_t));
      ts->B = bfam_malloc_aligned(ts->nStages*sizeof(bfam_long_real_t));
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
  bfam_free_aligned(ts->A);
  bfam_free_aligned(ts->B);
  bfam_free_aligned(ts->C);
  ts->nStages = 0;
  ts->t  = NAN;
  ts->dt = NAN;
  bfam_ts_free(&ts->p_ts);
}
