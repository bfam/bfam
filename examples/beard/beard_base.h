#ifndef BEARD_BASE_H
#define BEARD_BASE_H

#include <bfam.h>

typedef bfam_real_t (*beard_user_bc_t)(bfam_locidx_t uid, bfam_long_real_t t,
                                       const bfam_real_t *x,
                                       const bfam_real_t *n, bfam_real_t *TpP,
                                       bfam_real_t *TnP, bfam_real_t *vpP,
                                       bfam_real_t *vnP, void *user_data);

#endif
