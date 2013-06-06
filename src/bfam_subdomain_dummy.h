#ifndef BFAM_SUBDOMAIN_DUMMY_H
#define BFAM_SUBDOMAIN_DUMMY_H

#include <bfam_base.h>
#include <bfam_subdomain.h>


typedef struct bfam_subdomain_dummy
{
  bfam_subdomain_t base;
  int              N;    /**< order of scheme to test */
} bfam_subdomain_dummy_t;

/** create a dummy subdomain.
 *
 * \param [in] name name of this subdomain
 * \param [in] N    polynomial order of elements in each dimension
 *
 * \return Initialized dummy subdomain
 *
 */
bfam_subdomain_dummy_t*
bfam_subdomain_dummy_new(const char             *name,
                            const int               N);

/** initializes a dg quad subdomain
 *
 * \param [in,out] subdomain pointer to the subdomain to initialize
 * \param [in]     name      name of this subdomain
 * \param [in]     N         polynomial order of elements in each dimension
 *
 */
void
bfam_subdomain_dummy_init(bfam_subdomain_dummy_t *subdomain,
                             const char                *name,
                             const int                  N);

/** free up the memory allocated by the subdomain
 *
 * \param [in,out] subdomain subdomain to clean up
 *
 */
void
bfam_subdomain_dummy_free(bfam_subdomain_t *subdomain);

/** compute exact solution for dummy domain at time t
 *
 * \param [in] subdom   pointer to domain to calculate solution for
 * \param [in] t        time to compute exact solution at
 *
 */
bfam_long_real_t
bfam_subdomain_dummy_exact(bfam_subdomain_dummy_t *subdom, bfam_long_real_t t);

#endif
