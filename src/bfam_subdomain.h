#ifndef BFAM_SUBDOMAIN_H
#define BFAM_SUBDOMAIN_H

/**
 * base structure for all subdomains types. Any new subdomain should have this
 * as its first member with the name base, i.e.,
 * \code{.c}
 * typedef struct new_subdomain_type
 * {
 *   bfam_subdomain_t base;
 *   ...
 * }
 */
typedef struct bfam_subdomain
{
} bfam_subdomain_t;

/**
 * structure containing a domain (which is a collection of subdomains!)
 */
typedef struct domain 
{
  bfam_subdomain_t * domain; /**< array of pointers to subdomains */ int
  numSubdomains;             /**< number of subdomains that are currently in the domain */
  int sizeSubdomains;        /**< total number of subdomains the domain can
                               hold, i.e.  size of the array*/
} domain_t;

#endif
