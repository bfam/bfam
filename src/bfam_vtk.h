#ifndef BFAM_VTK_H
#define BFAM_VTK_H

#include <bfam_domain.h>

/** Write out vtk files for each domain.
 *
 * This loops through the subdomains and writes out a vtk file for each one.
 *
 * \param [in] domain     domain to output to vtk files
 * \param [in] match      type of match, \c BFAM_DOMAIN_OR will
 *                        match subdomains with any of the tags
 *                        and \c BFAM_DOMAIN_AND will match subdomains
 *                        with all of the tags.
 * \param [in] tags       \c NULL terminated array of the tags to match
 * \param [in] prefix     prefix for the vtk files
 * \param [in] scalars    \c NULL terminated array of scalars to match
 * \param [in] vectors    \c NULL terminated array of vector that will
 *                        be outputted
 * \param [in] components \c NULL terminated array of vector components
 *
 * \note That a domain will only output scalars and vectors that it contains.
 *
 */

void
bfam_vtk_write_file(bfam_domain_t *domain, bfam_domain_match_t match, const
    char **tags, const char *prefix, const char **scalars,
    const char **vectors, const char **components);

#endif
