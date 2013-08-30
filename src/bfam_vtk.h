#ifndef BFAM_VTK_H
#define BFAM_VTK_H

#include <bfam_base.h>
#include <bfam_domain.h>

/** Utility function for writing an empty vtu file
 *
 * \param [out] file            stream to write the vector to.
 * \param [in]  name            name of the vector.
 * \param [in]  writeBinary     boolean indicating if the data should be written
 *                              in binary.
 */
void
bfam_vtk_write_vtu_empty(FILE *file,int writeBinary);

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
 * \param [in] binary     boolean to indicate if the data should be
 *                        written in binary
 * \param [in] compress   boolean to indicate if the data should be compressed
 *                        when writing a binary data
 *
 *
 * \note That a domain will only output scalars and vectors that it contains.
 *
 */

void
bfam_vtk_write_struc_file(bfam_domain_t *domain, bfam_domain_match_t match,
    const char **tags, const char *prefix, const char **scalars,
    const char **vectors, const char **components, int binary,
    int compress);

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
 * \param [in] binary     boolean to indicate if the data should be
 *                        written in binary
 * \param [in] compress   boolean to indicate if the data should be compressed
 *                        when writing a binary data
 *
 *
 * \note That a domain will only output scalars and vectors that it contains.
 *
 */

void
bfam_vtk_write_file(bfam_domain_t *domain, bfam_domain_match_t match, const
    char **tags, const char *prefix, const char **scalars,
    const char **vectors, const char **components, int binary,
    int compress);


/** Utility function to write binary data in VTK format.
 *
 * Currently this is just a wrapper to call a similar function in libsc.
 *
 * \param [in]  compressed boolean specifying if the binary data should be
 *                         compressed.
 * \param [out] file       stream to write the data to.
 * \param [in]  data       data to write out.
 * \param [in]  size       size of the data in bytes.
 *
 * \returns 0 on success and -1 on file error.
 */
int
bfam_vtk_write_binary_data(int compressed, FILE *file, char *data, size_t size);


/** Utility function for writing a scalar data array.
 *
 * \param [out] file            stream to write the scalar to.
 * \param [in]  name            name of the scalar.
 * \param [in]  writeBinary     boolean indicating if the data should be written
 *                              in binary.
 * \param [in]  writeCompressed boolean indicating if the data should be
 *                              compressed.
 * \param [in]  Ntotal          length of the scalar.
 * \param [in]  s               scalar data.
 */
void
bfam_vtk_write_real_scalar_data_array(FILE* file, const char *name,
    int writeBinary, int writeCompressed, bfam_locidx_t Ntotal,
    const bfam_real_t *s);

/** Utility function for writing a vector data array.
 *
 * \param [out] file            stream to write the vector to.
 * \param [in]  name            name of the vector.
 * \param [in]  writeBinary     boolean indicating if the data should be written
 *                              in binary.
 * \param [in]  writeCompressed boolean indicating if the data should be
 *                              compressed.
 * \param [in]  Ntotal          length of the vector components.
 * \param [in]  v1              1st component of the vector.
 * \param [in]  v2              2nd component of the vector.
 * \param [in]  v3              3rd component of the vector.
 */
void
bfam_vtk_write_real_vector_data_array(FILE* file, const char *name,
    int writeBinary, int writeCompressed, bfam_locidx_t Ntotal,
    const bfam_real_t *v1, const bfam_real_t *v2, const bfam_real_t *v3);

#endif
