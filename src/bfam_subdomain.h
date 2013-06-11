#ifndef BFAM_SUBDOMAIN_H
#define BFAM_SUBDOMAIN_H

#include <bfam_base.h>
#include <bfam_mpicomm.h>
#include <bfam_critbit.h>
#include <bfam_dictionary.h>

struct bfam_subdomain;

/*
 * This is the field initialization function which will be called for each
 * subdomain.  This function fills the \a npoints values of \a field.
 *
 * \param [in]  npoints   this is the length of x, y, z, and field
 * \param [in]  x         x-coordinates of the points
 * \param [in]  y         y-coordinates of the points
 * \param [in]  z         z-coordinates of the points
 * \param [in]  s         pointer to the subdomain that the field is in
 * \param [in]  arg       user pointer
 * \param [out] field     the field values that need to be set
 *
 */
typedef void (*bfam_subdomain_init_field_t) (bfam_locidx_t npoints,
    bfam_real_t time, bfam_real_t *restrict x, bfam_real_t *restrict y,
    bfam_real_t *restrict z, struct bfam_subdomain *s, void *arg,
    bfam_real_t *restrict field);


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
  char*           name;     /**< Name of the subdomain */
  bfam_mpicomm_t* comm;     /**< communicator for this subdomain */
  bfam_critbit0_tree_t tags; /**< critbit for tags for the subdomain */
  bfam_dictionary_t fields; /**< a dictionary storing pointers to files */

  /* Function pointers that domain will need to call */
  void (*free)                (struct bfam_subdomain *thisSubdomain);

  /**< Write a vtk file */
  void (*vtk_write_vtu_piece) (struct bfam_subdomain *thisSubdomain,
                               FILE *file,
                               const char **scalars,
                               const char **vectors,
                               const char **components,
                               int writeBinary,
                               int writeCompressed,
                               int rank,
                               bfam_locidx_t id);

  /**< Add a field to the subdomain */
  int (*field_add) (struct bfam_subdomain *thisSubdomain, const char* name);

  /**< Initialize a field in the subdomain */
  void (*field_init) (struct bfam_subdomain *thisSubdomain, const char* name,
      bfam_real_t time, bfam_subdomain_init_field_t init_field, void *arg);
} bfam_subdomain_t;


/** initializes a subdomain
 *
 * There no new function for subdomains since these are really just a base class
 * and a concrete grid and physics type should be defined
 *
 * \param [in,out] thisSubdomain pointer to the subdomain
 * \param [in]     name Name of this subdomain
 */
void
bfam_subdomain_init(bfam_subdomain_t *subdomain,const char* name);

/** free up the memory allocated by the subdomain
 * 
 * \param [in,out] thisSubdomain subdomain to clean up
 */
void
bfam_subdomain_free(bfam_subdomain_t *thisSubdomain);

/** Add a tag to the subdomain
 *
 * \param [in,out] thisSubdomain subdomain to andd the tag to
 * \param [in]     tag           tag of the domain (\0 terminated string)
 *
 */
void
bfam_subdomain_add_tag(bfam_subdomain_t *thisSubdomain, const char* tag);

/** Remove a tag from the subdomain
 *
 * \param [in,out] thisSubdomain subdomain to remove the tag from
 * \param [in]     tag           tag of the domain (\0 terminated string)
 *
 * \returns It returns 1 if the tag was removed, 0 otherwise.
 */
int
bfam_subdomain_delete_tag(bfam_subdomain_t *thisSubdomain, const char* tag);

/** Check to see if a subdomain has a tag
 *
 * \param [in,out] thisSubdomain subdomain to search for the tag
 * \param [in]     tag           tag of the domain (\0 terminated string)
 *
 * \return nonzero iff \a thisSubdomain has the tag \a tag
 */
int
bfam_subdomain_has_tag(bfam_subdomain_t *thisSubdomain, const char* tag);

/** Add a field to the subdomain
 *
 * \param [in,out] thisSubdomain subdomain to search for the tag
 * \param [in]     name          name of the field to add to the subdomain
 *
 * \returns:
 *   $\cases{ 0 &if {\rm out of memory} \cr
 *            1 &if {\it name} {\rm was already a field} \cr
 *            2 &if {\it name} {\rm was added successfully}}$.
 */
int
bfam_subdomain_field_add(bfam_subdomain_t *thisSubdomain, const char* name);

/** Initialize a field in the subdomain
 *
 * \param [in,out] thisSubdomain subdomain to search for the tag
 * \param [in]     name          name of the field to initialize
 * \param [in]     time          time to pass to initilization function
 * \param [in]     init_field     initilization function
 * \param [in]     arg           user pointer
 */
void
bfam_subdomain_field_init(bfam_subdomain_t *thisSubdomain, const char* name,
      bfam_real_t time, bfam_subdomain_init_field_t init_field, void *arg);

#endif
