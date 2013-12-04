#ifndef BFAM_SUBDOMAIN_H
#define BFAM_SUBDOMAIN_H

#include <bfam_base.h>
#include <bfam_critbit.h>
#include <bfam_dictionary.h>

typedef struct bfam_subdomain_comm_args
{
  const char ** scalars_m;           /* \c NULL terminated array of scalars to
                                      * send */

  const char ** vectors_m;           /* \c NULL terminated array of vectors to
                                      * send
                                      * \note will suffix n for normal component
                                      *       and p[1-3] for perpendicular
                                      *       components */

  const char ** vector_components_m; /* \c NULL terminated array of vectors
                                      * components
                                      * \note must be three components per
                                      *       vector to send and a \c NULL entry
                                      *       will lead to 0 being used for that
                                      *       component */

  const char ** tensors_m;           /* \c NULL terminated array of tensor to send
                                      * \note will suffix n for normal component
                                      *       and p[1-3] for perpendicular
                                      *       components */

  const char ** tensor_components_m; /* \c NULL terminated array of symetric
                                      * tensor components in order
                                      * {11,22,33,12,13,23}
                                      * \note must be six components per tensor
                                      *       to send and \c NULL entry will
                                      *       lead to 0 being used for that
                                      *       component */

  const char ** face_scalars_m;      /* \c NULL terminated array of face scalars
                                      *  to send */
  const char ** scalars_p;
  const char ** vectors_p;
  const char ** vector_components_p;
  const char ** tensors_p;
  const char ** tensor_components_p;
  const char ** face_scalars_p;
} bfam_subdomain_comm_args_t;

typedef struct bfam_subdomain_face_map_entry
{
  bfam_locidx_t np; /* Neighbor's processor number */
  bfam_locidx_t ns; /* Neighbor's subdomain id */
  bfam_locidx_t nk; /* Neighbor's element number */
  int8_t        nf; /* Neighbor's face number */
  int8_t        nh; /* Neighbor's hanging number */

  bfam_locidx_t  s; /* Local subdomain id */
  bfam_locidx_t  k; /* Local element number */
  int8_t         f; /* Local face number */
  int8_t         h; /* Local hanging number */
  int8_t         o; /* Local orientation */

  bfam_locidx_t id; /* Interface id */

  bfam_locidx_t gi; /* Index variable */
  bfam_locidx_t  i; /* Index variable */
} bfam_subdomain_face_map_entry_t;

/*
 * Compare function which sorts a bfam_subdomain_face_map_entry_t array
 * in sending order.
 */
int
bfam_subdomain_face_send_cmp(const void *a, const void *b);

/*
 * Compare function which sorts a bfam_subdomain_face_map_entry_t array
 * in receiving order.
 */
int
bfam_subdomain_face_recv_cmp(const void *a, const void *b);

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
    const char *name,
    bfam_real_t time, bfam_real_t *restrict x, bfam_real_t *restrict y,
    bfam_real_t *restrict z, struct bfam_subdomain *s, void *arg,
    bfam_real_t *restrict field);


/**
 * base structure for to store glue data for a subdomain, i.e., plus and minus
 * sides stuff. It should also be included as first member with the name base
 * \code{.c}
 * typedef struct new_subdomain_glue_data_type
 * {
 *   bfam_subdomain_t base;
 *   ...
 * }
 */
typedef struct bfam_subdomain_glue_data
{
  bfam_locidx_t     rank; /* Rank of the subdomain on this side */
  bfam_locidx_t     id_s; /* Sort Id of the subdomain on this side */
  bfam_locidx_t       id; /* Id of the subdomain on this side */

  bfam_dictionary_t fields; /**< a dictionary storing glue fields */

  /* The following pointers should only be \ne NULL on the minus side */
  struct bfam_subdomain *sub_m;  /* Local neighboring subdomain */
} bfam_subdomain_glue_data_t;

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
  bfam_locidx_t     id;
  char*           name;     /**< Name of the subdomain */
  bfam_critbit0_tree_t tags; /**< critbit for tags for the subdomain */
  bfam_dictionary_t fields; /**< a dictionary storing pointers to fields */

  bfam_dictionary_t fields_face; /**< a dictionary storing face fields */

  /* glue quantities */
  bfam_subdomain_glue_data_t *glue_m;
  bfam_subdomain_glue_data_t *glue_p;

  /* This storage is depreciated and will be removed in the future */
  bfam_dictionary_t fields_m; /**< a dictionary storing minus fields */
  bfam_dictionary_t fields_p; /**< a dictionary storing plus fields */


  /* Function pointers that domain will need to call */
  void (*free)                (struct bfam_subdomain *thisSubdomain);

  /**< Write a vtk vtu file */
  int (*vtk_write_vtu_piece) (struct bfam_subdomain *thisSubdomain,
                               FILE *file,
                               bfam_real_t time,
                               const char **scalars,
                               const char **vectors,
                               const char **components,
                               int writeBinary,
                               int writeCompressed,
                               int rank,
                               bfam_locidx_t id,
                               int Np_write);

  /**< Write a vtk vts file */
  void (*vtk_write_vts_piece) (struct bfam_subdomain *thisSubdomain,
                               FILE *file,
                               const char **scalars,
                               const char **vectors,
                               const char **components,
                               int writeBinary,
                               int writeCompressed,
                               int rank);

  /**< vts filename */
  void (*vtk_write_vts_filename) (struct bfam_subdomain *thisSubdomain,
      void *buffer, size_t recv_sz);

  /**< write the filename suffix */
  void (*vtk_write_suffix) (struct bfam_subdomain *thisSubdomain,char * suffix,
      int len);

  /**< write the pvts filename */
  void (*vtk_write_pvts_pieces) (struct bfam_subdomain *thisSubdomain,
      const char* filename, FILE *file);

  /**< Add a field to the subdomain */
  int (*field_add) (struct bfam_subdomain *thisSubdomain, const char* name);

  /**< Add a field to the plus side of the subdomain */
  int (*field_plus_add) (struct bfam_subdomain *thisSubdomain,
                         const char* name);
  /**< Add a field to the minus side of the  subdomain */
  int (*field_minus_add) (struct bfam_subdomain *thisSubdomain,
                          const char* name);

  /**< Add a field to the faces of the subdomain */
  int (*field_face_add) (struct bfam_subdomain *thisSubdomain,
                         const char* name);


  /**< Initialize a field in the subdomain */
  void (*field_init) (struct bfam_subdomain *thisSubdomain, const char* name,
      bfam_real_t time, bfam_subdomain_init_field_t init_field, void *arg);

  /**< Glue grid communication info */
  void (*glue_comm_info) (struct bfam_subdomain *thisSubdomain, int *rank,
      bfam_locidx_t *s, int num_sort, size_t *send_sz, size_t *recv_sz,
      void *args);

  /**< Put data into the send buffer */
  void (*glue_put_send_buffer) (struct bfam_subdomain *thisSubdomain,
      void *buffer, size_t send_sz, void *args);

  /**< Get data from the recv buffer */
  void (*glue_get_recv_buffer) (struct bfam_subdomain *thisSubdomain,
      void *buffer, size_t recv_sz, void *args);
} bfam_subdomain_t;


/** initializes a subdomain
 *
 * There is no new function for subdomains since these are really just a base
 * class and a concrete grid and physics type should be defined
 *
 * \param [in,out] thisSubdomain pointer to the subdomain
 * \param [in]     id   Unique id number for this subdomain
 * \param [in]     name Name of this subdomain
 */
void
bfam_subdomain_init(bfam_subdomain_t *subdomain, bfam_locidx_t id,
    const char* name);

/** initializes a subdomain glue data
 *
 * There is no new function for subdomains glue data since these are really just
 * a base class and a concrete grid and physics type should be defined
 *
 * \param [in,out] thisGlue pointer to the subdomain glue data
 * \param [in]     id of this side
 * \param [in]     sort id of this side
 * \param [in]     mpirank for this glue subdomain data
 * \param [in]     pointer to the minus side subdomain (can be \c NULL);
 */
void
bfam_subdomain_glue_init(bfam_subdomain_glue_data_t *glue,
    const bfam_locidx_t rank, const bfam_locidx_t id, const bfam_locidx_t id_s,
    bfam_subdomain_t *sub_m);

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

/** Add a field to the plus side of the subdomain.
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
bfam_subdomain_field_plus_add(bfam_subdomain_t *thisSubdomain,
                              const char* name);

/** Add a field to the minus side of the subdomain.
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
bfam_subdomain_field_minus_add(bfam_subdomain_t *thisSubdomain,
                               const char* name);

/** Add a field to the faces of the subdomain.
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
bfam_subdomain_field_face_add(bfam_subdomain_t *thisSubdomain,
                              const char* name);

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
