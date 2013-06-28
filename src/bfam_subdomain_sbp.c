#include <bfam_subdomain_sbp.h>
#include <bfam_log.h>
#include <bfam_vtk.h>
#include <bfam_util.h>
#include <inttypes.h>

int bfam_subdomain_sbp_vtk_write_suffix(bfam_subdomain_t *thisSubdomain)
{
  return ((bfam_subdomain_sbp_t*)thisSubdomain)->loc_id;
}

bfam_subdomain_sbp_t*
bfam_subdomain_sbp_new(const bfam_locidx_t     id,
                            const bfam_locidx_t loc_id,
                            const bfam_locidx_t num_id,
                            const char             *name,
                            const int               dim,
                            const bfam_gloidx_t    *N,
                            const bfam_locidx_t    *Nl,
                            const bfam_locidx_t    *Nb,
                            const bfam_gloidx_t    *gx,
                            const bfam_long_real_t *c_x,
                            const bfam_long_real_t *c_y,
                            const bfam_long_real_t *c_z)
{
  bfam_subdomain_sbp_t *newSub = bfam_malloc(sizeof(bfam_subdomain_sbp_t));
  bfam_subdomain_sbp_init(newSub,id,loc_id,num_id,name,dim,N,Nl,Nb,gx,
                          c_x,c_y,c_z);
  return newSub;
}

void bfam_subdomain_sbp_vtk_vts_piece (struct bfam_subdomain *subdomain,
                                       FILE *file,
                                       const char **scalars,
                                       const char **vectors,
                                       const char **components,
                                       int writeBinary,
                                       int writeCompressed)
{
  bfam_subdomain_sbp_t* s = (bfam_subdomain_sbp_t*) subdomain;

  bfam_locidx_t Nl[3] = {0,0,0};
  bfam_locidx_t Nb[6] = {0,0,0,0,0,0};
  bfam_locidx_t Nt[3] = {0,0,0};
  bfam_locidx_t gx[3] = {0,0,0};
  for(int d = 0;d < s->dim; d++)
  {
    Nl[d] = s->Nl[d];
    Nb[2*d  ] = s->Nb[2*d  ];
    Nb[2*d+1] = s->Nb[2*d+1];
    Nt[d] = Nb[2*d] + Nl[d] + Nb[2*d+1];
    gx[d] = s->gx[d];
  }

  /* start StructuredGrid */
  fprintf(file, "  <StructuredGrid WholeExtent=\"");
  fprintf(file, "%" PRId64 " %" PRId64,
      (int64_t)gx[0],(int64_t)(gx[0]+Nl[0]));
  fprintf(file, " %" PRId64 " %" PRId64,
      (int64_t)gx[1],(int64_t)(gx[1]+Nl[1]));
  fprintf(file, " %" PRId64 " %" PRId64,
      (int64_t)gx[2],(int64_t)(gx[2]+Nl[2]));
  fprintf(file, "\">\n");

  /* write the piece data */
  fprintf(file, "    <Piece Extent=\"");
  fprintf(file, "%" PRId64 " %" PRId64,
      (int64_t)gx[0],(int64_t)(gx[0]+Nl[0]));
  fprintf(file, " %" PRId64 " %" PRId64,
      (int64_t)gx[1],(int64_t)(gx[1]+Nl[1]));
  fprintf(file, " %" PRId64 " %" PRId64,
      (int64_t)gx[2],(int64_t)(gx[2]+Nl[2]));
  fprintf(file, "\">\n");

  /*
   * Points
   */

  bfam_locidx_t sz = 1;
  for(int i = 0; i < s->dim; i++) sz *= (s->Nl[i]+1);
  bfam_real_t *storage0 = bfam_malloc_aligned(sz*sizeof(bfam_real_t));
  bfam_real_t *storage1 = bfam_malloc_aligned(sz*sizeof(bfam_real_t));
  bfam_real_t *storage2 = bfam_malloc_aligned(sz*sizeof(bfam_real_t));

  fprintf(file, "      <Points>\n");

  bfam_real_t *x = bfam_dictionary_get_value_ptr(&subdomain->fields,"_grid_x");
  bfam_real_t *y = bfam_dictionary_get_value_ptr(&subdomain->fields,"_grid_y");
  bfam_real_t *z = bfam_dictionary_get_value_ptr(&subdomain->fields,"_grid_z");

  if(x!=NULL)
  {
    int i = 0;
    for(int k = 0; k <= Nl[2]; k++)
      for(int j = 0; j <= Nl[1]; j++)
        memcpy(&storage0[i + (Nl[0]+1)*(j + k*(Nl[1]+1))],
            &x[i+Nb[0] + (Nt[0]+1)*(j+Nb[2] + (k+Nb[4])*(Nt[1]+1))],
            (Nl[0]+1)*sizeof(bfam_real_t));
  }
  else
  {
    memset(storage0,0,sz);
  }

  if(y!=NULL)
  {
    int i = 0;
    for(int k = 0; k <= Nl[2]; k++)
      for(int j = 0; j <= Nl[1]; j++)
        memcpy(&storage1[i + (Nl[0]+1)*(j + k*(Nl[1]+1))],
            &y[i+Nb[0] + (Nt[0]+1)*(j+Nb[2] + (k+Nb[4])*(Nt[1]+1))],
            (Nl[0]+1)*sizeof(bfam_real_t));
  }
  else
  {
    memset(storage1,0,sz);
  }

  if(z!=NULL)
  {
    int i = 0;
    for(int k = 0; k <= Nl[2]; k++)
      for(int j = 0; j <= Nl[1]; j++)
        memcpy(&storage2[i + (Nl[0]+1)*(j + k*(Nl[1]+1))],
            &z[i+Nb[0] + (Nt[0]+1)*(j+Nb[2] + (k+Nb[4])*(Nt[1]+1))],
            (Nl[0]+1)*sizeof(bfam_real_t));
  }
  else
  {
    memset(storage2,0,sz);
  }

  bfam_vtk_write_real_vector_data_array(file, "Position", writeBinary,
      writeCompressed, sz, storage0, storage1, storage2);

  fprintf(file, "      </Points>\n");

  /*
   * Cells
   */
  fprintf(file, "      <Cells>\n");
  fprintf(file, "      </Cells>\n");

  fprintf(file, "    </Piece>\n");

  /* close StructuredGrid */
  fprintf(file, "  </StructuredGrid>\n");

  bfam_free_aligned(storage0);
  bfam_free_aligned(storage1);
  bfam_free_aligned(storage2);
}

static int
bfam_subdomain_sbp_field_add(bfam_subdomain_t *subdomain, const char *name)
{
  bfam_subdomain_sbp_t *s = (bfam_subdomain_sbp_t*) subdomain;

  if(bfam_dictionary_get_value_ptr(&s->base.fields,name))
    return 1;

  size_t fieldSize = 0;
  for(int d = 0; d < s->dim; d++)
    fieldSize *= (s->Nl[d]+1+s->Nb[2*d]+s->Nb[2*d+1]);

  bfam_real_t *field = bfam_malloc_aligned(fieldSize);

  int rval = bfam_dictionary_insert_ptr(&s->base.fields, name, field);

  BFAM_ASSERT(rval != 1);

  if(rval == 0)
    bfam_free_aligned(field);

  return rval;
}

static void
bfam_subdomain_sbp_field_init(bfam_subdomain_t *subdomain,
    const char *name, bfam_real_t time, bfam_subdomain_init_field_t init_field,
    void *arg)
{
  bfam_subdomain_sbp_t *s = (bfam_subdomain_sbp_t*) subdomain;

  bfam_real_t *field = bfam_dictionary_get_value_ptr(&s->base.fields,name);

  BFAM_ABORT_IF(field==NULL, "Init: Field %s not found in subdomain %s",
      name, subdomain->name);

  // size_t fieldLength = s->Np*s->K;

  // init_field(fieldLength, time, s->x, s->y, s->z, subdomain, arg, field);
}

void
bfam_subdomain_sbp_init(bfam_subdomain_sbp_t *subdomain,
                            const bfam_locidx_t     id,
                            const bfam_locidx_t loc_id,
                            const bfam_locidx_t num_id,
                            const char             *name,
                            const int               dim,
                            const bfam_gloidx_t    *N,
                            const bfam_locidx_t    *Nl,
                            const bfam_locidx_t    *Nb,
                            const bfam_gloidx_t    *gx,
                            const bfam_long_real_t *c_x,
                            const bfam_long_real_t *c_y,
                            const bfam_long_real_t *c_z)
{
  bfam_subdomain_init(&subdomain->base, id, name);
  bfam_subdomain_add_tag(&subdomain->base, "_subdomain_sbp");
  subdomain->base.free = bfam_subdomain_sbp_free;
  subdomain->base.field_add = bfam_subdomain_sbp_field_add;
  subdomain->base.field_init = bfam_subdomain_sbp_field_init;
  subdomain->base.vtk_write_vts_piece = bfam_subdomain_sbp_vtk_vts_piece;
  subdomain->base.vtk_write_suffix = bfam_subdomain_sbp_vtk_write_suffix;

  subdomain->loc_id = loc_id;
  subdomain->num_id = num_id;
  subdomain->dim = dim;

  subdomain->N  = bfam_malloc(  dim*sizeof(bfam_gloidx_t));
  subdomain->Nl = bfam_malloc(  dim*sizeof(bfam_locidx_t));
  subdomain->Nb = bfam_malloc(2*dim*sizeof(bfam_locidx_t));
  subdomain->gx = bfam_malloc(  dim*sizeof(bfam_gloidx_t));
  bfam_gloidx_t gxb[dim];

  bfam_locidx_t Nltmp[dim];

  bfam_locidx_t sz = 1;
  for(int d = 0; d < dim; d++)
  {
    subdomain->N [  d  ] = N [  d  ];
    subdomain->Nl[  d  ] = Nl[  d  ];
    subdomain->Nb[2*d  ] = Nb[2*d  ];
    subdomain->Nb[2*d+1] = Nb[2*d+1];
    subdomain->gx[  d  ] = gx[  d  ];
    gxb[d] = gx[d]-Nb[2*d];

    Nltmp[d] = (Nb[2*d] + Nl[d] + Nb[2*d+1]);
    sz *= (Nltmp[d]+1);

    BFAM_ABORT_IF(N [d] < 0 || Nl[d] < 0 || Nb[2*d] < 0 || Nb[2*d+1] < 0 ||
                  gx[d] < 0 || gx[d] > N[d] || gx[d]+Nl[d] > N[d],
                  "%s: problem with a dimension %d", d, subdomain->base.name);
  }

  bfam_real_t *x = NULL;
  if(c_x != NULL)
  {
    x = bfam_malloc_aligned(sizeof(bfam_real_t)*sz);
    int val = bfam_dictionary_insert_ptr(&subdomain->base.fields, "_grid_x", x);
    BFAM_ABORT_IF(val != 2, "problem adding x to fields");
  }

  bfam_real_t *y = NULL;
  if(c_y != NULL)
  {
    BFAM_ASSERT(c_x != NULL);
    y = bfam_malloc_aligned(sizeof(bfam_real_t)*sz);
    int val = bfam_dictionary_insert_ptr(&subdomain->base.fields, "_grid_y", y);
    BFAM_ABORT_IF(val != 2, "problem adding y to fields");
  }

  bfam_real_t *z = NULL;
  if(c_z != NULL)
  {
    BFAM_ASSERT(c_y != NULL);
    z = bfam_malloc_aligned(sizeof(bfam_real_t)*sz);
    int val = bfam_dictionary_insert_ptr(&subdomain->base.fields, "_grid_z", z);
    BFAM_ABORT_IF(val != 2, "problem adding z to fields");
  }

  bfam_util_linear_blend(x,y,z,dim,N,Nltmp,gxb,c_x,c_y,c_z);
}

static int
bfam_subdomain_sbp_free_fields(const char * key, void *val,
    void *arg)
{
  bfam_free_aligned(val);

  return 1;
}

void
bfam_subdomain_sbp_free(bfam_subdomain_t *subdomain)
{
  bfam_subdomain_sbp_t *sub = (bfam_subdomain_sbp_t*) subdomain;

  bfam_dictionary_allprefixed_ptr(&sub->base.fields,"",
      &bfam_subdomain_sbp_free_fields,NULL);

  bfam_subdomain_free(subdomain);

  bfam_free(sub->N );
  bfam_free(sub->Nl);
  bfam_free(sub->Nb);
  bfam_free(sub->gx);
}
