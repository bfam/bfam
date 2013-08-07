#include <bfam_subdomain_sbp.h>
#include <bfam_log.h>
#include <bfam_vtk.h>
#include <bfam_util.h>
#include <inttypes.h>

void
bfam_subdomain_sbp_vtk_write_suffix(bfam_subdomain_t *thisSubdomain,
    char* suffix, int len)
{
  bfam_subdomain_sbp_t * s = (bfam_subdomain_sbp_t*) thisSubdomain;
  snprintf(suffix,len,"%05d",s->sub_ix[0]);
  for(int d = 1; d < s->dim;d++)
  {
    char tmpstr[BFAM_BUFSIZ];
    snprintf(tmpstr,BFAM_BUFSIZ,"_%05d",s->sub_ix[d]);
    strncat(suffix,tmpstr,len);
  }
}

bfam_subdomain_sbp_t*
bfam_subdomain_sbp_new(const bfam_locidx_t     id,
                            const bfam_locidx_t* sub_ix,
                            const bfam_locidx_t* patch_dist,
                            const char             *name,
                            const int               dim,
                            const bfam_gloidx_t    *N,
                            const bfam_locidx_t    *Nl,
                            const bfam_locidx_t    *buf_sz,
                            const bfam_gloidx_t    *gx,
                            const bfam_long_real_t *c_x,
                            const bfam_long_real_t *c_y,
                            const bfam_long_real_t *c_z)
{
  bfam_subdomain_sbp_t *newSub = bfam_malloc(sizeof(bfam_subdomain_sbp_t));
  bfam_subdomain_sbp_init(newSub,id,sub_ix,patch_dist,name,dim,N,Nl,buf_sz,gx,
                          c_x,c_y,c_z);
  return newSub;
}

void bfam_subdomain_sbp_vtk_vts_piece (struct bfam_subdomain *subdomain,
                                       FILE *file,
                                       const char **scalars,
                                       const char **vectors,
                                       const char **components,
                                       int writeBinary,
                                       int writeCompressed,
                                       int rank)
{
  bfam_subdomain_sbp_t* s = (bfam_subdomain_sbp_t*) subdomain;

  const char *format;

  if(writeBinary)
    format = "binary";
  else
    format = "ascii";


  bfam_locidx_t Nl[3] = {0,0,0};
  bfam_locidx_t buf_sz[6] = {0,0,0,0,0,0};
  bfam_locidx_t Nt[3] = {0,0,0};
  bfam_locidx_t gx[3] = {0,0,0};
  for(int d = 0;d < s->dim; d++)
  {
    Nl[d] = s->Nl[d];
    buf_sz[2*d  ] = s->buf_sz[2*d  ];
    buf_sz[2*d+1] = s->buf_sz[2*d+1];
    Nt[d] = buf_sz[2*d] + Nl[d] + buf_sz[2*d+1];
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

  bfam_locidx_t Ntot = 1;
  for(int i = 0; i < s->dim; i++) Ntot *= (s->Nl[i]+1);
  size_t sz = Ntot*sizeof(bfam_real_t);
  bfam_real_t *storage0 = bfam_malloc_aligned(sz);
  bfam_real_t *storage1 = bfam_malloc_aligned(sz);
  bfam_real_t *storage2 = bfam_malloc_aligned(sz);

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
            &x[i+buf_sz[0] + (Nt[0]+1)*(j+buf_sz[2]
              + (k+buf_sz[4])*(Nt[1]+1))],
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
            &y[i+buf_sz[0] + (Nt[0]+1)*(j+buf_sz[2]
              + (k+buf_sz[4])*(Nt[1]+1))],
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
            &z[i+buf_sz[0] + (Nt[0]+1)*(j+buf_sz[2]
              + (k+buf_sz[4])*(Nt[1]+1))],
            (Nl[0]+1)*sizeof(bfam_real_t));
  }
  else
  {
    memset(storage2,0,sz);
  }

  bfam_vtk_write_real_vector_data_array(file, "Position", writeBinary,
      writeCompressed, Ntot, storage0, storage1, storage2);

  fprintf(file, "      </Points>\n");

  /*
   * Cells
   */
  fprintf(file, "      <Cells>\n");
  fprintf(file, "      </Cells>\n");

  /*
   * PointData
   */
  char pointscalars[BFAM_BUFSIZ];
  bfam_util_strcsl(pointscalars, scalars);

  char pointvectors[BFAM_BUFSIZ];
  bfam_util_strcsl(pointvectors, vectors);

  fprintf(file, "      <PointData"
      " Scalars=\"mpirank,subdomain_id,local_subdomain_id,%s\">\n"
      " Vectors=\"%s\">\n", pointscalars,pointvectors);
  bfam_locidx_t *storage_locid =
                              bfam_malloc_aligned(Ntot*sizeof(bfam_locidx_t));

  /* mpi rank */
  fprintf(file, "        <DataArray type=\"%s\" Name=\"mpirank\""
           " format=\"%s\">\n", BFAM_LOCIDX_VTK, format);
  fprintf(file, "          ");
  if(writeBinary)
  {
    for(bfam_locidx_t i = 0; i < Ntot; ++i)
      storage_locid[i] = rank;

    int rval =
      bfam_vtk_write_binary_data(writeCompressed, file, (char*)storage_locid,
          sizeof(bfam_locidx_t)*Ntot);
    if(rval)
      BFAM_WARNING("Error encoding rank");
  }
  else
  {
    for(bfam_locidx_t i = 0;i < Ntot; ++i)
    {
      fprintf(file, " %6jd", (intmax_t)rank);
    }
  }
  fprintf(file, "\n");
  fprintf(file, "        </DataArray>\n");

  /* subdomain id */
  fprintf(file, "        <DataArray type=\"%s\" Name=\"subdomain_id\""
           " format=\"%s\">\n", BFAM_LOCIDX_VTK, format);
  fprintf(file, "          ");
  if(writeBinary)
  {
    for(bfam_locidx_t i = 0; i < Ntot; ++i)
      storage_locid[i] = s->base.id;

    int rval =
      bfam_vtk_write_binary_data(writeCompressed, file, (char*)storage_locid,
          sizeof(bfam_locidx_t)*Ntot);
    if(rval)
      BFAM_WARNING("Error encoding subdomain_id");
  }
  else
  {
    for(bfam_locidx_t i = 0;i < Ntot; ++i)
    {
      fprintf(file, " %6jd", (intmax_t)s->base.id);
    }
  }
  fprintf(file, "\n");
  fprintf(file, "        </DataArray>\n");

  /* local subdomain id */
  fprintf(file, "        <DataArray type=\"%s\" Name=\"local_subdomain_id\""
           " format=\"%s\">\n", BFAM_LOCIDX_VTK, format);
  fprintf(file, "          ");
  if(writeBinary)
  {
    for(bfam_locidx_t i = 0; i < Ntot; ++i)
      storage_locid[i] = s->patch_id;

    int rval =
      bfam_vtk_write_binary_data(writeCompressed, file, (char*)storage_locid,
          sizeof(bfam_locidx_t)*Ntot);
    if(rval)
      BFAM_WARNING("Error encoding local subdomain id");
  }
  else
  {
    for(bfam_locidx_t i = 0;i < Ntot; ++i)
    {
      fprintf(file, " %6jd", (intmax_t)s->patch_id);
    }
  }
  fprintf(file, "\n");
  fprintf(file, "        </DataArray>\n");

  if(scalars)
  {
    for(size_t s = 0; scalars[s]; ++s)
    {
      bfam_real_t *sdata = bfam_dictionary_get_value_ptr(&subdomain->fields,
          scalars[s]);
      BFAM_ABORT_IF(sdata == NULL, "VTK: Field %s not in subdomain %s",
          scalars[s], subdomain->name);
      int i = 0;
      for(int k = 0; k <= Nl[2]; k++)
        for(int j = 0; j <= Nl[1]; j++)
          memcpy(&storage0[i + (Nl[0]+1)*(j + k*(Nl[1]+1))],
              &sdata[i+buf_sz[0] + (Nt[0]+1)*(j+buf_sz[2]
                + (k+buf_sz[4])*(Nt[1]+1))],
              (Nl[0]+1)*sizeof(bfam_real_t));
      bfam_vtk_write_real_scalar_data_array(file, scalars[s],
          writeBinary, writeCompressed, Ntot, storage0);
    }
  }

  if(vectors)
  {
    for(size_t v = 0; vectors[v]; ++v)
    {

      bfam_real_t *v1 =
        bfam_dictionary_get_value_ptr(&subdomain->fields, components[3*v+0]);
      bfam_real_t *v2 =
        bfam_dictionary_get_value_ptr(&subdomain->fields, components[3*v+1]);
      bfam_real_t *v3 =
        bfam_dictionary_get_value_ptr(&subdomain->fields, components[3*v+2]);

      if(v1!=NULL)
      {
        int i = 0;
        for(int k = 0; k <= Nl[2]; k++)
          for(int j = 0; j <= Nl[1]; j++)
            memcpy(&storage0[i + (Nl[0]+1)*(j + k*(Nl[1]+1))],
                &v1[i+buf_sz[0] + (Nt[0]+1)*(j+buf_sz[2]
                  + (k+buf_sz[4])*(Nt[1]+1))],
                (Nl[0]+1)*sizeof(bfam_real_t));
      }
      else
      {
        memset(storage0,0,sz);
      }

      if(v2!=NULL)
      {
        int i = 0;
        for(int k = 0; k <= Nl[2]; k++)
          for(int j = 0; j <= Nl[1]; j++)
            memcpy(&storage1[i + (Nl[0]+1)*(j + k*(Nl[1]+1))],
                &v2[i+buf_sz[0] + (Nt[0]+1)*(j+buf_sz[2]
                  + (k+buf_sz[4])*(Nt[1]+1))],
                (Nl[0]+1)*sizeof(bfam_real_t));
      }
      else
      {
        memset(storage1,0,sz);
      }

      if(v3!=NULL)
      {
        int i = 0;
        for(int k = 0; k <= Nl[2]; k++)
          for(int j = 0; j <= Nl[1]; j++)
            memcpy(&storage2[i + (Nl[0]+1)*(j + k*(Nl[1]+1))],
                &v3[i+buf_sz[0] + (Nt[0]+1)*(j+buf_sz[2]
                  + (k+buf_sz[4])*(Nt[1]+1))],
                (Nl[0]+1)*sizeof(bfam_real_t));
      }
      else
      {
        memset(storage2,0,sz);
      }

      bfam_vtk_write_real_vector_data_array(file, vectors[v],
          writeBinary, writeCompressed, Ntot, storage0, storage1, storage2);
    }
  }

  fprintf(file, "      </PointData>\n");

  fprintf(file, "    </Piece>\n");

  /* close StructuredGrid */
  fprintf(file, "  </StructuredGrid>\n");

  bfam_free_aligned(storage_locid);
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

  size_t fieldSize = 1;
  for(int d = 0; d < s->dim; d++)
    fieldSize *= (s->Nl[d]+1+s->buf_sz[2*d]+s->buf_sz[2*d+1]);

  bfam_real_t *field = bfam_malloc_aligned(fieldSize*sizeof(bfam_real_t));

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

  size_t fieldSize = 1;
  for(int d = 0; d < s->dim; d++)
    fieldSize *= (s->Nl[d]+1+s->buf_sz[2*d]+s->buf_sz[2*d+1]);

  bfam_real_t *restrict x =
    bfam_dictionary_get_value_ptr(&subdomain->fields, "_grid_x");
  bfam_real_t *restrict y =
    bfam_dictionary_get_value_ptr(&subdomain->fields, "_grid_y");
  bfam_real_t *restrict z =
    bfam_dictionary_get_value_ptr(&subdomain->fields, "_grid_z");

  init_field(fieldSize, name, time, x, y, z, subdomain, arg, field);
}

void
bfam_subdomain_sbp_init(bfam_subdomain_sbp_t *subdomain,
                            const bfam_locidx_t     id,
                            const bfam_locidx_t* sub_ix,
                            const bfam_locidx_t* patch_dist,
                            const char             *name,
                            const int               dim,
                            const bfam_gloidx_t    *N,
                            const bfam_locidx_t    *Nl,
                            const bfam_locidx_t    *buf_sz,
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

  subdomain->sub_ix = bfam_malloc(dim*sizeof(bfam_locidx_t));
  subdomain->patch_dist  = bfam_malloc(dim*sizeof(bfam_locidx_t));
  subdomain->patch_id = 0;
  subdomain->num_id = 1;
  subdomain->dim = dim;

  subdomain->N  = bfam_malloc(  dim*sizeof(bfam_gloidx_t));
  subdomain->Nl = bfam_malloc(  dim*sizeof(bfam_locidx_t));
  subdomain->buf_sz = bfam_malloc(2*dim*sizeof(bfam_locidx_t));
  subdomain->gx = bfam_malloc(  dim*sizeof(bfam_gloidx_t));
  bfam_gloidx_t gxb[dim];

  bfam_locidx_t Nltmp[dim];

  bfam_locidx_t sz = 1;
  for(int d = 0; d < dim; d++)
  {
    subdomain->N [  d  ] = N [  d  ];
    subdomain->Nl[  d  ] = Nl[  d  ];
    subdomain->buf_sz[2*d  ] = buf_sz[2*d  ];
    subdomain->buf_sz[2*d+1] = buf_sz[2*d+1];
    subdomain->gx[  d  ] = gx[  d  ];
    gxb[d] = gx[d]-buf_sz[2*d];

    Nltmp[d] = (buf_sz[2*d] + Nl[d] + buf_sz[2*d+1]);
    sz *= (Nltmp[d]+1);

    subdomain->sub_ix[d] = sub_ix[d];
    subdomain->patch_dist [d] = patch_dist[d];
    subdomain->patch_id += sub_ix[d]*subdomain->num_id;
    subdomain->num_id *= (patch_dist[d]+1);

    BFAM_ABORT_IF(N [d] < 0 || Nl[d] < 0 ||
                  buf_sz[2*d] < 0 || buf_sz[2*d+1] < 0 ||
                  gx[d] < 0 || gx[d] > N[d] || gx[d]+Nl[d] > N[d],
                  "%s: problem with a dimension %d", d, subdomain->base.name);
  }

  bfam_real_t *x = NULL;
  if(c_x != NULL)
  {
    int rval = bfam_subdomain_sbp_field_add(&subdomain->base, "_grid_x");
    BFAM_ASSERT(rval==2);
    x = (bfam_real_t*) bfam_dictionary_get_value_ptr(&subdomain->base.fields,
        "_grid_x");
  }

  bfam_real_t *y = NULL;
  if(c_y != NULL)
  {
    BFAM_ASSERT(c_x != NULL);
    int rval = bfam_subdomain_sbp_field_add(&subdomain->base, "_grid_y");
    BFAM_ASSERT(rval==2);
    y = (bfam_real_t*) bfam_dictionary_get_value_ptr(&subdomain->base.fields,
        "_grid_y");
  }

  bfam_real_t *z = NULL;
  if(c_z != NULL)
  {
    BFAM_ASSERT(c_y != NULL);
    int rval = bfam_subdomain_sbp_field_add(&subdomain->base, "_grid_z");
    BFAM_ASSERT(rval==2);
    z = (bfam_real_t*) bfam_dictionary_get_value_ptr(&subdomain->base.fields,
        "_grid_z");
  }

  bfam_util_linear_blend(x,y,z,dim,N,Nltmp,gxb,c_x,c_y,c_z);
}

static int
bfam_subdomain_sbp_free_fields(const char * key, void *val,
    void *arg)
{
  if(val != NULL) bfam_free_aligned(val);

  return 1;
}

void
bfam_subdomain_sbp_free(bfam_subdomain_t *subdomain)
{
  bfam_subdomain_sbp_t *sub = (bfam_subdomain_sbp_t*) subdomain;

  bfam_dictionary_allprefixed_ptr(&sub->base.fields,"",
      &bfam_subdomain_sbp_free_fields,NULL);
  bfam_dictionary_allprefixed_ptr(&sub->base.fields_p,"",
      &bfam_subdomain_sbp_free_fields,NULL);
  bfam_dictionary_allprefixed_ptr(&sub->base.fields_m,"",
      &bfam_subdomain_sbp_free_fields,NULL);
  bfam_dictionary_allprefixed_ptr(&sub->base.fields_face,"",
      &bfam_subdomain_sbp_free_fields,NULL);

  bfam_subdomain_free(subdomain);

  bfam_free(sub->N );
  bfam_free(sub->Nl);
  bfam_free(sub->buf_sz);
  bfam_free(sub->gx);
  bfam_free(sub->sub_ix);
  bfam_free(sub->patch_dist);
}

bfam_subdomain_sbp_intra_glue_t*
bfam_subdomain_sbp_intra_glue_new(const bfam_locidx_t              id,
                                 const char                      *name,
                                 const bfam_locidx_t              rank_m,
                                 const bfam_locidx_t              rank_p,
                                 bfam_subdomain_sbp_t            *sub_m,
                                 const int                        face)
{
  bfam_subdomain_sbp_intra_glue_t* sub =
    bfam_malloc(sizeof(bfam_subdomain_sbp_intra_glue_t));

  bfam_subdomain_sbp_intra_glue_init(sub,id,name,rank_m,rank_p,sub_m,face);

  return sub;
}

static void
bfam_subdomain_sbp_intra_glue_comm_info(bfam_subdomain_t *thisSubdomain,
    int *rank, bfam_locidx_t *sort, int sort_num,
    size_t *send_sz, size_t *recv_sz)
{
  BFAM_ASSERT(sort_num > 5);
  bfam_subdomain_sbp_intra_glue_t *sub =
    (bfam_subdomain_sbp_intra_glue_t*) thisSubdomain;
  bfam_subdomain_sbp_t *sub_m = sub->sub_m;

  *rank     = sub->rank_p;
  sort[0] = sub->id_m;
  sort[1] = sub->id_p;
  sort[2] = sub->patch_id_m;
  sort[3] = sub->patch_id_p;
  sort[4] = sub->fce_m;
  sort[5] = sub->fce_p;

  size_t s_num = 1;
  size_t r_num = 1;
  for(int d = 0;d < sub_m->dim;d++)
  {
    s_num *= (sub->send_ix[2*d+1]+1-sub->send_ix[2*d]);
    r_num *= (sub->recv_ix[2*d+1]+1-sub->recv_ix[2*d]);
  }

#ifdef BFAM_DEBUG
  if(sub_m->dim == 2)
  {
    BFAM_LDEBUG(" r_m %3d r_p %3d Nl[%3d %3d] buf_sz[%3d %3d %3d %3d]"
        " send_ix[%3d %3d %3d %3d] size %3d",
        sub->rank_m,sub->rank_p,
        sub_m->Nl[0],sub_m->Nl[1],
        sub_m->buf_sz[0],sub_m->buf_sz[1],sub_m->buf_sz[2],sub_m->buf_sz[3],
        sub->send_ix[0],sub->send_ix[1],sub->send_ix[2],sub->send_ix[3],
        s_num);
  }
  if(sub_m->dim == 3)
  {
    BFAM_LDEBUG(" r_m %3d r_p %3d Nl[%3d %3d %3d]"
        " buf_sz[%3d %3d %3d %3d %3d %3d]"
        " send_ix[%3d %3d %3d %3d %3d %3d] size %3d",
        sub->rank_m,sub->rank_p,
        sub_m->Nl[0],sub_m->Nl[1],sub_m->Nl[2],
        sub_m->buf_sz[0],sub_m->buf_sz[1],sub_m->buf_sz[2],
        sub_m->buf_sz[3],sub_m->buf_sz[4],sub_m->buf_sz[5],
        sub->send_ix[0],sub->send_ix[1],sub->send_ix[2],
        sub->send_ix[3],sub->send_ix[4],sub->send_ix[5],
        r_num);
  }
#endif


  *send_sz  = sub->base.fields_m.num_entries*s_num*sizeof(bfam_real_t);
  *recv_sz  = sub->base.fields_p.num_entries*r_num*sizeof(bfam_real_t);
}

static int
bfam_subdomain_sbp_intra_glue_field_minus_add(bfam_subdomain_t *subdomain,
    const char *name)
{
  bfam_subdomain_sbp_intra_glue_t *s =
    (bfam_subdomain_sbp_intra_glue_t*) subdomain;

  if(bfam_dictionary_get_value_ptr(&s->base.fields_m,name))
    return 1;

  int rval = bfam_dictionary_insert_ptr(&s->base.fields_m, name, NULL);

  BFAM_ASSERT(rval != 1);

  return rval;
}


static int
bfam_subdomain_sbp_intra_glue_field_plus_add(bfam_subdomain_t *subdomain,
    const char *name)
{
  bfam_subdomain_sbp_intra_glue_t *s =
    (bfam_subdomain_sbp_intra_glue_t*) subdomain;

  if(bfam_dictionary_get_value_ptr(&s->base.fields_p,name))
    return 1;

  int rval = bfam_dictionary_insert_ptr(&s->base.fields_p, name, NULL);

  BFAM_ASSERT(rval != 1);

  return rval;
}

typedef struct bfam_subdomain_sbp_intra_get_put_data
{
  bfam_subdomain_sbp_intra_glue_t *sub;
  bfam_real_t *buffer;
  size_t size;
  size_t field;
} bfam_subdomain_sbp_intra_get_put_data_t;

static int
bfam_subdomain_sbp_intra_glue_get_fields_m(const char * key, void *val,
    void *arg)
{
  bfam_subdomain_sbp_intra_get_put_data_t *data =
    (bfam_subdomain_sbp_intra_get_put_data_t*) arg;

  bfam_subdomain_sbp_intra_glue_t *sub = data->sub;

  const size_t buffer_offset = data->field * sub->field_size_m;
  bfam_real_t *restrict send_field = data->buffer + buffer_offset;

  bfam_subdomain_sbp_t *sub_m = data->sub->sub_m;

  const bfam_real_t *restrict sub_m_field =
    bfam_dictionary_get_value_ptr(&sub->sub_m->base.fields, key);

  BFAM_ASSERT(sub_m->dim < 4);

  BFAM_ASSERT( send_field != NULL);
  BFAM_ASSERT(sub_m_field != NULL);

  BFAM_ASSUME_ALIGNED(sub_m_field, 32);

  /* length of line being copied */
  size_t line_size =  sub->send_ix[1]+1-sub->send_ix[0];
  if(sub_m->dim == 3)
  {
    bfam_locidx_t num_x = sub_m->buf_sz[0] + sub_m->Nl[0]+1 + sub_m->buf_sz[1];
    bfam_locidx_t num_y = sub_m->buf_sz[2] + sub_m->Nl[1]+1 + sub_m->buf_sz[3];
    const int i = sub->send_ix[0];
    for(int j = sub->send_ix[2];j < sub->send_ix[3]+1;j++)
    {
      for(int k = sub->send_ix[4];k < sub->send_ix[5]+1;k++)
      {
        memcpy(send_field,sub_m_field+(i+num_x*(j+num_y*k)),
            line_size * sizeof(bfam_real_t));
        send_field += line_size;
      }
    }
  }
  else if(sub_m->dim == 2)
  {
    bfam_locidx_t num_x = sub_m->buf_sz[0] + sub_m->Nl[0]+1 + sub_m->buf_sz[1];
    const int i = sub->send_ix[0];
    for(int j = sub->send_ix[2];j < sub->send_ix[3]+1;j++)
    {
      memcpy(send_field,sub_m_field+(i+num_x*j),
          line_size * sizeof(bfam_real_t));
      send_field += line_size;
    }
  }
  else if(sub_m->dim == 1)
  {
    const int i = sub->send_ix[0];
    memcpy(send_field,sub_m_field+i,line_size * sizeof(bfam_real_t));
  }
  else BFAM_ABORT("must have dim 1, 2, or 3");

  data->field++;
  return 0;
}


void
bfam_subdomain_sbp_intra_glue_put_send_buffer(bfam_subdomain_t *thisSubdomain,
    void *buffer, size_t send_sz)
{
  bfam_subdomain_sbp_intra_get_put_data_t data;

  data.sub    = (bfam_subdomain_sbp_intra_glue_t*) thisSubdomain;
  data.buffer = (bfam_real_t*) buffer;
  data.size   = send_sz;
  data.field  = 0;

  /*
   * Fill fields_m and the send buffer from sub_m.
   */
  bfam_dictionary_allprefixed_ptr(&data.sub->base.fields_m, "",
      &bfam_subdomain_sbp_intra_glue_get_fields_m, &data);
}

static int
bfam_subdomain_sbp_intra_glue_put_fields_p(const char * key, void *val,
    void *arg)
{
  bfam_subdomain_sbp_intra_get_put_data_t *data =
    (bfam_subdomain_sbp_intra_get_put_data_t*) arg;

  bfam_subdomain_sbp_intra_glue_t *sub = data->sub;

  const size_t buffer_offset = data->field * sub->field_size_p;
  bfam_real_t *restrict recv_field = data->buffer + buffer_offset;

  bfam_subdomain_sbp_t *sub_m = data->sub->sub_m;

  bfam_real_t *restrict sub_m_field =
    bfam_dictionary_get_value_ptr(&sub->sub_m->base.fields, key);

  BFAM_ASSERT(sub_m->dim < 4);

  BFAM_ASSERT( recv_field != NULL);
  BFAM_ASSERT(sub_m_field != NULL);

  BFAM_ASSUME_ALIGNED(sub_m_field, 32);

  /* length of line being copied */
  size_t line_size =  sub->recv_ix[1]+1-sub->recv_ix[0];
  if(sub_m->dim == 3)
  {
    bfam_locidx_t num_x = sub_m->buf_sz[0] + sub_m->Nl[0]+1 + sub_m->buf_sz[1];
    bfam_locidx_t num_y = sub_m->buf_sz[2] + sub_m->Nl[1]+1 + sub_m->buf_sz[3];
    const int i = sub->recv_ix[0];
    for(int j = sub->recv_ix[2];j < sub->recv_ix[3]+1;j++)
    {
      for(int k = sub->recv_ix[4];k < sub->recv_ix[5]+1;k++)
      {
        memcpy(sub_m_field+(i+num_x*(j+num_y*k)),recv_field,
            line_size * sizeof(bfam_real_t));
        recv_field += line_size;
      }
    }
  }
  else if(sub_m->dim == 2)
  {
    bfam_locidx_t num_x = sub_m->buf_sz[0] + sub_m->Nl[0]+1 + sub_m->buf_sz[1];
    const int i = sub->recv_ix[0];
    for(int j = sub->recv_ix[2];j < sub->recv_ix[3]+1;j++)
    {
      memcpy(sub_m_field+(i+num_x*j),recv_field,
          line_size * sizeof(bfam_real_t));
      recv_field += line_size;
    }
  }
  else if(sub_m->dim == 1)
  {
    const int i = sub->recv_ix[0];
    memcpy(sub_m_field+i,recv_field,line_size * sizeof(bfam_real_t));
  }
  else BFAM_ABORT("must have dim 1, 2, or 3");

  data->field++;
  return 0;
}

void
bfam_subdomain_sbp_intra_glue_get_recv_buffer(bfam_subdomain_t *thisSubdomain,
    void *buffer, size_t recv_sz)
{
  bfam_subdomain_sbp_intra_get_put_data_t data;

  data.sub    = (bfam_subdomain_sbp_intra_glue_t*) thisSubdomain;
  data.buffer = (bfam_real_t*) buffer;
  data.size   = recv_sz;
  data.field  = 0;

  /*
   * Fill fields_p from receive buffer
   */
  bfam_dictionary_allprefixed_ptr(&data.sub->base.fields_m, "",
      &bfam_subdomain_sbp_intra_glue_put_fields_p, &data);
}


void
bfam_subdomain_sbp_intra_glue_init(bfam_subdomain_sbp_intra_glue_t* sub,
                                 const bfam_locidx_t              id,
                                 const char                      *name,
                                 const bfam_locidx_t              rank_m,
                                 const bfam_locidx_t              rank_p,
                                 bfam_subdomain_sbp_t            *sub_m,
                                 const int                        face)
{
  bfam_subdomain_init(&sub->base, id, name);
  bfam_subdomain_add_tag(&sub->base, "_subdomain_sbp");
  bfam_subdomain_add_tag(&sub->base, "_subdomain_sbp_intra_glue");

  sub->base.glue_comm_info = bfam_subdomain_sbp_intra_glue_comm_info;
  sub->base.glue_put_send_buffer =
    bfam_subdomain_sbp_intra_glue_put_send_buffer;
  sub->base.glue_get_recv_buffer =
    bfam_subdomain_sbp_intra_glue_get_recv_buffer;
  sub->base.field_minus_add = bfam_subdomain_sbp_intra_glue_field_minus_add;
  sub->base.field_plus_add  = bfam_subdomain_sbp_intra_glue_field_plus_add;
  sub->base.free = bfam_subdomain_sbp_intra_glue_free;

  sub->rank_m = rank_m;
  sub->rank_p = rank_p;
  sub->sub_m  = sub_m;
  sub->face   = face;

  sub->id_m  = sub_m->base.id;
  sub->patch_id_m = sub_m->patch_id;
  sub->fce_m = face;

  /* initialize all the indices to full size */
  sub->send_ix = bfam_malloc(2*sub_m->dim*sizeof(bfam_locidx_t));
  sub->recv_ix = bfam_malloc(2*sub_m->dim*sizeof(bfam_locidx_t));
  sub->field_size_m = 1;
  sub->field_size_p = 1;
  for(int d = 0;d < sub_m->dim;d++)
  {
    sub->send_ix[2*d  ] = sub_m->buf_sz[2*d];
    sub->send_ix[2*d+1] = sub_m->buf_sz[2*d]+sub_m->Nl[d];
    sub->recv_ix[2*d  ] = sub_m->buf_sz[2*d];
    sub->recv_ix[2*d+1] = sub_m->buf_sz[2*d]+sub_m->Nl[d];

    sub->field_size_m *= (sub->send_ix[2*d+1]+1-sub->send_ix[2*d]);
    sub->field_size_p *= (sub->recv_ix[2*d+1]+1-sub->recv_ix[2*d]);
  }

  sub->id_p  = sub_m->base.id; /* assume we number in the same order! */
  switch(face)
  {
    case 0 :
      sub->patch_id_p  = sub_m->patch_id-1;

      sub->field_size_m /= (sub->send_ix[1]+1-sub->send_ix[0]);
      sub->field_size_p /= (sub->recv_ix[1]+1-sub->recv_ix[0]);
      sub->field_size_m *= sub_m->buf_sz[0];
      sub->field_size_p *= sub_m->buf_sz[0];

      sub->recv_ix[0] = 0;
      sub->recv_ix[1] = sub_m->buf_sz[0]-1 + sub->recv_ix[0];
      sub->send_ix[0] = sub_m->buf_sz[0]   + sub->recv_ix[0] ;
      sub->send_ix[1] = sub_m->buf_sz[0]   + sub->recv_ix[1] ;
      break;
    case 1 :
      sub->patch_id_p  = sub_m->patch_id+1;

      sub->field_size_m /= (sub->send_ix[1]+1-sub->send_ix[0]);
      sub->field_size_p /= (sub->recv_ix[1]+1-sub->recv_ix[0]);
      sub->field_size_m *= sub_m->buf_sz[1];
      sub->field_size_p *= sub_m->buf_sz[1];

      sub->recv_ix[0] = sub_m->buf_sz[0] + sub_m->Nl[0]+1;
      sub->recv_ix[1] = sub_m->buf_sz[1]-1 + sub->recv_ix[0];
      sub->send_ix[0] = sub->recv_ix[0] - sub_m->buf_sz[1];
      sub->send_ix[1] = sub->recv_ix[1] - sub_m->buf_sz[1];
      break;
    case 2 :
      sub->patch_id_p  = sub_m->patch_id-(sub_m->patch_dist[0]+1);

      sub->field_size_m /= (sub->send_ix[3]+1-sub->send_ix[2]);
      sub->field_size_p /= (sub->recv_ix[3]+1-sub->recv_ix[2]);
      sub->field_size_m *= sub_m->buf_sz[2];
      sub->field_size_p *= sub_m->buf_sz[2];

      sub->recv_ix[2] = 0;
      sub->recv_ix[3] = sub_m->buf_sz[2]-1 + sub->recv_ix[2];
      sub->send_ix[2] = sub_m->buf_sz[2]   + sub->recv_ix[2] ;
      sub->send_ix[3] = sub_m->buf_sz[2]   + sub->recv_ix[3] ;
      break;
    case 3 :
      sub->patch_id_p  = sub_m->patch_id+(sub_m->patch_dist[0]+1);

      sub->field_size_m /= (sub->send_ix[3]+1-sub->send_ix[2]);
      sub->field_size_p /= (sub->recv_ix[3]+1-sub->recv_ix[2]);
      sub->field_size_m *= sub_m->buf_sz[3];
      sub->field_size_p *= sub_m->buf_sz[3];

      sub->recv_ix[2] = sub_m->buf_sz[2] + sub_m->Nl[1]+1;
      sub->recv_ix[3] = sub_m->buf_sz[3]-1 + sub->recv_ix[2];
      sub->send_ix[2] = sub->recv_ix[2] - sub_m->buf_sz[3];
      sub->send_ix[3] = sub->recv_ix[3] - sub_m->buf_sz[3];
      break;
    case 4 :
      sub->patch_id_p  =
        sub_m->patch_id-(sub_m->patch_dist[0]+1)*(sub_m->patch_dist[1]+1);

      sub->field_size_m /= (sub->send_ix[5]+1-sub->send_ix[4]);
      sub->field_size_p /= (sub->recv_ix[5]+1-sub->recv_ix[4]);
      sub->field_size_m *= sub_m->buf_sz[4];
      sub->field_size_p *= sub_m->buf_sz[4];

      sub->recv_ix[4] = 0;
      sub->recv_ix[5] = sub_m->buf_sz[4]-1 + sub->recv_ix[4];
      sub->send_ix[4] = sub_m->buf_sz[4]   + sub->recv_ix[4] ;
      sub->send_ix[5] = sub_m->buf_sz[4]   + sub->recv_ix[5] ;

      break;
    case 5 :
      sub->patch_id_p  =
        sub_m->patch_id+(sub_m->patch_dist[0]+1)*(sub_m->patch_dist[1]+1);

      sub->field_size_m /= (sub->send_ix[5]+1-sub->send_ix[4]);
      sub->field_size_p /= (sub->recv_ix[5]+1-sub->recv_ix[4]);
      sub->field_size_m *= sub_m->buf_sz[5];
      sub->field_size_p *= sub_m->buf_sz[5];

      sub->recv_ix[4] = sub_m->buf_sz[4] + sub_m->Nl[2]+1;
      sub->recv_ix[5] = sub_m->buf_sz[5]-1 + sub->recv_ix[4];
      sub->send_ix[4] = sub->recv_ix[4] - sub_m->buf_sz[5];
      sub->send_ix[5] = sub->recv_ix[5] - sub_m->buf_sz[5];
      break;
  }

  if(0 == face%2)
    sub->fce_p = face+1;
  else
    sub->fce_p = face-1;
}

void
bfam_subdomain_sbp_intra_glue_free(bfam_subdomain_t *thisSubdomain)
{
  bfam_subdomain_sbp_intra_glue_t* sub =
    (bfam_subdomain_sbp_intra_glue_t*) thisSubdomain;

  bfam_dictionary_allprefixed_ptr(&sub->base.fields,"",
      &bfam_subdomain_sbp_free_fields,NULL);
  bfam_dictionary_allprefixed_ptr(&sub->base.fields_p,"",
      &bfam_subdomain_sbp_free_fields,NULL);
  bfam_dictionary_allprefixed_ptr(&sub->base.fields_m,"",
      &bfam_subdomain_sbp_free_fields,NULL);
  bfam_dictionary_allprefixed_ptr(&sub->base.fields_face,"",
      &bfam_subdomain_sbp_free_fields,NULL);

  bfam_subdomain_free(thisSubdomain);

  bfam_free(sub->send_ix);
  bfam_free(sub->recv_ix);
}

bfam_subdomain_sbp_inter_glue_t*
bfam_subdomain_sbp_inter_glue_new(const bfam_locidx_t              id,
                                 const char                      *name,
                                 const bfam_locidx_t              rank_m,
                                 const bfam_locidx_t              rank_p,
                                 bfam_subdomain_sbp_t            *sub_m,
                                 const bfam_gloidx_t             *fc_gx_rg,
                                 const int                        face,
                                 const int                        orient,
                                 bfam_locidx_t                    id_p,
                                 bfam_locidx_t                    patch_id_p,
                                 bfam_locidx_t                    face_p)
{
  bfam_subdomain_sbp_inter_glue_t *newSub =
    bfam_malloc(sizeof(bfam_subdomain_sbp_inter_glue_t));
  bfam_subdomain_sbp_inter_glue_init(newSub,id,name,rank_m,rank_p,
      sub_m,fc_gx_rg,face, orient,id_p,patch_id_p,face_p);
  return newSub;
}

static int
bfam_subdomain_sbp_inter_glue_field_minus_add(bfam_subdomain_t *subdomain,
    const char *name)
{
  bfam_subdomain_sbp_inter_glue_t *s =
    (bfam_subdomain_sbp_inter_glue_t*) subdomain;

  if(bfam_dictionary_get_value_ptr(&s->base.fields_m,name))
    return 1;

  bfam_real_t *field =
    bfam_malloc_aligned(s->field_size_m*sizeof(bfam_real_t));
  int rval = bfam_dictionary_insert_ptr(&s->base.fields_m, name, field);

  BFAM_ASSERT(rval != 1);

  return rval;
}


static int
bfam_subdomain_sbp_inter_glue_field_plus_add(bfam_subdomain_t *subdomain,
    const char *name)
{
  bfam_subdomain_sbp_inter_glue_t *s =
    (bfam_subdomain_sbp_inter_glue_t*) subdomain;

  if(bfam_dictionary_get_value_ptr(&s->base.fields_p,name))
    return 1;

  bfam_real_t *field = bfam_malloc_aligned(s->field_size_p*sizeof(bfam_real_t));
  int rval = bfam_dictionary_insert_ptr(&s->base.fields_p, name, field);

  BFAM_ASSERT(rval != 1);

  return rval;
}

static void
bfam_subdomain_sbp_inter_glue_comm_info(bfam_subdomain_t *thisSubdomain,
    int *rank, bfam_locidx_t *sort, int sort_num,
    size_t *send_sz, size_t *recv_sz)
{
  BFAM_ASSERT(sort_num > 5);
  bfam_subdomain_sbp_inter_glue_t *sub =
    (bfam_subdomain_sbp_inter_glue_t*) thisSubdomain;

  *rank     = sub->rank_p;
  sort[0] = sub->id_m;
  sort[1] = sub->id_p;
  sort[2] = sub->patch_id_m;
  sort[3] = sub->patch_id_p;
  sort[4] = sub->fce_m;
  sort[5] = sub->fce_p;

  *send_sz  =
    sub->base.fields_m.num_entries*sub->field_size_m*sizeof(bfam_real_t);
  *recv_sz  =
    sub->base.fields_p.num_entries*sub->field_size_p*sizeof(bfam_real_t);
}

typedef struct bfam_subdomain_sbp_inter_get_put_data
{
  bfam_subdomain_sbp_inter_glue_t *sub;
  bfam_real_t *buffer;
  size_t size;
  size_t field;
} bfam_subdomain_sbp_inter_get_put_data_t;

static int
bfam_subdomain_sbp_inter_glue_get_fields_m(const char * key, void *val,
    void *arg)
{
  bfam_subdomain_sbp_inter_get_put_data_t *data =
    (bfam_subdomain_sbp_inter_get_put_data_t*) arg;

  bfam_subdomain_sbp_inter_glue_t *sub = data->sub;

  const size_t buffer_offset = data->field * sub->field_size_m;
  bfam_real_t *restrict send_field = data->buffer + buffer_offset;

  bfam_subdomain_sbp_t *sub_m = data->sub->sub_m;

  const bfam_real_t *restrict sub_m_field =
    bfam_dictionary_get_value_ptr(&sub->sub_m->base.fields, key);
  bfam_real_t *restrict glue_field = val;

  BFAM_ASSERT(sub_m->dim < 4);

  BFAM_ASSERT( glue_field != NULL);
  BFAM_ASSERT( send_field != NULL);
  BFAM_ASSERT(sub_m_field != NULL);

  BFAM_ASSUME_ALIGNED(sub_m_field, 32);
  BFAM_ASSUME_ALIGNED( glue_field, 32);

  /* length of line being copied */
  size_t line_size = sub->vl_ix_rg[1]+1-sub->vl_ix_rg[0];
  if(sub_m->dim == 3)
  {
    BFAM_ABORT("3D not implemented");
  }
  else if(sub_m->dim == 2)
  {
    bfam_locidx_t num_x = sub_m->buf_sz[0] + sub_m->Nl[0]+1 + sub_m->buf_sz[1];
    const int i = sub->vl_ix_rg[0];
    for(int j = sub->vl_ix_rg[2],count=0;j < sub->vl_ix_rg[3]+1;j++,count++)
      memcpy(glue_field+line_size*count,sub_m_field+(i+num_x*j),
          line_size*sizeof(bfam_real_t));
    if(sub->orient==0)
    {
      memcpy(send_field,glue_field,sub->field_size_m*sizeof(bfam_real_t));
    }
    else
    {
      for(bfam_locidx_t itor=0,jtor=sub->field_size_m-1;
          itor < (bfam_locidx_t)sub->field_size_m; itor++,jtor--)
        send_field[itor] = glue_field[jtor];
    }
  }
  else if(sub_m->dim == 1)
  {
    BFAM_ABORT("1D not implemented");
  }
  else BFAM_ABORT("must have dim 1, 2, or 3");

  data->field++;
  return 0;
}

void
bfam_subdomain_sbp_inter_glue_put_send_buffer(bfam_subdomain_t *thisSubdomain,
    void *buffer, size_t send_sz)
{
  bfam_subdomain_sbp_inter_get_put_data_t data;

  data.sub    = (bfam_subdomain_sbp_inter_glue_t*) thisSubdomain;
  data.buffer = (bfam_real_t*) buffer;
  data.size   = send_sz;
  data.field  = 0;

  /*
   * Fill fields_m and the send buffer from sub_m.
   */
  bfam_dictionary_allprefixed_ptr(&data.sub->base.fields_m, "",
      &bfam_subdomain_sbp_inter_glue_get_fields_m, &data);
}

static int
bfam_subdomain_sbp_inter_glue_put_fields_p(const char * key, void *val,
    void *arg)
{
  bfam_subdomain_sbp_inter_get_put_data_t *data =
    (bfam_subdomain_sbp_inter_get_put_data_t*) arg;

  bfam_subdomain_sbp_inter_glue_t *sub = data->sub;

  const size_t buffer_offset = data->field * sub->field_size_p;
  bfam_real_t *restrict recv_field = data->buffer + buffer_offset;

  bfam_real_t *restrict glue_field = val;

  BFAM_ASSUME_ALIGNED(glue_field, 32);

  memcpy(glue_field, recv_field, sub->field_size_p * sizeof(bfam_real_t));

  data->field++;
  return 0;
}

void
bfam_subdomain_sbp_inter_glue_get_recv_buffer(bfam_subdomain_t *thisSubdomain,
    void *buffer, size_t recv_sz)
{
  bfam_subdomain_sbp_inter_get_put_data_t data;

  data.sub    = (bfam_subdomain_sbp_inter_glue_t*) thisSubdomain;
  data.buffer = (bfam_real_t*) buffer;
  data.size   = recv_sz;
  data.field  = 0;

  /*
   * Fill fields_p from the receive buffer
   */
  bfam_dictionary_allprefixed_ptr(&data.sub->base.fields_p, "",
      &bfam_subdomain_sbp_inter_glue_put_fields_p, &data);
}


void
bfam_subdomain_sbp_inter_glue_init(bfam_subdomain_sbp_inter_glue_t* sub,
                                 const bfam_locidx_t              id,
                                 const char                      *name,
                                 const bfam_locidx_t              rank_m,
                                 const bfam_locidx_t              rank_p,
                                 bfam_subdomain_sbp_t            *sub_m,
                                 const bfam_gloidx_t             *fc_gx_rg,
                                 const int                        face,
                                 const int                        orient,
                                 bfam_locidx_t                    id_p,
                                 bfam_locidx_t                    patch_id_p,
                                 bfam_locidx_t                    face_p)
{
  bfam_subdomain_init(&sub->base, id, name);
  bfam_subdomain_add_tag(&sub->base, "_subdomain_sbp");
  bfam_subdomain_add_tag(&sub->base, "_subdomain_sbp_inter_glue");

  sub->base.glue_comm_info = bfam_subdomain_sbp_inter_glue_comm_info;
  sub->base.glue_put_send_buffer =
    bfam_subdomain_sbp_inter_glue_put_send_buffer;
  sub->base.glue_get_recv_buffer =
    bfam_subdomain_sbp_inter_glue_get_recv_buffer;
  sub->base.field_minus_add = bfam_subdomain_sbp_inter_glue_field_minus_add;
  sub->base.field_plus_add  = bfam_subdomain_sbp_inter_glue_field_plus_add;
  sub->base.free = bfam_subdomain_sbp_inter_glue_free;

  sub->rank_m = rank_m;
  sub->rank_p = rank_p;
  sub->sub_m  = sub_m;
  sub->face   = face;
  sub->orient = orient;

  sub->id_m  = sub_m->base.id;
  sub->patch_id_m = sub_m->patch_id;
  sub->fce_m = face;

  sub->id_p  = id_p;
  sub->patch_id_p = patch_id_p;
  sub->fce_p = face_p;

  sub->field_size_m = 1;
  sub->field_size_p = 1;

  sub->fc_gx_rg = bfam_malloc((sub_m->dim-1)*2*sizeof(bfam_gloidx_t));
  sub->vl_ix_rg = bfam_malloc( sub_m->dim   *2*sizeof(bfam_gloidx_t));
  for(int d = 0; d < sub_m->dim-1;d++)
  {
    sub->fc_gx_rg[2*d  ] = fc_gx_rg[2*d  ];
    sub->fc_gx_rg[2*d+1] = fc_gx_rg[2*d+1];

    // THIS WILL HAVE TO CHANGE FOR NON-CONFORMING!
    sub->field_size_m *= (fc_gx_rg[2*d+1]+1-fc_gx_rg[2*d]);
    sub->field_size_p *= (fc_gx_rg[2*d+1]+1-fc_gx_rg[2*d]);
  }

  switch(face)
  {
    case 0:
      BFAM_ASSERT(sub_m->dim > 0);
      sub->vl_ix_rg[0] = sub_m->buf_sz[0];
      sub->vl_ix_rg[1] = sub->vl_ix_rg[0];
      if(sub_m->dim>1)
      {
        sub->vl_ix_rg[2] = sub_m->buf_sz[2] + fc_gx_rg[0] - sub_m->gx[1];
        sub->vl_ix_rg[3] = sub_m->buf_sz[2] + fc_gx_rg[1] - sub_m->gx[1];
      }
      if(sub_m->dim>2)
      {
        sub->vl_ix_rg[4] = sub_m->buf_sz[4] + fc_gx_rg[2] - sub_m->gx[2];
        sub->vl_ix_rg[5] = sub_m->buf_sz[4] + fc_gx_rg[3] - sub_m->gx[2];
      }
      break;

    case 1:
      BFAM_ASSERT(sub_m->dim > 0);
      sub->vl_ix_rg[0] = sub_m->buf_sz[0] + sub_m->Nl[0];
      sub->vl_ix_rg[1] = sub->vl_ix_rg[0];
      if(sub_m->dim>1)
      {
        sub->vl_ix_rg[2] = sub_m->buf_sz[2] + fc_gx_rg[0] - sub_m->gx[1];
        sub->vl_ix_rg[3] = sub_m->buf_sz[2] + fc_gx_rg[1] - sub_m->gx[1];
      }
      if(sub_m->dim>2)
      {
        sub->vl_ix_rg[4] = sub_m->buf_sz[4] + fc_gx_rg[2] - sub_m->gx[2];
        sub->vl_ix_rg[5] = sub_m->buf_sz[4] + fc_gx_rg[3] - sub_m->gx[2];
      }
      break;

    case 2:
      BFAM_ASSERT(sub_m->dim > 1);
      sub->vl_ix_rg[0] = sub_m->buf_sz[0] + fc_gx_rg[0] - sub_m->gx[0];
      sub->vl_ix_rg[1] = sub_m->buf_sz[0] + fc_gx_rg[1] - sub_m->gx[0];
      sub->vl_ix_rg[2] = sub_m->buf_sz[2];
      sub->vl_ix_rg[3] = sub->vl_ix_rg[2];
      if(sub_m->dim>2)
      {
        sub->vl_ix_rg[4] = sub_m->buf_sz[4] + fc_gx_rg[2] - sub_m->gx[2];
        sub->vl_ix_rg[5] = sub_m->buf_sz[4] + fc_gx_rg[3] - sub_m->gx[2];
      }
      break;

    case 3:
      BFAM_ASSERT(sub_m->dim > 1);
      sub->vl_ix_rg[0] = sub_m->buf_sz[0] + fc_gx_rg[0] - sub_m->gx[0];
      sub->vl_ix_rg[1] = sub_m->buf_sz[0] + fc_gx_rg[1] - sub_m->gx[0];
      sub->vl_ix_rg[2] = sub_m->buf_sz[2] + sub_m->Nl[1];
      sub->vl_ix_rg[3] = sub->vl_ix_rg[2];
      if(sub_m->dim>2)
      {
        sub->vl_ix_rg[4] = sub_m->buf_sz[4] + fc_gx_rg[2] - sub_m->gx[2];
        sub->vl_ix_rg[5] = sub_m->buf_sz[4] + fc_gx_rg[3] - sub_m->gx[2];
      }
      break;

    case 4:
      BFAM_ASSERT(sub_m->dim > 2);
      sub->vl_ix_rg[0] = sub_m->buf_sz[0] + fc_gx_rg[0] - sub_m->gx[0];
      sub->vl_ix_rg[1] = sub_m->buf_sz[0] + fc_gx_rg[1] - sub_m->gx[0];
      sub->vl_ix_rg[2] = sub_m->buf_sz[2] + fc_gx_rg[2] - sub_m->gx[1];
      sub->vl_ix_rg[3] = sub_m->buf_sz[2] + fc_gx_rg[3] - sub_m->gx[1];
      sub->vl_ix_rg[4] = sub_m->buf_sz[4];
      sub->vl_ix_rg[5] = sub->vl_ix_rg[4];
      break;

    case 5:
      BFAM_ASSERT(sub_m->dim > 2);
      sub->vl_ix_rg[0] = sub_m->buf_sz[0] + fc_gx_rg[0] - sub_m->gx[0];
      sub->vl_ix_rg[1] = sub_m->buf_sz[0] + fc_gx_rg[1] - sub_m->gx[0];
      sub->vl_ix_rg[2] = sub_m->buf_sz[2] + fc_gx_rg[2] - sub_m->gx[1];
      sub->vl_ix_rg[3] = sub_m->buf_sz[2] + fc_gx_rg[3] - sub_m->gx[1];
      sub->vl_ix_rg[4] = sub_m->buf_sz[4] + sub_m->Nl[2];
      sub->vl_ix_rg[5] = sub->vl_ix_rg[4];
      break;

    default:
      BFAM_ABORT("Not implemented for dim > 3");
  }
}

void
bfam_subdomain_sbp_inter_glue_free(bfam_subdomain_t *thisSubdomain)
{
  bfam_subdomain_sbp_inter_glue_t* sub =
    (bfam_subdomain_sbp_inter_glue_t*) thisSubdomain;

  bfam_dictionary_allprefixed_ptr(&sub->base.fields,"",
      &bfam_subdomain_sbp_free_fields,NULL);
  bfam_dictionary_allprefixed_ptr(&sub->base.fields_p,"",
      &bfam_subdomain_sbp_free_fields,NULL);
  bfam_dictionary_allprefixed_ptr(&sub->base.fields_m,"",
      &bfam_subdomain_sbp_free_fields,NULL);
  bfam_dictionary_allprefixed_ptr(&sub->base.fields_face,"",
      &bfam_subdomain_sbp_free_fields,NULL);

  bfam_subdomain_free(thisSubdomain);

  bfam_free(sub->fc_gx_rg);
  bfam_free(sub->vl_ix_rg);
}
