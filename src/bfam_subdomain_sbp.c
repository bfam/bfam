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
                            const bfam_locidx_t* sub_N,
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
  bfam_subdomain_sbp_init(newSub,id,sub_ix,sub_N,name,dim,N,Nl,Nb,gx,
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
  bfam_locidx_t *storage_locid = bfam_malloc_aligned(Ntot*sizeof(bfam_locidx_t));

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
      storage_locid[i] = s->loc_id;

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
      fprintf(file, " %6jd", (intmax_t)s->loc_id);
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
              &sdata[i+Nb[0] + (Nt[0]+1)*(j+Nb[2] + (k+Nb[4])*(Nt[1]+1))],
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
                &v1[i+Nb[0] + (Nt[0]+1)*(j+Nb[2] + (k+Nb[4])*(Nt[1]+1))],
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
                &v2[i+Nb[0] + (Nt[0]+1)*(j+Nb[2] + (k+Nb[4])*(Nt[1]+1))],
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
                &v3[i+Nb[0] + (Nt[0]+1)*(j+Nb[2] + (k+Nb[4])*(Nt[1]+1))],
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
    fieldSize *= (s->Nl[d]+1+s->Nb[2*d]+s->Nb[2*d+1]);

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
    fieldSize *= (s->Nl[d]+1+s->Nb[2*d]+s->Nb[2*d+1]);

  bfam_real_t *restrict x =
    bfam_dictionary_get_value_ptr(&subdomain->fields, "_grid_x");
  bfam_real_t *restrict y =
    bfam_dictionary_get_value_ptr(&subdomain->fields, "_grid_y");
  bfam_real_t *restrict z =
    bfam_dictionary_get_value_ptr(&subdomain->fields, "_grid_z");

  init_field(fieldSize, time, x, y, z, subdomain, arg, field);
}

void
bfam_subdomain_sbp_init(bfam_subdomain_sbp_t *subdomain,
                            const bfam_locidx_t     id,
                            const bfam_locidx_t* sub_ix,
                            const bfam_locidx_t* sub_N,
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

  subdomain->sub_ix = bfam_malloc(dim*sizeof(bfam_locidx_t));
  subdomain->sub_N  = bfam_malloc(dim*sizeof(bfam_locidx_t));
  subdomain->loc_id = 0;
  subdomain->num_id = 1;
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

    subdomain->sub_ix[d] = sub_ix[d];
    subdomain->sub_N [d] = sub_N[d];
    subdomain->loc_id += sub_ix[d]*subdomain->num_id;
    subdomain->num_id *= (sub_N[d]+1);

    BFAM_ABORT_IF(N [d] < 0 || Nl[d] < 0 || Nb[2*d] < 0 || Nb[2*d+1] < 0 ||
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
  bfam_free(sub->sub_ix);
  bfam_free(sub->sub_N);
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
  sort[2] = sub->loc_m;
  sort[3] = sub->loc_p;
  sort[4] = sub->fce_m;
  sort[5] = sub->fce_p;

  size_t sr_num = 1;
  for(int d = 0;d < sub_m->dim;d++)
    sr_num *= (sub->ix[2*d+1]+1-sub->ix[2*d]);

#ifdef BFAM_DEBUG
  if(sub_m->dim == 2)
  {
    BFAM_LDEBUG(" r_m %3d r_p %3d Nl[%3d %3d] Nb[%3d %3d %3d %3d]"
        " ix[%3d %3d %3d %3d] size %3d",
        sub->rank_m,sub->rank_p,
        sub_m->Nl[0],sub_m->Nl[1],
        sub_m->Nb[0],sub_m->Nb[1],sub_m->Nb[2],sub_m->Nb[3],
        sub->ix[0],sub->ix[1],sub->ix[2],sub->ix[3],
        sr_num);
  }
  if(sub_m->dim == 3)
  {
    BFAM_LDEBUG(" r_m %3d r_p %3d Nl[%3d %3d %3d] Nb[%3d %3d %3d %3d %3d %3d]"
        " ix[%3d %3d %3d %3d %3d %3d] size %3d",
        sub->rank_m,sub->rank_p,
        sub_m->Nl[0],sub_m->Nl[1],sub_m->Nl[2],
        sub_m->Nb[0],sub_m->Nb[1],sub_m->Nb[2],
        sub_m->Nb[3],sub_m->Nb[4],sub_m->Nb[5],
        sub->ix[0],sub->ix[1],sub->ix[2],
        sub->ix[3],sub->ix[4],sub->ix[5],
        sr_num);
  }
#endif


  *send_sz  = sub->base.fields_p.num_entries*sr_num*sizeof(bfam_real_t);
  *recv_sz  = sub->base.fields_m.num_entries*sr_num*sizeof(bfam_real_t);
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

  bfam_subdomain_sbp_t *sub_m = data->sub->sub_m;

  BFAM_ASSERT(sub_m->dim < 4);

  // FINISH

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
  // sub->base.glue_put_send_buffer =
  //   bfam_subdomain_sbp_intra_glue_put_send_buffer;
  // sub->base.glue_get_recv_buffer =
  //   bfam_subdomain_sbp_intra_glue_get_recv_buffer;
  sub->base.field_minus_add = bfam_subdomain_sbp_intra_glue_field_minus_add;
  sub->base.field_plus_add  = bfam_subdomain_sbp_intra_glue_field_plus_add;
  sub->base.free = bfam_subdomain_sbp_intra_glue_free;

  sub->rank_m = rank_m;
  sub->rank_p = rank_p;
  sub->sub_m  = sub_m;
  sub->face   = face;

  sub->id_m  = sub_m->base.id;
  sub->loc_m = sub_m->loc_id;
  sub->fce_m = face;

  /* initialize all the indices to full size */
  sub->ix = bfam_malloc(2*sub_m->dim*sizeof(bfam_locidx_t));
  for(int d = 0;d < sub_m->dim;d++)
  {
    sub->ix[2*d  ] = sub_m->Nb[2*d];
    sub->ix[2*d+1] = sub_m->Nb[2*d]+sub_m->Nl[d];
  }

  sub->id_p  = sub_m->base.id; /* assume we number in the same order! */
  switch(face)
  {
    case 0 :
      sub->loc_p  = sub_m->loc_id-1;
      sub->ix[0] = 0;
      sub->ix[1] = sub_m->Nb[0]-1;
      break;
    case 1 :
      sub->loc_p  = sub_m->loc_id+1;
      sub->ix[0] = sub_m->Nb[0] + sub_m->Nl[0]+1;
      sub->ix[1] = sub_m->Nb[0] + sub_m->Nl[0]+1 + sub_m->Nb[1]-1;
      break;
    case 2 :
      sub->loc_p  = sub_m->loc_id-(sub_m->sub_N[0]+1);
      sub->ix[2] = 0;
      sub->ix[3] = sub_m->Nb[2]-1;
      break;
    case 3 :
      sub->loc_p  = sub_m->loc_id+(sub_m->sub_N[0]+1);
      sub->ix[2] = sub_m->Nb[2] + sub_m->Nl[1]+1;
      sub->ix[3] = sub_m->Nb[2] + sub_m->Nl[1]+1 + sub_m->Nb[3]-1;
      break;
    case 4 :
      sub->loc_p  = sub_m->loc_id-(sub_m->sub_N[0]+1)*(sub_m->sub_N[1]+1);
      sub->ix[4] = 0;
      sub->ix[5] = sub_m->Nb[4]-1;
      break;
    case 5 :
      sub->loc_p  = sub_m->loc_id+(sub_m->sub_N[0]+1)*(sub_m->sub_N[1]+1);
      sub->ix[4] = sub_m->Nb[4] + sub_m->Nl[2]+1;
      sub->ix[5] = sub_m->Nb[4] + sub_m->Nl[2]+1 + sub_m->Nb[5]-1;
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
  bfam_subdomain_free(thisSubdomain);
  bfam_free(sub->ix);
}

