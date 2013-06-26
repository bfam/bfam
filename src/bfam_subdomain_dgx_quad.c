#include <bfam_subdomain_dgx_quad.h>
#include <bfam_jacobi.h>
#include <bfam_log.h>
#include <bfam_util.h>
#include <bfam_vtk.h>

static void
bfam_subdomain_dgx_quad_vtk_write_vtu_piece(bfam_subdomain_t *subdomain,
    FILE *file, const char **scalars, const char **vectors,
    const char **components, int writeBinary, int writeCompressed,
    int rank, bfam_locidx_t id)
{
  bfam_subdomain_dgx_quad_t *s = (bfam_subdomain_dgx_quad_t*) subdomain;

  const char *format;

  if(writeBinary)
    format = "binary";
  else
    format = "ascii";

  const bfam_locidx_t  K  = s->K;
  const int            N  = s->N;
  const int            Np = s->Np;

  const int Ncorners = s->Ncorners;

  const bfam_locidx_t Ncells = K * N * N;
  const bfam_locidx_t Ntotal = K * Np;

  bfam_real_t *restrict x =
    bfam_dictionary_get_value_ptr(&subdomain->fields, "_grid_x");
  bfam_real_t *restrict y =
    bfam_dictionary_get_value_ptr(&subdomain->fields, "_grid_y");
  bfam_real_t *restrict z =
    bfam_dictionary_get_value_ptr(&subdomain->fields, "_grid_z");

  fprintf(file,
           "    <Piece NumberOfPoints=\"%jd\" NumberOfCells=\"%jd\">\n",
           (intmax_t) Ntotal, (intmax_t) Ncells);

  /*
   * Points
   */
  fprintf (file, "      <Points>\n");

  bfam_vtk_write_real_vector_data_array(file, "Position", writeBinary,
      writeCompressed, Ntotal, x, y, z);

  fprintf(file, "      </Points>\n");

  /*
   * Cells
   */
  fprintf(file, "      <Cells>\n");

  /*
   * Connectivity
   */
  fprintf(file, "        <DataArray type=\"%s\" Name=\"connectivity\""
          " format=\"%s\">\n", BFAM_LOCIDX_VTK, format);
  if(writeBinary)
  {
    size_t cellsSize = Ncells*Ncorners*sizeof(bfam_locidx_t);
    bfam_locidx_t *cells = bfam_malloc_aligned(cellsSize);

    for(bfam_locidx_t k = 0, i = 0; k < K; ++k)
    {
      for(int m = 0; m < N; ++m)
      {
        for(int n = 0; n < N; ++n)
        {
          cells[i++] = Np * k + (N + 1) * (m + 0) + (n + 0);
          cells[i++] = Np * k + (N + 1) * (m + 0) + (n + 1);
          cells[i++] = Np * k + (N + 1) * (m + 1) + (n + 0);
          cells[i++] = Np * k + (N + 1) * (m + 1) + (n + 1);
        }
      }
    }

    fprintf(file, "          ");
    int rval =
      bfam_vtk_write_binary_data(writeCompressed, file, (char*)cells,
          cellsSize);
    fprintf(file, "\n");
    if(rval)
      BFAM_WARNING("Error encoding cells");

    bfam_free_aligned(cells);
  }
  else
  {
    for(bfam_locidx_t k = 0; k < K; ++k)
      for(int m = 0; m < N; ++m)
        for(int n = 0; n < N; ++n)
          fprintf(file,
                   "          %8jd %8jd %8jd %8jd\n",
                   (intmax_t) Np * k + (N + 1) * (m + 0) + (n + 0),
                   (intmax_t) Np * k + (N + 1) * (m + 0) + (n + 1),
                   (intmax_t) Np * k + (N + 1) * (m + 1) + (n + 0),
                   (intmax_t) Np * k + (N + 1) * (m + 1) + (n + 1));
  }
  fprintf(file, "        </DataArray>\n");

  /*
   * Offsets
   */
  fprintf (file, "        <DataArray type=\"%s\" Name=\"offsets\""
           " format=\"%s\">\n", BFAM_LOCIDX_VTK, format);
  fprintf(file, "          ");
  if(writeBinary)
  {
    size_t offsetsSize = Ncells*sizeof(bfam_locidx_t);
    bfam_locidx_t *offsets = bfam_malloc_aligned(offsetsSize);

    for(bfam_locidx_t i = 1; i <= Ncells; ++i)
      offsets[i - 1] = Ncorners * i;

    int rval =
      bfam_vtk_write_binary_data(writeCompressed, file, (char*)offsets,
          offsetsSize);
    if(rval)
      BFAM_WARNING("Error encoding offsets");

    bfam_free_aligned(offsets);
  }
  else
  {
    for(bfam_locidx_t i = 1, sk = 1; i <= Ncells; ++i, ++sk)
    {
      fprintf(file, " %8jd", (intmax_t) (Ncorners * i));
      if(!(sk % 20) && i != Ncells)
        fprintf(file, "\n          ");
    }
  }
  fprintf(file, "\n");
  fprintf(file, "        </DataArray>\n");

  /*
   * Types
   */
  fprintf(file, "        <DataArray type=\"UInt8\" Name=\"types\""
           " format=\"%s\">\n", format);
  fprintf(file, "          ");
  if(writeBinary)
  {
    size_t typesSize = Ncells*sizeof(uint8_t);
    uint8_t *types = bfam_malloc_aligned(typesSize);

    for(bfam_locidx_t i = 0; i < Ncells; ++i)
      types[i] = 8; /* VTK_PIXEL */

    int rval =
      bfam_vtk_write_binary_data(writeCompressed, file, (char*)types,
          typesSize);
    if(rval)
      BFAM_WARNING("Error encoding types");

    bfam_free_aligned(types);
  }
  else
  {
    for(bfam_locidx_t i = 0, sk = 1; i < Ncells; ++i, ++sk)
    {
      fprintf(file, " 8"); /* VTK_PIXEL */
      if (!(sk % 20) && i != (Ncells - 1))
        fprintf(file, "\n         ");
    }
  }
  fprintf(file, "\n");
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "      </Cells>\n");

  /*
   * Cell Data
   */
  fprintf(file, "      <CellData Scalars=\"mpirank,subdomain_id\">\n");
  fprintf(file, "        <DataArray type=\"%s\" Name=\"mpirank\""
           " format=\"%s\">\n", BFAM_LOCIDX_VTK, format);
  fprintf(file, "          ");
  if(writeBinary)
  {
    size_t ranksSize = Ncells*sizeof(bfam_locidx_t);
    bfam_locidx_t *ranks = bfam_malloc_aligned(ranksSize);

    for(bfam_locidx_t i = 0; i < Ncells; ++i)
      ranks[i] = rank;

    int rval =
      bfam_vtk_write_binary_data(writeCompressed, file, (char*)ranks,
          ranksSize);
    if(rval)
      BFAM_WARNING("Error encoding ranks");

    bfam_free_aligned(ranks);
  }
  else
  {
    for(bfam_locidx_t i = 0, sk = 1; i < Ncells; ++i, ++sk)
    {
      fprintf(file, " %6jd", (intmax_t)rank);
      if (!(sk % 8) && i != (Ncells - 1))
        fprintf(file, "\n         ");
    }
  }
  fprintf(file, "\n");
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "        <DataArray type=\"%s\" Name=\"subdomain_id\""
           " format=\"%s\">\n", BFAM_LOCIDX_VTK, format);
  fprintf(file, "          ");
  if(writeBinary)
  {
    size_t idsSize = Ncells*sizeof(bfam_locidx_t);
    bfam_locidx_t *ids = bfam_malloc_aligned(idsSize);

    for(bfam_locidx_t i = 0; i < Ncells; ++i)
      ids[i] = id;

    int rval =
      bfam_vtk_write_binary_data(writeCompressed, file, (char*)ids,
          idsSize);
    if(rval)
      BFAM_WARNING("Error encoding ids");

    bfam_free_aligned(ids);
  }
  else
  {
    for(bfam_locidx_t i = 0, sk = 1; i < Ncells; ++i, ++sk)
    {
      fprintf(file, " %6jd", (intmax_t)id);
      if (!(sk % 8) && i != (Ncells - 1))
        fprintf(file, "\n         ");
    }
  }
  fprintf(file, "\n");
  fprintf(file, "        </DataArray>\n");

  fprintf(file, "      </CellData>\n");

  char pointscalars[BFAM_BUFSIZ];
  bfam_util_strcsl(pointscalars, scalars);

  char pointvectors[BFAM_BUFSIZ];
  bfam_util_strcsl(pointvectors, vectors);

  fprintf(file, "      <PointData Scalars=\"%s\" Vectors=\"%s\">\n",
      pointscalars, pointvectors);

  if(scalars)
  {
    for(size_t s = 0; scalars[s]; ++s)
    {
      bfam_real_t *sdata = bfam_dictionary_get_value_ptr(&subdomain->fields,
          scalars[s]);
      BFAM_ABORT_IF(sdata == NULL, "VTK: Field %s not in subdomain %s",
          scalars[s], subdomain->name);
      bfam_vtk_write_real_scalar_data_array(file, scalars[s],
          writeBinary, writeCompressed, Ntotal, sdata);
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

      BFAM_ABORT_IF(v1 == NULL, "VTK: Field %s not in subdomain %s",
          components[3*v+0], subdomain->name);
      BFAM_ABORT_IF(v2 == NULL, "VTK: Field %s not in subdomain %s",
          components[3*v+1], subdomain->name);
      BFAM_ABORT_IF(v3 == NULL, "VTK: Field %s not in subdomain %s",
          components[3*v+2], subdomain->name);

      bfam_vtk_write_real_vector_data_array(file, vectors[v],
          writeBinary, writeCompressed, Ntotal, v1, v2, v3);
    }
  }

  fprintf(file, "      </PointData>\n");
  fprintf(file, "    </Piece>\n");
}


bfam_subdomain_dgx_quad_t*
bfam_subdomain_dgx_quad_new(const bfam_locidx_t     id,
                            const char             *name,
                            const int               N,
                            const bfam_locidx_t     Nv,
                            const bfam_long_real_t *VX,
                            const bfam_long_real_t *VY,
                            const bfam_long_real_t *VZ,
                            const bfam_locidx_t     K,
                            const bfam_locidx_t    *EToV,
                            const bfam_locidx_t    *EToE,
                            const int8_t           *EToF)
{
  bfam_subdomain_dgx_quad_t* newSubdomain =
    bfam_malloc(sizeof(bfam_subdomain_dgx_quad_t));

  bfam_subdomain_dgx_quad_init(newSubdomain, id, name, N, Nv, VX, VY, VZ, K,
                               EToV, EToE, EToF);

  return newSubdomain;
}

static int
bfam_subdomain_dgx_quad_field_add(bfam_subdomain_t *subdomain, const char *name)
{
  bfam_subdomain_dgx_quad_t *s = (bfam_subdomain_dgx_quad_t*) subdomain;

  if(bfam_dictionary_get_value_ptr(&s->base.fields,name))
    return 1;

  size_t fieldSize = s->Np*s->K*sizeof(bfam_real_t);
  bfam_real_t *field = bfam_malloc_aligned(fieldSize);

  int rval = bfam_dictionary_insert_ptr(&s->base.fields, name, field);

  BFAM_ASSERT(rval != 1);

  if(rval == 0)
    bfam_free_aligned(field);

  return rval;
}

static void
bfam_subdomain_dgx_quad_field_init(bfam_subdomain_t *subdomain,
    const char *name, bfam_real_t time, bfam_subdomain_init_field_t init_field,
    void *arg)
{
  bfam_subdomain_dgx_quad_t *s = (bfam_subdomain_dgx_quad_t*) subdomain;

  bfam_real_t *field = bfam_dictionary_get_value_ptr(&s->base.fields,name);

  BFAM_ABORT_IF(field==NULL, "Init: Field %s not found in subdomain %s",
      name, subdomain->name);

  size_t fieldLength = s->Np*s->K;

  bfam_real_t *restrict x =
    bfam_dictionary_get_value_ptr(&subdomain->fields, "_grid_x");
  bfam_real_t *restrict y =
    bfam_dictionary_get_value_ptr(&subdomain->fields, "_grid_y");
  bfam_real_t *restrict z =
    bfam_dictionary_get_value_ptr(&subdomain->fields, "_grid_z");

  init_field(fieldLength, time, x, y, z, subdomain, arg, field);
}

void
bfam_subdomain_dgx_quad_init(bfam_subdomain_dgx_quad_t       *subdomain,
                             const bfam_locidx_t              id,
                             const char                      *name,
                             const int                        N,
                             const bfam_locidx_t              Nv,
                             const bfam_long_real_t          *VX,
                             const bfam_long_real_t          *VY,
                             const bfam_long_real_t          *VZ,
                             const bfam_locidx_t              K,
                             const bfam_locidx_t             *EToV,
                             const bfam_locidx_t             *EToE,
                             const int8_t                    *EToF)
{
  bfam_subdomain_init(&subdomain->base, id, name);
  bfam_subdomain_add_tag(&subdomain->base, "_subdomain_dgx_quad");
  subdomain->base.free = bfam_subdomain_dgx_quad_free;
  subdomain->base.vtk_write_vtu_piece =
    bfam_subdomain_dgx_quad_vtk_write_vtu_piece;
  subdomain->base.field_add = bfam_subdomain_dgx_quad_field_add;
  subdomain->base.field_init = bfam_subdomain_dgx_quad_field_init;

  const int Np = (N+1)*(N+1);
  const int Nfp = N+1;
  const int Nfaces = 4;
  const int Ncorners = 4;
  const int Nh = 3;
  const int No = 2;

  const int Nrp = N+1;
  bfam_long_real_t *lr, *lw;
  lr = bfam_malloc_aligned(Nrp*sizeof(bfam_long_real_t));
  lw = bfam_malloc_aligned(Nrp*sizeof(bfam_long_real_t));

  bfam_jacobi_gauss_lobatto_quadrature(0, 0, N, lr, lw);

  bfam_long_real_t *lx, *ly, *lz;
  lx = bfam_malloc_aligned(K*Np*sizeof(bfam_long_real_t));
  ly = bfam_malloc_aligned(K*Np*sizeof(bfam_long_real_t));
  lz = bfam_malloc_aligned(K*Np*sizeof(bfam_long_real_t));


  for(bfam_locidx_t k = 0; k < K; ++k)
  {
    bfam_locidx_t va = EToV[Ncorners * k + 0];
    bfam_locidx_t vb = EToV[Ncorners * k + 1];
    bfam_locidx_t vc = EToV[Ncorners * k + 2];
    bfam_locidx_t vd = EToV[Ncorners * k + 3];

    for(int n = 0; n < Nrp; ++n)
    {
      for(int m = 0; m < Nrp; ++m)
      {
        bfam_long_real_t wa, wb, wc, wd;

        wa = (1-lr[m])*(1-lr[n]);
        wb = (1+lr[m])*(1-lr[n]);
        wc = (1-lr[m])*(1+lr[n]);
        wd = (1+lr[m])*(1+lr[n]);

        lx[Np*k + Nrp*n + m] = (wa*VX[va]+wb*VX[vb]+wc*VX[vc]+wd*VX[vd])/4;
        ly[Np*k + Nrp*n + m] = (wa*VY[va]+wb*VY[vb]+wc*VY[vc]+wd*VY[vd])/4;
        lz[Np*k + Nrp*n + m] = (wa*VZ[va]+wb*VZ[vb]+wc*VZ[vc]+wd*VZ[vd])/4;
      }
    }
  }

  bfam_long_real_t *restrict V =
    bfam_malloc_aligned(Nrp*Nrp*sizeof(bfam_long_real_t));

  bfam_jacobi_p_vandermonde(0, 0, N, Nrp, lr, V);

  bfam_long_real_t *restrict D =
    bfam_malloc_aligned(Nrp*Nrp*sizeof(bfam_long_real_t));

  bfam_jacobi_p_differentiation(0, 0, N, Nrp, lr, V, D);

  /*
   * Set subdomain values
   */
  subdomain->N = N;
  subdomain->Np = Np;
  subdomain->Nfp = Nfp;
  subdomain->Nfaces = Nfaces;
  subdomain->Ncorners = Ncorners;
  subdomain->Nh = Nh;
  subdomain->No = No;

  subdomain->r = bfam_malloc_aligned(Nrp*sizeof(bfam_real_t));
  subdomain->w = bfam_malloc_aligned(Nrp*sizeof(bfam_real_t));
  for(int n = 0; n<Nrp; ++n)
  {
    subdomain->r[n] = (bfam_real_t) lr[n];
    subdomain->w[n] = (bfam_real_t) lw[n];
  }

  subdomain->K = K;

  int rval;
  rval = bfam_subdomain_dgx_quad_field_add(&subdomain->base, "_grid_x");
  BFAM_ASSERT(rval == 2);
  rval = bfam_subdomain_dgx_quad_field_add(&subdomain->base, "_grid_y");
  BFAM_ASSERT(rval == 2);
  rval = bfam_subdomain_dgx_quad_field_add(&subdomain->base, "_grid_z");
  BFAM_ASSERT(rval == 2);

  bfam_real_t *restrict x =
    bfam_dictionary_get_value_ptr(&subdomain->base.fields, "_grid_x");
  bfam_real_t *restrict y =
    bfam_dictionary_get_value_ptr(&subdomain->base.fields, "_grid_y");
  bfam_real_t *restrict z =
    bfam_dictionary_get_value_ptr(&subdomain->base.fields, "_grid_z");

  for(int n = 0; n < K*Np; ++n)
  {
    x[n] = (bfam_real_t) lx[n];
    y[n] = (bfam_real_t) ly[n];
    z[n] = (bfam_real_t) lz[n];
  }

  subdomain->Dr = bfam_malloc_aligned(Nrp * Nrp * sizeof(bfam_real_t));
  for(int n = 0; n < Nrp*Nrp; ++n)
    subdomain->Dr[n] = (bfam_real_t) D[n];

  subdomain->fmask = bfam_malloc_aligned(Nfaces * sizeof(int*));

  for(int f = 0; f < Nfaces; ++f)
    subdomain->fmask[f] = bfam_malloc_aligned(Nfp * sizeof(int));

  /*
   * Face 0 -x
   */
  for(int i = 0; i < N+1; ++i)
    subdomain->fmask[0][i] = i*(N+1);

  /*
   * Face 1 +x
   */
  for(int i = 0; i < N+1; ++i)
    subdomain->fmask[1][i] = (i+1)*(N+1)-1;

  /*
   * Face 2 -y
   */
  for(int i = 0; i < N+1; ++i)
    subdomain->fmask[2][i] = i;

  /*
   * Face 3 +y
   */
  for(int i = 0; i < N+1; ++i)
    subdomain->fmask[3][i] = (N+1)*N + i;

  bfam_free_aligned(D);
  bfam_free_aligned(V);

  bfam_free_aligned(lr);
  bfam_free_aligned(lw);

  bfam_free_aligned(lx);
  bfam_free_aligned(ly);
  bfam_free_aligned(lz);
}

static int
bfam_subdomain_dgx_quad_free_fields(const char * key, void *val,
    void *arg)
{
  bfam_free_aligned(val);

  return 1;
}

void
bfam_subdomain_dgx_quad_free(bfam_subdomain_t *thisSubdomain)
{
  bfam_subdomain_dgx_quad_t *sub = (bfam_subdomain_dgx_quad_t*) thisSubdomain;

  bfam_dictionary_allprefixed_ptr(&sub->base.fields,"",
      &bfam_subdomain_dgx_quad_free_fields,NULL);
  bfam_dictionary_allprefixed_ptr(&sub->base.fields_p,"",
      &bfam_subdomain_dgx_quad_free_fields,NULL);
  bfam_dictionary_allprefixed_ptr(&sub->base.fields_m,"",
      &bfam_subdomain_dgx_quad_free_fields,NULL);


  bfam_subdomain_free(thisSubdomain);

  bfam_free_aligned(sub->Dr);

  bfam_free_aligned(sub->r);
  bfam_free_aligned(sub->w);

  for(int f = 0; f < sub->Nfaces; ++f)
    bfam_free_aligned(sub->fmask[f]);
  bfam_free_aligned(sub->fmask);
}

static int
bfam_subdomain_dgx_quad_glue_field_minus_add(bfam_subdomain_t *subdomain,
    const char *name)
{
  bfam_subdomain_dgx_quad_glue_t *s =
    (bfam_subdomain_dgx_quad_glue_t*) subdomain;

  if(bfam_dictionary_get_value_ptr(&s->base.fields_m,name))
    return 1;

  size_t fieldSize = s->Np*s->K*sizeof(bfam_real_t);
  bfam_real_t *field = bfam_malloc_aligned(fieldSize);

  int rval = bfam_dictionary_insert_ptr(&s->base.fields_m, name, field);

  BFAM_ASSERT(rval != 1);

  if(rval == 0)
    bfam_free_aligned(field);

  return rval;
}

static int
bfam_subdomain_dgx_quad_glue_field_plus_add(bfam_subdomain_t *subdomain,
    const char *name)
{
  bfam_subdomain_dgx_quad_glue_t *s =
    (bfam_subdomain_dgx_quad_glue_t*) subdomain;

  if(bfam_dictionary_get_value_ptr(&s->base.fields_p,name))
    return 1;

  size_t fieldSize = s->Np*s->K*sizeof(bfam_real_t);
  bfam_real_t *field = bfam_malloc_aligned(fieldSize);

  int rval = bfam_dictionary_insert_ptr(&s->base.fields_p, name, field);

  BFAM_ASSERT(rval != 1);

  if(rval == 0)
    bfam_free_aligned(field);

  return rval;
}

static void
bfam_subdomain_dgx_quad_glue_comm_info(bfam_subdomain_t *thisSubdomain,
    int *rank, bfam_locidx_t *my_id, bfam_locidx_t *neigh_id,
    size_t *send_sz, size_t *recv_sz)
{
  bfam_subdomain_dgx_quad_glue_t *sub =
    (bfam_subdomain_dgx_quad_glue_t*) thisSubdomain;

  *rank     = sub->rank_p;
  *neigh_id = sub->id_p;
  *my_id    = sub->id_m;

  const size_t send_num = sub->base.fields_m.num_entries * sub->K * sub->Np;
  const size_t recv_num = sub->base.fields_p.num_entries * sub->K * sub->Np;

  *send_sz  = send_num*sizeof(bfam_real_t);
  *recv_sz  = recv_num*sizeof(bfam_real_t);

  BFAM_LDEBUG(" rank %3d   ms %3jd   ns %3jd   send_sz %3zd   recv_sz %3zd",
      *rank, (intmax_t) *my_id, (intmax_t) *neigh_id, *send_sz, *recv_sz);
}

typedef struct bfam_subdomain_dgx_quad_get_put_data
{
  bfam_subdomain_dgx_quad_glue_t *sub;
  bfam_real_t *buffer;
  size_t size;
  size_t field;
} bfam_subdomain_dgx_quad_get_put_data_t;

static int
bfam_subdomain_dgx_quad_glue_get_fields_m(const char * key, void *val,
    void *arg)
{
  bfam_subdomain_dgx_quad_get_put_data_t *data =
    (bfam_subdomain_dgx_quad_get_put_data_t*) arg;

  bfam_subdomain_dgx_quad_glue_t *sub = data->sub;
  const bfam_locidx_t K = sub->K;
  const int Np = sub->Np;

  const bfam_locidx_t *restrict EToEp = sub->EToEp;
  const bfam_locidx_t *restrict EToEm = sub->EToEm;
  const int8_t        *restrict EToFm = sub->EToFm;
  const int8_t        *restrict EToHm = sub->EToHm;
  const int8_t        *restrict EToOm = sub->EToOm;

  const int sub_m_Np = sub->sub_m->Np;
  const int sub_m_Nfp = sub->sub_m->Nfp;

  const size_t buffer_offset = data->field * Np * K;

  bfam_real_t *restrict send_field = data->buffer + buffer_offset;

  BFAM_ASSERT((data->field+1) * Np * K * sizeof(bfam_real_t) <= data->size);

  const bfam_real_t *restrict sub_m_field =
    bfam_dictionary_get_value_ptr(&sub->sub_m->base.fields, key);

  bfam_real_t *restrict glue_field = val;

  BFAM_ASSERT( send_field != NULL);
  BFAM_ASSERT(sub_m_field != NULL);
  BFAM_ASSERT( glue_field != NULL);

  BFAM_ASSUME_ALIGNED(sub_m_field, 32);
  BFAM_ASSUME_ALIGNED( glue_field, 32);

  BFAM_ASSUME_ALIGNED(EToEp, 32);
  BFAM_ASSUME_ALIGNED(EToEm, 32);
  BFAM_ASSUME_ALIGNED(EToFm, 32);
  BFAM_ASSUME_ALIGNED(EToHm, 32);
  BFAM_ASSUME_ALIGNED(EToOm, 32);

  for(bfam_locidx_t k = 0; k < K; ++k)
  {
    BFAM_ASSERT(EToEp[k] < sub->K);
    BFAM_ASSERT(EToEm[k] < sub->sub_m->K);
    BFAM_ASSERT(EToFm[k] < sub->sub_m->Nfaces);
    BFAM_ASSERT(EToHm[k] < sub->sub_m->Nh);
    BFAM_ASSERT(EToOm[k] < sub->sub_m->No);

    bfam_real_t *restrict send_elem = send_field + EToEp[k] * Np;

    bfam_locidx_t *restrict fmask = sub->sub_m->fmask[EToFm[k]];

    const bfam_real_t *restrict sub_m_elem = sub_m_field + EToEm[k] * sub_m_Np;

    bfam_real_t *restrict glue_elem = glue_field + k * Np;

    /*
     * Decide which interpolation operation to use.
     */
    const bfam_real_t *restrict interpolation = sub->interpolation[EToHm[k]];
    BFAM_ASSUME_ALIGNED(interpolation, 32);

    /*
     * Interpolate.
     */
    if(interpolation)
    {
      /*
       * XXX: Replace with something faster; this will also have to change
       * for 3D.
       */
      for(int n = 0; n < Np; ++n)
        glue_elem[n] = 0;
      for(int j = 0; j < sub_m_Nfp; ++j)
        for(int i = 0; i < Np; ++i)
          glue_elem[i] += interpolation[j * Np + i] * sub_m_elem[fmask[j]];
    }
    else
    {
      for(int i = 0; i < Np; ++i)
        glue_elem[i] = sub_m_elem[fmask[i]];
    }

    /*
     * Copy data to send buffer based on orientation.
     */
    if(EToOm[k])
    {
      for(int n = 0; n < Np; ++n)
        send_elem[n] = glue_elem[Np-1-n];
    }
    else
    {
      memcpy(send_elem, glue_elem, Np * sizeof(bfam_real_t));
    }
  }

  ++data->field;
  return 0;
}

void
bfam_subdomain_dgx_quad_glue_put_send_buffer(bfam_subdomain_t *thisSubdomain,
    void *buffer, size_t send_sz)
{
  bfam_subdomain_dgx_quad_get_put_data_t data;

  data.sub    = (bfam_subdomain_dgx_quad_glue_t*) thisSubdomain;
  data.buffer = (bfam_real_t*) buffer;
  data.size   = send_sz;
  data.field  = 0;

  BFAM_ASSERT(send_sz == data.sub->base.fields_m.num_entries * data.sub->K *
      data.sub->Np * sizeof(bfam_real_t));

  /*
   * Fill fields_m and the send buffer from sub_m.
   */
  bfam_dictionary_allprefixed_ptr(&data.sub->base.fields_m, "",
      &bfam_subdomain_dgx_quad_glue_get_fields_m, &data);
}

static int
bfam_subdomain_dgx_quad_glue_put_fields_p(const char * key, void *val,
    void *arg)
{
  bfam_subdomain_dgx_quad_get_put_data_t *data =
    (bfam_subdomain_dgx_quad_get_put_data_t*) arg;

  const bfam_locidx_t K = data->sub->K;
  const int Np = data->sub->Np;

  const size_t buffer_offset = data->field * Np * K;
  const bfam_real_t *restrict recv_field = data->buffer + buffer_offset;

  BFAM_ASSERT((data->field+1) * Np * K * sizeof(bfam_real_t) <= data->size);

  bfam_real_t *restrict glue_field = val;

  BFAM_ASSUME_ALIGNED(glue_field, 32);

  memcpy(glue_field, recv_field, K * Np * sizeof(bfam_real_t));

  ++data->field;
  return 0;
}

void
bfam_subdomain_dgx_quad_glue_get_recv_buffer(bfam_subdomain_t *thisSubdomain,
    void *buffer, size_t recv_sz)
{
  bfam_subdomain_dgx_quad_get_put_data_t data;

  data.sub    = (bfam_subdomain_dgx_quad_glue_t*) thisSubdomain;
  data.buffer = (bfam_real_t*) buffer;
  data.size   = recv_sz;
  data.field  = 0;

  BFAM_ASSERT(recv_sz == data.sub->base.fields_m.num_entries * data.sub->K *
      data.sub->Np * sizeof(bfam_real_t));

  /*
   * Fill fields_p from the recv buffer.
   */
  bfam_dictionary_allprefixed_ptr(&data.sub->base.fields_p, "",
      &bfam_subdomain_dgx_quad_glue_put_fields_p, &data);

  BFAM_ASSERT(recv_sz == data.sub->base.fields_p.num_entries * data.sub->K *
      data.sub->Np * sizeof(bfam_real_t));
}


bfam_subdomain_dgx_quad_glue_t*
bfam_subdomain_dgx_quad_glue_new(const bfam_locidx_t              id,
                                 const char                      *name,
                                 const int                        N_m,
                                 const int                        N_p,
                                 const bfam_locidx_t              rank_m,
                                 const bfam_locidx_t              rank_p,
                                 const bfam_locidx_t              id_m,
                                 const bfam_locidx_t              id_p,
                                 bfam_subdomain_dgx_quad_t       *sub_m,
                                 bfam_locidx_t                   *ktok_m,
                                 const bfam_locidx_t              K,
                                 bfam_subdomain_face_map_entry_t *mapping)
{
  bfam_subdomain_dgx_quad_glue_t* newSubdomain =
    bfam_malloc(sizeof(bfam_subdomain_dgx_quad_glue_t));

  bfam_subdomain_dgx_quad_glue_init(newSubdomain, id, name, N_m, N_p, rank_m,
      rank_p, id_m, id_p, sub_m, ktok_m, K, mapping);

  return newSubdomain;
}

void
bfam_subdomain_dgx_quad_glue_init(bfam_subdomain_dgx_quad_glue_t  *subdomain,
                                  const bfam_locidx_t              id,
                                  const char                      *name,
                                  const int                        N_m,
                                  const int                        N_p,
                                  const bfam_locidx_t              rank_m,
                                  const bfam_locidx_t              rank_p,
                                  const bfam_locidx_t              id_m,
                                  const bfam_locidx_t              id_p,
                                  bfam_subdomain_dgx_quad_t       *sub_m,
                                  bfam_locidx_t                   *ktok_m,
                                  const bfam_locidx_t              K,
                                  bfam_subdomain_face_map_entry_t *mapping)
{
  bfam_subdomain_init(&subdomain->base, id, name);
  bfam_subdomain_add_tag(&subdomain->base, "_subdomain_dgx_quad");
  subdomain->base.glue_comm_info = bfam_subdomain_dgx_quad_glue_comm_info;
  subdomain->base.glue_put_send_buffer =
    bfam_subdomain_dgx_quad_glue_put_send_buffer;
  subdomain->base.glue_get_recv_buffer =
    bfam_subdomain_dgx_quad_glue_get_recv_buffer;
  subdomain->base.field_minus_add =
    bfam_subdomain_dgx_quad_glue_field_minus_add;
  subdomain->base.field_plus_add = bfam_subdomain_dgx_quad_glue_field_plus_add;
  subdomain->base.free = bfam_subdomain_dgx_quad_glue_free;

  const int N   = BFAM_MAX(N_m, N_p);
  const int Nrp = N+1;
  const int sub_m_Nrp = N_m+1;

  const int num_interp = 3;

  bfam_long_real_t **lr;
  lr = bfam_malloc_aligned(num_interp*sizeof(bfam_long_real_t*));

  for(int i = 0; i < num_interp; ++i)
    lr[i] = bfam_malloc_aligned(Nrp*sizeof(bfam_long_real_t));

  bfam_long_real_t *restrict lw;
  lw = bfam_malloc_aligned(Nrp*sizeof(bfam_long_real_t));

  bfam_jacobi_gauss_lobatto_quadrature(0, 0, N, lr[0], lw);

  const bfam_long_real_t half = BFAM_LONG_REAL(0.5);
  for(int n = 0; n < Nrp; ++n)
  {
    lr[1][n] = -half + half * lr[0][n];
    lr[2][n] =  half + half * lr[0][n];
  }

  bfam_long_real_t *restrict sub_m_lr, *restrict sub_m_lw;
  sub_m_lr = bfam_malloc_aligned(sub_m_Nrp*sizeof(bfam_long_real_t));
  sub_m_lw = bfam_malloc_aligned(sub_m_Nrp*sizeof(bfam_long_real_t));

  bfam_jacobi_gauss_lobatto_quadrature(0, 0, N_m, sub_m_lr, sub_m_lw);

  bfam_long_real_t *restrict sub_m_V;
  sub_m_V = bfam_malloc_aligned(sub_m_Nrp*sub_m_Nrp*sizeof(bfam_long_real_t));

  bfam_jacobi_p_vandermonde(0, 0, N_m, N_m+1, sub_m_lr, sub_m_V);

  bfam_long_real_t **interpolation =
    bfam_malloc_aligned(num_interp*sizeof(bfam_long_real_t*));

  for(int i = 0; i < num_interp; ++i)
  {
    interpolation[i] =
      bfam_malloc_aligned(Nrp*sub_m_Nrp*sizeof(bfam_long_real_t));

    bfam_jacobi_p_interpolation(0, 0, N_m, Nrp, lr[i], sub_m_V,
        interpolation[i]);
  }

  subdomain->N_m = N_m;
  subdomain->N_p = N_p;

  subdomain->rank_m = rank_m;
  subdomain->rank_p = rank_p;

  subdomain->id_m = id_m;
  subdomain->id_p = id_p;

  subdomain->s_m = imaxabs(id_m)-1;
  subdomain->s_p = imaxabs(id_p)-1;

  subdomain->sub_m = sub_m;

  subdomain->N        = N;
  subdomain->Np       = Nrp;
  subdomain->Nfp      = 1;
  subdomain->Nfaces   = 2;
  subdomain->Ncorners = 2;

  subdomain->r = bfam_malloc_aligned(Nrp*sizeof(bfam_real_t));
  subdomain->w = bfam_malloc_aligned(Nrp*sizeof(bfam_real_t));
  for(int n = 0; n < Nrp; ++n)
  {
    subdomain->r[n] = (bfam_real_t) lr[0][n];
    subdomain->w[n] = (bfam_real_t) lw[n];
  }

  subdomain->K = K;

  subdomain->num_interp = num_interp;

  subdomain->interpolation =
    bfam_malloc_aligned(subdomain->num_interp * sizeof(bfam_real_t*));

  for(int i = 0; i < subdomain->num_interp; ++i)
  {
    if(i == 0 && N_m == N)
    {
      /*
       * Identity interpolation operator
       */
      subdomain->interpolation[i] = NULL;
    }
    else
    {
      subdomain->interpolation[i] =
        bfam_malloc_aligned(Nrp * sub_m_Nrp * sizeof(bfam_real_t));
      for(int n = 0; n < Nrp*sub_m_Nrp; ++n)
        subdomain->interpolation[i][n] = (bfam_real_t) interpolation[i][n];
    }
  }

  subdomain->EToEp = bfam_malloc_aligned(K*sizeof(bfam_locidx_t));
  subdomain->EToEm = bfam_malloc_aligned(K*sizeof(bfam_locidx_t));
  subdomain->EToFm = bfam_malloc_aligned(K*sizeof(int8_t));
  subdomain->EToHm = bfam_malloc_aligned(K*sizeof(int8_t));
  subdomain->EToOm = bfam_malloc_aligned(K*sizeof(int8_t));

  qsort(mapping, K, sizeof(bfam_subdomain_face_map_entry_t),
      bfam_subdomain_face_send_cmp);

  for(bfam_locidx_t k = 0; k < K; ++k)
    mapping[k].i = k;

  qsort(mapping, K, sizeof(bfam_subdomain_face_map_entry_t),
      bfam_subdomain_face_recv_cmp);

  for(bfam_locidx_t k = 0; k < K; ++k)
  {
    subdomain->EToEm[k] = ktok_m[mapping[k].k];
    subdomain->EToFm[k] = mapping[k].f;
    subdomain->EToHm[k] = mapping[k].h;
    subdomain->EToOm[k] = mapping[k].o;

    subdomain->EToEp[k] = mapping[k].i;
  }

#ifdef BFAM_DEBUG
  for(bfam_locidx_t k = 0; k < K; ++k)
    BFAM_ASSERT(mapping[k].s  == subdomain->s_m &&
                mapping[k].ns == subdomain->s_p);
#endif

  for(int i = 0; i < num_interp; ++i)
  {
    bfam_free_aligned(lr[i]);
    bfam_free_aligned(interpolation[i]);
  }

  bfam_free_aligned(lr);
  bfam_free_aligned(lw);

  bfam_free_aligned(interpolation);

  bfam_free_aligned(sub_m_lr);
  bfam_free_aligned(sub_m_lw);
  bfam_free_aligned(sub_m_V);
}

void
bfam_subdomain_dgx_quad_glue_free(bfam_subdomain_t *subdomain)
{
  bfam_subdomain_dgx_quad_glue_t *sub =
    (bfam_subdomain_dgx_quad_glue_t*) subdomain;

  bfam_dictionary_allprefixed_ptr(&sub->base.fields,"",
      &bfam_subdomain_dgx_quad_free_fields,NULL);
  bfam_dictionary_allprefixed_ptr(&sub->base.fields_p,"",
      &bfam_subdomain_dgx_quad_free_fields,NULL);
  bfam_dictionary_allprefixed_ptr(&sub->base.fields_m,"",
      &bfam_subdomain_dgx_quad_free_fields,NULL);

  bfam_free_aligned(sub->r);
  bfam_free_aligned(sub->w);

  for(int i = 0; i < sub->num_interp; ++i)
    if(sub->interpolation[i])
      bfam_free_aligned(sub->interpolation[i]);
  bfam_free_aligned(sub->interpolation);

  bfam_free_aligned(sub->EToEp);
  bfam_free_aligned(sub->EToEm);
  bfam_free_aligned(sub->EToFm);
  bfam_free_aligned(sub->EToHm);
  bfam_free_aligned(sub->EToOm);

  bfam_subdomain_free(subdomain);
}
