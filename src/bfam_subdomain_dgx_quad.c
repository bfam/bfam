#include <bfam_subdomain_dgx_quad.h>
#include <bfam_jacobi.h>
#include <bfam_log.h>

static int
bfam_vtk_write_binary_data(int compressed, FILE *file, char *data, size_t size)
{
  if(compressed)
    return sc_vtk_write_compressed(file, data, size);
  else
    return sc_vtk_write_binary(file, data, size);
}

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

  fprintf(file,
           "    <Piece NumberOfPoints=\"%jd\" NumberOfCells=\"%jd\">\n",
           (intmax_t) Ntotal, (intmax_t) Ncells);

  /*
   * Points
   */
  fprintf (file, "      <Points>\n");

  fprintf (file, "        <DataArray type=\"%s\" Name=\"Position\""
           " NumberOfComponents=\"3\" format=\"%s\">\n",
           BFAM_REAL_VTK, format);
  if(writeBinary)
  {
    size_t xyzSize = 3*Ntotal*sizeof(bfam_real_t);
    bfam_real_t *xyz = bfam_malloc_aligned(xyzSize);

    for(bfam_locidx_t n = 0; n < Ntotal; ++n)
    {
      xyz[3*n + 0] = s->x[n];
      xyz[3*n + 1] = s->y[n];
      xyz[3*n + 2] = s->z[n];
    }

    fprintf(file, "          ");
    int rval =
      bfam_vtk_write_binary_data(writeCompressed, file, (char*)xyz, xyzSize);
    fprintf(file, "\n");
    if(rval)
      BFAM_WARNING("Error encoding xyz");

    bfam_free_aligned(xyz);
  }
  else
  {
    for(bfam_locidx_t n = 0; n < Ntotal; ++n)
    {
      fprintf(file, "         %"BFAM_REAL_FMTe
                            " %"BFAM_REAL_FMTe
                            " %"BFAM_REAL_FMTe
                            "\n",
                            s->x[n], s->y[n], s->z[n]);
    }
  }

  fprintf(file, "        </DataArray>\n");
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

  fprintf(file, "    </Piece>\n");
}


bfam_subdomain_dgx_quad_t*
bfam_subdomain_dgx_quad_new(const char             *name,
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

  bfam_subdomain_dgx_quad_init(newSubdomain, name, N, Nv, VX, VY, VZ, K, EToV,
                               EToE, EToF);

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

  init_field(fieldLength, time, s->x, s->y, s->z, subdomain, arg, field);
}

void
bfam_subdomain_dgx_quad_init(bfam_subdomain_dgx_quad_t       *subdomain,
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
  bfam_subdomain_init(&subdomain->base, name);
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

  /*
   * Set subdomain values
   */
  subdomain->N = N;
  subdomain->Np = Np;
  subdomain->Nfp = Nfp;
  subdomain->Nfaces = Nfaces;
  subdomain->Ncorners = Ncorners;

  subdomain->r = bfam_malloc_aligned(Nrp*sizeof(bfam_real_t));
  subdomain->w = bfam_malloc_aligned(Nrp*sizeof(bfam_real_t));
  for(int n = 0; n<Nrp; ++n)
  {
    subdomain->r[n] = (bfam_real_t) lr[n];
    subdomain->w[n] = (bfam_real_t) lw[n];
  }

  subdomain->K = K;

  subdomain->x = bfam_malloc_aligned(K*Np*sizeof(bfam_real_t));
  subdomain->y = bfam_malloc_aligned(K*Np*sizeof(bfam_real_t));
  subdomain->z = bfam_malloc_aligned(K*Np*sizeof(bfam_real_t));
  for(int n = 0; n < K*Np; ++n)
  {
    subdomain->x[n] = (bfam_real_t) lx[n];
    subdomain->y[n] = (bfam_real_t) ly[n];
    subdomain->z[n] = (bfam_real_t) lz[n];
  }

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

  bfam_subdomain_free(thisSubdomain);

  bfam_free_aligned(sub->r);
  bfam_free_aligned(sub->w);

  bfam_free_aligned(sub->x);
  bfam_free_aligned(sub->y);
  bfam_free_aligned(sub->z);
}
