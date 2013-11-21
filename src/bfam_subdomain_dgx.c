#include <bfam_subdomain_dgx.h>
#include <bfam_jacobi.h>
#include <bfam_kron.h>
#include <bfam_log.h>
#include <bfam_util.h>
#include <bfam_vtk.h>

#ifndef BFAM_DGX_DIMENSION
#define BFAM_DGX_DIMENSION
#define USE_GENERIC_DGX_DIMENSION
#else
#define DIM (BFAM_DGX_DIMENSION)
#endif

static void
bfam_subdomain_dgx_vtk_interp(bfam_locidx_t K,
    int N_d,       bfam_real_t * restrict d,
    int N_s, const bfam_real_t * restrict s,
    const bfam_real_t *restrict interp)
{
  BFAM_ABORT("interp not implemented");
}

static int
bfam_subdomain_dgx_vtk_write_vtu_piece(bfam_subdomain_t *subdomain,
    FILE *file, bfam_real_t time, const char **scalars, const char **vectors,
    const char **components, int writeBinary, int writeCompressed,
    int rank, bfam_locidx_t id, int Np_write)
{
  bfam_subdomain_dgx_t *sub = (bfam_subdomain_dgx_t*) subdomain;

  const char *format;

  if(writeBinary)
    format = "binary";
  else
    format = "ascii";

  const bfam_locidx_t  K  = sub->K;
  int              N_vtk  = sub->N;
  int              Np_vtk = sub->Np;
  bfam_real_t *interp = NULL;

  bfam_real_t *restrict stor1 = NULL;
  bfam_real_t *restrict stor2 = NULL;
  bfam_real_t *restrict stor3 = NULL;

  BFAM_INFO("Np_write = %d", Np_write);

  if(Np_write > 0)
  {
    BFAM_ABORT_IF_NOT(Np_write > 1, "Np_write = %d is not valid",
        Np_write);

    N_vtk  = Np_write - 1;
    Np_vtk = Np_write*Np_write;

    interp = bfam_malloc_aligned(sizeof(bfam_real_t)*(sub->N+1)*(N_vtk+1));

    bfam_long_real_t *calc_interp =
      bfam_malloc_aligned(sizeof(bfam_long_real_t)*(sub->N+1)*(N_vtk+1));
    bfam_long_real_t *lr =
      bfam_malloc_aligned(sizeof(bfam_long_real_t)*Np_write);

    for(int r = 0; r < Np_write; r++)
      lr[r] = -1 + 2*(bfam_long_real_t)r/(Np_write-1);

    bfam_jacobi_p_interpolation(0, 0, sub->N, Np_write, lr, sub->V,
        calc_interp);

    for(int n = 0; n < (sub->N+1)*(N_vtk+1); n++)
      interp[n] = (bfam_real_t)calc_interp[n];

    stor1 = bfam_malloc_aligned(sizeof(bfam_real_t)*Np_vtk*K);
    stor2 = bfam_malloc_aligned(sizeof(bfam_real_t)*Np_vtk*K);
    stor3 = bfam_malloc_aligned(sizeof(bfam_real_t)*Np_vtk*K);

    bfam_free_aligned(lr);
    bfam_free_aligned(calc_interp);
  }

  const int Ncorners = sub->Ng[sub->numg-1];

  const bfam_locidx_t Ncells = K * N_vtk * N_vtk;
  const bfam_locidx_t Ntotal = K * Np_vtk;

  bfam_real_t *restrict x =
    bfam_dictionary_get_value_ptr(&subdomain->fields, "_grid_x0");
  bfam_real_t *restrict y =
    bfam_dictionary_get_value_ptr(&subdomain->fields, "_grid_x1");
  bfam_real_t *restrict z =
    bfam_dictionary_get_value_ptr(&subdomain->fields, "_grid_x2");

  if(interp == NULL)
  {
    stor1 = x;
    stor2 = y;
    stor3 = z;
  }
  else
  {
    bfam_subdomain_dgx_vtk_interp(K,N_vtk,stor1,sub->N,x,interp);
    bfam_subdomain_dgx_vtk_interp(K,N_vtk,stor2,sub->N,y,interp);
    bfam_subdomain_dgx_vtk_interp(K,N_vtk,stor3,sub->N,z,interp);
  }

  fprintf(file,
           "    <Piece NumberOfPoints=\"%jd\" NumberOfCells=\"%jd\">\n",
           (intmax_t) Ntotal, (intmax_t) Ncells);

  /*
   * Points
   */
  fprintf (file, "      <Points>\n");

  bfam_vtk_write_real_vector_data_array(file, "Position", writeBinary,
      writeCompressed, Ntotal, stor1, stor2, stor3);

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

    if(sub->dim == 1)
      for(bfam_locidx_t k = 0, i = 0; k < K; ++k)
      {
        for(int n = 0; n < N_vtk; ++n)
        {
          cells[i++] = Np_vtk * k + (n + 0);
          cells[i++] = Np_vtk * k + (n + 1);
        }
      }
    else if(sub->dim == 2)
      for(bfam_locidx_t k = 0, i = 0; k < K; ++k)
      {
        for(int m = 0; m < N_vtk; ++m)
        {
          for(int n = 0; n < N_vtk; ++n)
          {
            cells[i++] = Np_vtk * k + (N_vtk + 1) * (m + 0) + (n + 0);
            cells[i++] = Np_vtk * k + (N_vtk + 1) * (m + 0) + (n + 1);
            cells[i++] = Np_vtk * k + (N_vtk + 1) * (m + 1) + (n + 0);
            cells[i++] = Np_vtk * k + (N_vtk + 1) * (m + 1) + (n + 1);
          }
        }
      }
    else if(sub->dim == 3)
      for(bfam_locidx_t k = 0, i = 0; k < K; ++k)
      {
        for(int l = 0; l < N_vtk; ++l)
        {
          for(int m = 0; m < N_vtk; ++m)
          {
            for(int n = 0; n < N_vtk; ++n)
            {
              cells[i++] = Np_vtk*k +
                (N_vtk+1)*(N_vtk+1)*(l+0) + (N_vtk+1)*(m+0) + (n+0);
              cells[i++] = Np_vtk*k +
                (N_vtk+1)*(N_vtk+1)*(l+0) + (N_vtk+1)*(m+0) + (n+1);
              cells[i++] = Np_vtk*k +
                (N_vtk+1)*(N_vtk+1)*(l+0) + (N_vtk+1)*(m+1) + (n+0);
              cells[i++] = Np_vtk*k +
                (N_vtk+1)*(N_vtk+1)*(l+0) + (N_vtk+1)*(m+1) + (n+1);
              cells[i++] = Np_vtk*k +
                (N_vtk+1)*(N_vtk+1)*(l+1) + (N_vtk+1)*(m+0) + (n+0);
              cells[i++] = Np_vtk*k +
                (N_vtk+1)*(N_vtk+1)*(l+1) + (N_vtk+1)*(m+0) + (n+1);
              cells[i++] = Np_vtk*k +
                (N_vtk+1)*(N_vtk+1)*(l+1) + (N_vtk+1)*(m+1) + (n+0);
              cells[i++] = Np_vtk*k +
                (N_vtk+1)*(N_vtk+1)*(l+1) + (N_vtk+1)*(m+1) + (n+1);
            }
          }
        }
      }
    else BFAM_ABORT("not implemented for dim = %d", sub->dim);

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
    if(sub->dim == 1)
      for(bfam_locidx_t k = 0; k < K; ++k)
        for(int n = 0; n < N_vtk; ++n)
          fprintf(file,
              "          %8jd %8jd\n",
              (intmax_t) Np_vtk * k + (n + 0),
              (intmax_t) Np_vtk * k + (n + 1));
    if(sub->dim == 2)
      for(bfam_locidx_t k = 0; k < K; ++k)
        for(int m = 0; m < N_vtk; ++m)
          for(int n = 0; n < N_vtk; ++n)
            fprintf(file,
                "          %8jd %8jd %8jd %8jd\n",
                (intmax_t) Np_vtk * k + (N_vtk + 1) * (m + 0) + (n + 0),
                (intmax_t) Np_vtk * k + (N_vtk + 1) * (m + 0) + (n + 1),
                (intmax_t) Np_vtk * k + (N_vtk + 1) * (m + 1) + (n + 0),
                (intmax_t) Np_vtk * k + (N_vtk + 1) * (m + 1) + (n + 1));
    if(sub->dim == 3)
      for(bfam_locidx_t k = 0; k < K; ++k)
        for(int l = 0; l < N_vtk; ++l)
          for(int m = 0; m < N_vtk; ++m)
            for(int n = 0; n < N_vtk; ++n)
              fprintf(file,
                  "          %8jd %8jd %8jd %8jd %8jd %8jd %8jd %8jd\n",
                  (intmax_t) Np_vtk*k
                      +(N_vtk+1)*(N_vtk+1)*(l+0) + (N_vtk+1)*(m+0) + (n+0),
                  (intmax_t) Np_vtk*k
                      +(N_vtk+1)*(N_vtk+1)*(l+0) + (N_vtk+1)*(m+0) + (n+1),
                  (intmax_t) Np_vtk*k
                      +(N_vtk+1)*(N_vtk+1)*(l+0) + (N_vtk+1)*(m+1) + (n+0),
                  (intmax_t) Np_vtk*k
                      +(N_vtk+1)*(N_vtk+1)*(l+0) + (N_vtk+1)*(m+1) + (n+1),
                  (intmax_t) Np_vtk*k
                      +(N_vtk+1)*(N_vtk+1)*(l+1) + (N_vtk+1)*(m+0) + (n+0),
                  (intmax_t) Np_vtk*k
                      +(N_vtk+1)*(N_vtk+1)*(l+1) + (N_vtk+1)*(m+0) + (n+1),
                  (intmax_t) Np_vtk*k
                      +(N_vtk+1)*(N_vtk+1)*(l+1) + (N_vtk+1)*(m+1) + (n+0),
                  (intmax_t) Np_vtk*k
                      +(N_vtk+1)*(N_vtk+1)*(l+1) + (N_vtk+1)*(m+1) + (n+1));
    else BFAM_ABORT("not implemented for dim = %d", sub->dim);
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
  fprintf(file, "      <CellData Scalars=\"time,mpirank,subdomain_id\">\n");
  fprintf(file, "        <DataArray type=\"%s\" Name=\"time\""
           " format=\"%s\">\n", BFAM_REAL_VTK, format);
  fprintf(file, "          ");
  if(writeBinary)
  {
    size_t timesize = Ncells*sizeof(bfam_real_t);
    bfam_real_t *times = bfam_malloc_aligned(timesize);

    for(bfam_locidx_t i = 0; i < Ncells; ++i)
      times[i] = time;

    int rval =
      bfam_vtk_write_binary_data(writeCompressed, file, (char*)times,
          timesize);
    if(rval)
      BFAM_WARNING("Error encoding times");

    bfam_free_aligned(times);
  }
  else
  {
    for(bfam_locidx_t i = 0, sk = 1; i < Ncells; ++i, ++sk)
    {
      fprintf(file, " %"BFAM_REAL_FMTe, time);
      if (!(sk % 8) && i != (Ncells - 1))
        fprintf(file, "\n         ");
    }
  }
  fprintf(file, "\n");
  fprintf(file, "        </DataArray>\n");
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
      if(interp == NULL)
      {
        stor1 = sdata;
      }
      else
      {
        bfam_subdomain_dgx_vtk_interp(K,N_vtk,stor1,sub->N,sdata,interp);
      }

      bfam_vtk_write_real_scalar_data_array(file, scalars[s],
          writeBinary, writeCompressed, Ntotal, stor1);
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
      if(interp == NULL)
      {
        stor1 = v1;
        stor2 = v2;
        stor3 = v3;
      }
      else
      {
        bfam_subdomain_dgx_vtk_interp(K,N_vtk,stor1,sub->N,v1,interp);
        bfam_subdomain_dgx_vtk_interp(K,N_vtk,stor2,sub->N,v2,interp);
        bfam_subdomain_dgx_vtk_interp(K,N_vtk,stor3,sub->N,v3,interp);
      }

      bfam_vtk_write_real_vector_data_array(file, vectors[v],
          writeBinary, writeCompressed, Ntotal, stor1, stor2, stor3);
    }
  }

  fprintf(file, "      </PointData>\n");
  fprintf(file, "    </Piece>\n");

  if(interp != NULL)
  {
    bfam_free_aligned(interp);
    bfam_free_aligned(stor1);
    bfam_free_aligned(stor2);
    bfam_free_aligned(stor3);
  }
  return 1;
}

static int
bfam_subdomain_dgx_field_add(bfam_subdomain_t *subdomain, const char *name)
{
  bfam_subdomain_dgx_t *s = (bfam_subdomain_dgx_t*) subdomain;

  if(bfam_dictionary_get_value_ptr(&s->base.fields,name))
    return 1;

  size_t fieldSize = s->Np*s->K*sizeof(bfam_real_t);
  bfam_real_t *field = bfam_malloc_aligned(fieldSize);
#ifdef BFAM_DEBUG
  for(int i = 0; i < s->Np*s->K;i++) field[i] = bfam_real_nan("");
#endif

  int rval = bfam_dictionary_insert_ptr(&s->base.fields, name, field);

  BFAM_ASSERT(rval != 1);

  if(rval == 0)
    bfam_free_aligned(field);

  return rval;
}

static inline int***
BFAM_APPEND_EXPAND(bfam_subdomain_dgx_gmask_set_,BFAM_DGX_DIMENSION)
         (const int numg, const int N, int *Np, int *Ng, int *Ngp, int inDIM)
{
#ifdef USE_GENERIC_DGX_DIMENSION
  BFAM_WARNING("Using generic bfam_subdomain_dgx_gmask_set");
  const int DIM = inDIM;
#endif
  BFAM_ABORT_IF(DIM > 3 || DIM < 0,
      "bfam_subdomain_dgx_gmask_set cannot handle dim = %d",DIM);

  if(DIM == 0)
  {
    *Np = 1;
    return NULL;
  }


  /* this could probably be made generic for arbitrary dimensions, but until
   * that's needed... */
  switch(DIM)
  {
    case 1:
      *Np = N+1;

      /* just corners */
      Ng [0]  = 2;
      Ngp[0]  = 1;
      break;

    case 2:
      *Np = (N+1)*(N+1);

      /* edges */
      Ng [0]  = 4;
      Ngp[0]  = N+1;

      /* corners */
      Ng [1]  = 4;
      Ngp[1]  = 1;
      break;

    case 3:
      *Np = (N+1)*(N+1)*(N+1);

      /* faces */
      Ng [0]  = 6;
      Ngp[0]  = (N+1)*(N+1);

      /* edges */
      Ng [1]  = 12;
      Ngp[1]  = N+1;

      /* corners */
      Ng [2]  = 8;
      Ngp[2]  = 1;
      break;

    default:
      BFAM_ABORT("cannot handle dim = %d",DIM);
  }

  int ***gmask = bfam_malloc_aligned(numg * sizeof(int**));
  for(int g = 0; g < numg; g++)
  {
    gmask[g] = bfam_malloc_aligned(Ng[g] * sizeof(int*));
    for(int i = 0; i < Ng[g]; i++)
      gmask[g][i] = bfam_malloc_aligned(Ngp[g] * sizeof(int));
  }

  switch(DIM)
  {
    case 1:
      gmask[0][0][0] = 0;
      gmask[0][1][0] = N;
      break;

    case 2:
      /* edges */
      for(int i = 0; i < N+1; ++i) gmask[0][0][i] = i*(N+1);
      for(int i = 0; i < N+1; ++i) gmask[0][1][i] = (i+1)*(N+1)-1;
      for(int i = 0; i < N+1; ++i) gmask[0][2][i] = i;
      for(int i = 0; i < N+1; ++i) gmask[0][3][i] = (N+1)*N + i;

      /* corners */
      for(int j = 0; j < 2; ++j)
        for(int i = 0; i < 2; ++i)
          gmask[1][i+j*2][0] = i*N + j*(N+1);
      break;

    case 3:
      /* This could all probably be cleaned up... */

      /* faces */
      {
        int n,i,j,k,f=-1;

        n = 0; i = 0; f++;
        for(k = 0; k < N+1; k++) for(j = 0; j < N+1; j++)
          gmask[0][f][n++] = i+j*(N+1)+k*(N+1)*(N+1);

        n = 0; i = N; f++;
        for(k = 0; k < N+1; k++) for(j = 0; j < N+1; j++)
          gmask[0][f][n++] = i+j*(N+1)+k*(N+1)*(N+1);

        n = 0; j = 0; f++;
        for(k = 0; k < N+1; k++) for(i = 0; i < N+1; i++)
          gmask[0][f][n++] = i+j*(N+1)+k*(N+1)*(N+1);

        n = 0; j = N; f++;
        for(k = 0; k < N+1; k++) for(i = 0; i < N+1; i++)
          gmask[0][f][n++] = i+j*(N+1)+k*(N+1)*(N+1);

        n = 0; k = 0; f++;
        for(j = 0; j < N+1; j++) for(i = 0; i < N+1; i++)
          gmask[0][f][n++] = i+j*(N+1)+k*(N+1)*(N+1);

        n = 0; k = N; f++;
        for(j = 0; j < N+1; j++) for(i = 0; i < N+1; i++)
          gmask[0][f][n++] = i+j*(N+1)+k*(N+1)*(N+1);
      }

      /* edges */
      {
        int n,i,j,k,e = 0;

        for(k = 0; k < N+1;k+=N)
          for(j = 0; j < N+1;j+=N)
          {
            n = 0;
            for(i = 0; i < N+1; i++)
              gmask[1][e][n++] = i+j*(N+1)+k*(N+1)*(N+1);
            e++;
          }
        for(k = 0; k < N+1;k+=N)
          for(i = 0; i < N+1;i+=N)
          {
            n = 0;
            for(j = 0; j < N+1; j++)
              gmask[1][e][n++] = i+j*(N+1)+k*(N+1)*(N+1);
            e++;
          }
        for(j = 0; j < N+1;j+=N)
          for(i = 0; i < N+1;i+=N)
          {
            n = 0;
            for(k = 0; k < N+1; k++)
              gmask[1][e][n++] = i+j*(N+1)+k*(N+1)*(N+1);
            e++;
          }
      }

      /* corners */
      for(int k = 0, c = 0; k < N+1; k+=N)
        for(int j = 0; j < N+1;j+=N)
          for(int i = 0; i < N+1;i+=N)
            gmask[2][c++][0] = i+j*(N+1)+k*(N+1)*(N+1);

      break;

    default:
      BFAM_ABORT("cannot handle dim = %d",DIM);
  }

  return gmask;
}


void
BFAM_APPEND_EXPAND(bfam_subdomain_dgx_init_,BFAM_DGX_DIMENSION)(
                              bfam_subdomain_dgx_t *subdomain,
                        const bfam_locidx_t         id,
                        const char                 *name,
                        const int                   N,
                        const bfam_locidx_t         Nv,
                        const int                   num_Vi,
                        const bfam_long_real_t    **Vi,
                        const bfam_locidx_t         K,
                        const bfam_locidx_t        *EToV,
                        const bfam_locidx_t        *EToE,
                        const int8_t               *EToF,
                        const int                   inDIM)
{
#ifdef USE_GENERIC_DGX_DIMENSION
  BFAM_WARNING("Using generic bfam_subdomain_dgx_init");
  const int DIM = inDIM;
#endif
  BFAM_ASSERT(DIM == inDIM);
  BFAM_ABORT_IF(DIM < 0, "dimension %d is not possible in bfam",DIM);
  BFAM_ABORT_IF(DIM == 0 && N != 0,
                "if DIM < 1 then N must be zero (i.e., constant");

  bfam_subdomain_init(&subdomain->base, id, name);
  bfam_subdomain_add_tag(&subdomain->base, "_subdomain_dgx");
  char dim_str[BFAM_BUFSIZ];
  snprintf(dim_str,BFAM_BUFSIZ,"_dimension_%d",DIM);
  bfam_subdomain_add_tag(&subdomain->base, dim_str);
  subdomain->dim = DIM;

  subdomain->base.free =
              BFAM_APPEND_EXPAND(bfam_subdomain_dgx_free_,BFAM_DGX_DIMENSION);
  subdomain->base.vtk_write_vtu_piece =
    bfam_subdomain_dgx_vtk_write_vtu_piece;
  subdomain->base.field_add = bfam_subdomain_dgx_field_add;
  // subdomain->base.field_face_add = bfam_subdomain_dgx_quad_field_face_add;
  // subdomain->base.field_init = bfam_subdomain_dgx_quad_field_init;

  subdomain->numg = DIM;
  const int numg = subdomain->numg;

  int *Ng  = NULL;
  if(numg > 0) Ng  = bfam_malloc_aligned(sizeof(int)*numg);
  subdomain->Ng = Ng;

  int *Ngp  = NULL;
  if(numg > 0) Ngp = bfam_malloc_aligned(sizeof(int)*numg);
  subdomain->Ngp = Ngp;

  subdomain->gmask =
    BFAM_APPEND_EXPAND(bfam_subdomain_dgx_gmask_set_,BFAM_DGX_DIMENSION)(
                       numg, N, &subdomain->Np, Ng, Ngp, DIM);

  const int Np = subdomain->Np;

  if(DIM > 0)
  {
    const int Nrp = N+1;
    bfam_long_real_t *lr, *lw;
    lr = bfam_malloc_aligned(Nrp*sizeof(bfam_long_real_t));
    lw = bfam_malloc_aligned(Nrp*sizeof(bfam_long_real_t));

    bfam_jacobi_gauss_lobatto_quadrature(0, 0, N, lr, lw);

    bfam_long_real_t **lxi =
      bfam_malloc_aligned(num_Vi*sizeof(bfam_long_real_t*));

    for(int i = 0;i < num_Vi;i++)
      lxi[i] = bfam_malloc_aligned(K*Np*sizeof(bfam_long_real_t));

    /* Loop over all the elements and set up the grid*/
    int Ncorners = Ng[numg-1];
    for(bfam_locidx_t k = 0; k < K; ++k)
    {
      const bfam_locidx_t *v = EToV+Ncorners*k;
      bfam_long_real_t w[Ncorners];

      if(DIM == 1)
        for(int n = 0; n < Nrp; ++n)
        {
          int offset = n;
          w[0] = 1-lr[n];
          w[1] = 1+lr[n];
          for(int i = 0; i < num_Vi;i++)
          {
            lxi[i][Np*k + offset] = 0;
            for(int c = 0; c < Ncorners; c++)
              lxi[i][Np*k + offset] += w[c]*Vi[i][v[c]];
            lxi[i][Np*k + offset] /= Ncorners;
          }
        }
      else if(DIM == 2)
        for(int n = 0; n < Nrp; ++n)
          for(int m = 0; m < Nrp; ++m)
          {
            int offset = n*Nrp+m;
            w[0] = (1-lr[m])*(1-lr[n]);
            w[1] = (1+lr[m])*(1-lr[n]);
            w[2] = (1-lr[m])*(1+lr[n]);
            w[3] = (1+lr[m])*(1+lr[n]);
            for(int i = 0; i < num_Vi;i++)
            {
              lxi[i][Np*k + offset] = 0;
              for(int c = 0; c < Ncorners; c++)
                lxi[i][Np*k + offset] += w[c]*Vi[i][v[c]];
              lxi[i][Np*k + offset] /= Ncorners;
            }
        }
      else if(DIM == 3)
        for(int n = 0; n < Nrp; ++n)
          for(int m = 0; m < Nrp; ++m)
            for(int l = 0; l < Nrp; ++l)
            {
              int offset = n*Nrp*Nrp+l*Nrp+m;
              w[0] = (1-lr[m])*(1-lr[n])*(1-lr[l]);
              w[1] = (1+lr[m])*(1-lr[n])*(1-lr[l]);
              w[2] = (1-lr[m])*(1+lr[n])*(1-lr[l]);
              w[3] = (1+lr[m])*(1+lr[n])*(1-lr[l]);
              w[4] = (1-lr[m])*(1-lr[n])*(1+lr[l]);
              w[5] = (1+lr[m])*(1-lr[n])*(1+lr[l]);
              w[6] = (1-lr[m])*(1+lr[n])*(1+lr[l]);
              w[7] = (1+lr[m])*(1+lr[n])*(1+lr[l]);
              for(int i = 0; i < num_Vi;i++)
              {
                lxi[i][Np*k + offset] = 0;
                for(int c = 0; c < Ncorners; c++)
                  lxi[i][Np*k + offset] += w[c]*Vi[i][v[c]];
                lxi[i][Np*k + offset] /= Ncorners;
              }
            }
      else BFAM_ABORT("not setup of dim = %d",DIM);

    }

    subdomain->V = bfam_malloc_aligned(Nrp*Nrp*sizeof(bfam_long_real_t));
    bfam_long_real_t *restrict V = subdomain->V;

    bfam_jacobi_p_vandermonde(0, 0, N, Nrp, lr, V);

    bfam_long_real_t *restrict D =
      bfam_malloc_aligned(Nrp*Nrp*sizeof(bfam_long_real_t));

    bfam_jacobi_p_differentiation(0, 0, N, Nrp, lr, V, D);

    bfam_long_real_t *restrict M =
      bfam_malloc_aligned(Nrp*Nrp*sizeof(bfam_long_real_t));

    bfam_jacobi_p_mass(0, 0, N, V, M);


    bfam_long_real_t **lJrx =
      bfam_malloc_aligned(num_Vi*DIM*sizeof(bfam_long_real_t*));
    for(int n = 0; n < num_Vi*DIM; n++)
      lJrx[n] = bfam_malloc_aligned(K*Np*sizeof(bfam_long_real_t));

    bfam_long_real_t *lJ = bfam_malloc_aligned(K*Np*sizeof(bfam_long_real_t));

    /* Ng[0] = number of faces, Ngp[0] = number of face points */
    bfam_long_real_t **ln =
      bfam_malloc_aligned(num_Vi*Ng[0]*sizeof(bfam_long_real_t*));
    for(int n = 0; n < num_Vi*Ng[0]; n++)
      ln[n] = bfam_malloc_aligned(K*Ng[0]*Ngp[0]*sizeof(bfam_long_real_t));

    bfam_long_real_t *lsJ =
      bfam_malloc_aligned(K*Ng[0]*Ngp[0]*sizeof(bfam_long_real_t));

    /*
    bfam_subdomain_dgx_geo(N, K, subdomain->gmask, lxi, D, lJrx, lJ, ln, lsJ,
        DIM);
    */

    subdomain->r = bfam_malloc_aligned(Nrp*sizeof(bfam_real_t));
    subdomain->w = bfam_malloc_aligned(Nrp*sizeof(bfam_real_t));
    subdomain->wi = bfam_malloc_aligned(Nrp*sizeof(bfam_real_t));

    for(int n = 0; n<Nrp; ++n)
    {
      subdomain->r[n]  = (bfam_real_t) lr[n];
      subdomain->w[n]  = (bfam_real_t) lw[n];
      subdomain->wi[n] = (bfam_real_t) (1.0l/lw[n]);
    }

    subdomain->K = K;

    /* store the volume stuff */
    subdomain->N = N;
    /* store the grid */
    for(int i = 0; i < num_Vi; i++)
    {
      char name[BFAM_BUFSIZ];
      snprintf(name,BFAM_BUFSIZ,"_grid_x%d",i);
      int rval = bfam_subdomain_dgx_field_add(&subdomain->base, name);
      BFAM_ABORT_IF_NOT(rval == 2, "Error adding %s",name);
      bfam_real_t *restrict xi =
        bfam_dictionary_get_value_ptr(&subdomain->base.fields, name);
      for(int n = 0; n < K*Np; ++n)
      {
        xi[n] = (bfam_real_t) lxi[i][n];
      }
    }

    /* store the metric stuff */

    /* store the face stuff */

    subdomain->Dr = bfam_malloc_aligned(Nrp * Nrp * sizeof(bfam_real_t));
    for(int n = 0; n < Nrp*Nrp; ++n) subdomain->Dr[n] = (bfam_real_t) D[n];

    //JK subdomain->vmapP = bfam_malloc_aligned(K*Nfp*Nfaces*sizeof(bfam_locidx_t));
    //JK subdomain->vmapM = bfam_malloc_aligned(K*Nfp*Nfaces*sizeof(bfam_locidx_t));

    //JK bfam_subdomain_dgx_quad_buildmaps(K, Np, Nfp, Nfaces, EToE, EToF,
    //JK     subdomain->fmask, subdomain->vmapP, subdomain->vmapM);

    /* free stuff */
    bfam_free_aligned(lsJ);
    for(int n = 0; n < num_Vi*Ng[0]; n++)
      bfam_free_aligned(ln[n]);
    bfam_free_aligned(ln);

    bfam_free_aligned(lJ);

    for(int n = 0; n < num_Vi*DIM; n++)
      bfam_free_aligned(lJrx[n]);
    bfam_free_aligned(lJrx);

    bfam_free_aligned(D);
    bfam_free_aligned(M);

    bfam_free_aligned(lr);
    bfam_free_aligned(lw);

    for(int i = 0;i < num_Vi;i++) bfam_free_aligned(lxi[i]);
    bfam_free_aligned(lxi);
  }
  else
  {
    subdomain->Dr = NULL;
    subdomain->V  = NULL;
    subdomain->r  = NULL;
    subdomain->w  = NULL;
    subdomain->wi = NULL;
  }
}

bfam_subdomain_dgx_t*
BFAM_APPEND_EXPAND(bfam_subdomain_dgx_new_,BFAM_DGX_DIMENSION)(
                       const bfam_locidx_t      id,
                       const char              *name,
                       const int                N,
                       const bfam_locidx_t      Nv,
                       const int                num_Vi,
                       const bfam_long_real_t **Vi,
                       const bfam_locidx_t      K,
                       const bfam_locidx_t     *EToV,
                       const bfam_locidx_t     *EToE,
                       const int8_t            *EToF,
                       const int                inDIM)
{
#ifdef USE_GENERIC_DGX_DIMENSION
  BFAM_WARNING("Using generic bfam_subdomain_dgx_new");
  const int DIM = inDIM;
#endif
  BFAM_ASSERT(DIM == inDIM);

  bfam_subdomain_dgx_t* newSubdomain =
    bfam_malloc(sizeof(bfam_subdomain_dgx_t));

  BFAM_APPEND_EXPAND(bfam_subdomain_dgx_init_,BFAM_DGX_DIMENSION)(
                         newSubdomain, id, name, N, Nv, num_Vi, Vi, K, EToV,
                         EToE, EToF,DIM);
  return newSubdomain;
}

static int
bfam_subdomain_dgx_free_fields(const char * key, void *val,
    void *arg)
{
  bfam_free_aligned(val);

  return 1;
}

void
BFAM_APPEND_EXPAND(bfam_subdomain_dgx_free_,BFAM_DGX_DIMENSION)(
    bfam_subdomain_t *thisSubdomain)
{
  bfam_subdomain_dgx_t *sub = (bfam_subdomain_dgx_t*) thisSubdomain;

  bfam_dictionary_allprefixed_ptr(&sub->base.fields,"",
      &bfam_subdomain_dgx_free_fields,NULL);
  bfam_dictionary_allprefixed_ptr(&sub->base.fields_p,"",
      &bfam_subdomain_dgx_free_fields,NULL);
  bfam_dictionary_allprefixed_ptr(&sub->base.fields_m,"",
      &bfam_subdomain_dgx_free_fields,NULL);
  bfam_dictionary_allprefixed_ptr(&sub->base.fields_face,"",
      &bfam_subdomain_dgx_free_fields,NULL);

  bfam_subdomain_free(thisSubdomain);

  if(sub->Dr) bfam_free_aligned(sub->Dr); sub->Dr = NULL;
  if(sub->V ) bfam_free_aligned(sub->V ); sub->V  = NULL;

  if(sub->r ) bfam_free_aligned(sub->r ); sub->r  = NULL;
  if(sub->w ) bfam_free_aligned(sub->w ); sub->w  = NULL;
  if(sub->wi) bfam_free_aligned(sub->wi); sub->wi = NULL;

  // for(int f = 0; f < sub->Nfaces; ++f)
  //    bfam_free_aligned(sub->fmask[f]);
  // bfam_free_aligned(sub->fmask);

  if(sub->gmask)
  {
    for(int g = 0; g < sub->numg; g++)
    {
      for(int i = 0; i < sub->Ng[g]; i++)
        bfam_free_aligned(sub->gmask[g][i]);
      bfam_free_aligned(sub->gmask[g]);
    }
    bfam_free_aligned(sub->gmask);
    sub->gmask = NULL;
  }
  if(sub->Ng ) bfam_free_aligned(sub->Ng ); sub->Ng  = NULL;
  if(sub->Ngp) bfam_free_aligned(sub->Ngp); sub->Ngp = NULL;

  // bfam_free_aligned(sub->vmapP);
  // bfam_free_aligned(sub->vmapM);
}
