#include <bfam_base.h>
#include <bfam_vtk.h>
#include <bfam_log.h>
#include <bfam_util.h>
#include <sc.h>
#include <sc_io.h>

#define BFAM_VTK_VTU_FORMAT "%s_%05d.vtu"
#define BFAM_VTK_VTS_FORMAT "%s_%s_%s.vts"
#define BFAM_VTK_PVTS_FORMAT "%s_%s.pvts"

void bfam_vtk_write_vtu_empty(FILE *file, int writeBinary)
{
  const char *format;

  if (writeBinary)
    format = "binary";
  else
    format = "ascii";

  fprintf(file, "    <Piece NumberOfPoints=\"%jd\" NumberOfCells=\"%jd\">\n",
          (intmax_t)0, (intmax_t)0);

  /*
   * Points
   */
  fprintf(file, "      <Points>\n");
  fprintf(file, "        <DataArray type=\"%s\" Name=\"%s\""
                " NumberOfComponents=\"3\" format=\"%s\">\n",
          BFAM_REAL_VTK, "Position", format);
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
                " format=\"%s\">\n",
          BFAM_LOCIDX_VTK, format);
  fprintf(file, "        </DataArray>\n");

  /*
   * Offsets
   */
  fprintf(file, "        <DataArray type=\"%s\" Name=\"offsets\""
                " format=\"%s\">\n",
          BFAM_LOCIDX_VTK, format);
  fprintf(file, "        </DataArray>\n");

  /*
   * Types
   */
  fprintf(file, "        <DataArray type=\"UInt8\" Name=\"types\""
                " format=\"%s\">\n",
          format);
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "      </Cells>\n");
  fprintf(file, "    </Piece>\n");
}

static void bfam_vtk_write_file_pvtu(int size, const char *directory,
                                     const char *prefix, const char **scalars,
                                     const char **vectors,
                                     const char **components, int binary,
                                     int compress)
{
  const int endian = bfam_endian();

  char filename[BFAM_BUFSIZ];
  snprintf(filename, BFAM_BUFSIZ, "%s.pvtu", prefix);

  const char *format;
  if (binary)
    format = "binary";
  else
    format = "ascii";

  BFAM_VERBOSE("Writing file: '%s'", filename);
  FILE *file = fopen(filename, "w");

  if (file == NULL)
  {
    BFAM_LERROR("Could not open %s for output!\n", filename);
    return;
  }

  fprintf(file, "<?xml version=\"1.0\"?>\n");
  fprintf(file, "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\"");

  if (binary && compress)
    fprintf(file, " compressor=\"vtkZLibDataCompressor\"");

  if (endian == BFAM_BIG_ENDIAN)
    fprintf(file, " byte_order=\"BigEndian\">\n");
  else
    fprintf(file, " byte_order=\"LittleEndian\">\n");

  fprintf(file, "  <PUnstructuredGrid GhostLevel=\"0\">\n");
  fprintf(file, "    <PPoints>\n");
  fprintf(file, "      <PDataArray type=\"%s\" Name=\"Position\""
                " NumberOfComponents=\"3\" format=\"%s\"/>\n",
          BFAM_REAL_VTK, format);
  fprintf(file, "    </PPoints>\n");

  fprintf(file, "    <PCellData Scalars=\"time,mpirank,subdomain_id\">\n");
  fprintf(file, "      <PDataArray type=\"%s\" Name=\"time\""
                " format=\"%s\"/>\n",
          BFAM_REAL_VTK, format);
  fprintf(file, "      <PDataArray type=\"%s\" Name=\"mpirank\""
                " format=\"%s\"/>\n",
          BFAM_LOCIDX_VTK, format);
  fprintf(file, "      <PDataArray type=\"%s\" Name=\"subdomain_id\""
                " format=\"%s\"/>\n",
          BFAM_LOCIDX_VTK, format);
  fprintf(file, "    </PCellData>\n");

  char pointscalars[BFAM_BUFSIZ];
  bfam_util_strcsl(pointscalars, scalars);

  char pointvectors[BFAM_BUFSIZ];
  bfam_util_strcsl(pointvectors, vectors);

  fprintf(file, "    <PPointData Scalars=\"%s\" Vectors=\"%s\">\n",
          pointscalars, pointvectors);

  if (scalars)
    for (size_t s = 0; scalars[s]; ++s)
      fprintf(file, "      <PDataArray type=\"%s\" Name=\"%s\""
                    " format=\"%s\"/>\n",
              BFAM_REAL_VTK, scalars[s], format);

  if (vectors)
    for (size_t v = 0; vectors[v]; ++v)
      fprintf(file, "      <PDataArray type=\"%s\" Name=\"%s\""
                    " NumberOfComponents=\"3\" format=\"%s\"/>\n",
              BFAM_REAL_VTK, vectors[v], format);

  fprintf(file, "    </PPointData>\n");

  for (int s = 0; s < size; ++s)
  {
    if (directory)
      fprintf(file, "    <Piece Source=\"%s/" BFAM_VTK_VTU_FORMAT "\"/>\n",
              directory, prefix, s);
    else
      fprintf(file, "    <Piece Source=\"" BFAM_VTK_VTU_FORMAT "\"/>\n", prefix,
              s);
  }
  fprintf(file, "  </PUnstructuredGrid>\n");
  fprintf(file, "</VTKFile>\n");

  if (ferror(file))
  {
    BFAM_LERROR("Error writing to %s\n", filename);
  }

  if (fclose(file))
  {
    BFAM_LERROR("Error closing %s\n", filename);
  }
}

void bfam_vtk_write_file(bfam_domain_t *domain, bfam_domain_match_t match,
                         const char **tags, const char *directory,
                         const char *prefix, bfam_real_t time,
                         const char **scalars, const char **vectors,
                         const char **components, int binary, int compress,
                         int Np_write)
{
  const int endian = bfam_endian();

  bfam_locidx_t numElements = domain->numSubdomains;

  int rank, size;
  BFAM_MPI_CHECK(MPI_Comm_rank(domain->comm, &rank));
  BFAM_MPI_CHECK(MPI_Comm_size(domain->comm, &size));

  if (rank == 0)
    bfam_vtk_write_file_pvtu(size, directory, prefix, scalars, vectors,
                             components, binary, compress);

  bfam_subdomain_t **subdomains =
      bfam_malloc(numElements * sizeof(bfam_subdomain_t *));
  bfam_locidx_t numSubdomains;

  bfam_domain_get_subdomains(domain, match, tags, numElements, subdomains,
                             &numSubdomains);

  char filename[BFAM_BUFSIZ];
  if (directory)
    snprintf(filename, BFAM_BUFSIZ, "%s/" BFAM_VTK_VTU_FORMAT, directory,
             prefix, rank);
  else
    snprintf(filename, BFAM_BUFSIZ, BFAM_VTK_VTU_FORMAT, prefix, rank);

  BFAM_VERBOSE("Writing file: '%s'", filename);
  FILE *file = fopen(filename, "w");

  if (file == NULL)
  {
    BFAM_LERROR("Could not open %s for output!\n", filename);
    return;
  }

  fprintf(file, "<?xml version=\"1.0\"?>\n");
  fprintf(file, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"");

  if (binary && compress)
    fprintf(file, " compressor=\"vtkZLibDataCompressor\"");

  if (endian == BFAM_BIG_ENDIAN)
    fprintf(file, " byte_order=\"BigEndian\">\n");
  else
    fprintf(file, " byte_order=\"LittleEndian\">\n");

  fprintf(file, "  <UnstructuredGrid>\n");

  int files_written = 0;
  for (bfam_locidx_t s = 0; s < numSubdomains; ++s)
  {
    bfam_subdomain_t *subdomain = subdomains[s];

    if (subdomain->vtk_write_vtu_piece)
    {
      files_written += subdomain->vtk_write_vtu_piece(
          subdomains[s], file, time, scalars, vectors, components, binary,
          compress, rank, s, Np_write);
    }
    else
    {
      BFAM_WARNING("Subdomain: %s does not implement vtk_write_vtu_piece",
                   subdomain->name);
    }
  }

  if (files_written == 0)
    bfam_vtk_write_vtu_empty(file, binary);

  fprintf(file, "  </UnstructuredGrid>\n");
  fprintf(file, "</VTKFile>\n");

  if (ferror(file))
  {
    BFAM_LERROR("Error writing to %s\n", filename);
  }

  if (fclose(file))
  {
    BFAM_LERROR("Error closing %s\n", filename);
  }

  bfam_free(subdomains);
}

int bfam_vtk_write_binary_data(int compressed, FILE *file, char *data,
                               size_t size)
{
  if (compressed)
    return sc_vtk_write_compressed(file, data, size);
  else
    return sc_vtk_write_binary(file, data, size);
}

void bfam_vtk_write_real_scalar_data_array(FILE *file, const char *name,
                                           int writeBinary, int writeCompressed,
                                           bfam_locidx_t Ntotal,
                                           const bfam_real_t *s)
{
  const char *format;

  if (writeBinary)
    format = "binary";
  else
    format = "ascii";

  fprintf(file, "        <DataArray type=\"%s\" Name=\"%s\" format=\"%s\">\n",
          BFAM_REAL_VTK, name, format);
  if (writeBinary)
  {
    fprintf(file, "          ");
    int rval = bfam_vtk_write_binary_data(writeCompressed, file, (char *)s,
                                          Ntotal * sizeof(bfam_real_t));
    fprintf(file, "\n");
    if (rval)
      BFAM_WARNING("Error encoding %s", name);
  }
  else
  {
    for (bfam_locidx_t n = 0; n < Ntotal; ++n)
      fprintf(file, "         %" BFAM_REAL_FMTe "\n", s[n]);
  }

  fprintf(file, "        </DataArray>\n");
}

void bfam_vtk_write_real_vector_data_array(FILE *file, const char *name,
                                           int writeBinary, int writeCompressed,
                                           bfam_locidx_t Ntotal,
                                           const bfam_real_t *v1,
                                           const bfam_real_t *v2,
                                           const bfam_real_t *v3)
{
  const char *format;

  if (writeBinary)
    format = "binary";
  else
    format = "ascii";

  fprintf(file, "        <DataArray type=\"%s\" Name=\"%s\""
                " NumberOfComponents=\"3\" format=\"%s\">\n",
          BFAM_REAL_VTK, name, format);
  if (writeBinary)
  {
    size_t vSize = 3 * Ntotal * sizeof(bfam_real_t);
    bfam_real_t *v = bfam_malloc_aligned(vSize);

    for (bfam_locidx_t n = 0; n < Ntotal; ++n)
    {
      v[3 * n + 0] = (v1 != NULL) ? v1[n] : 0;
      v[3 * n + 1] = (v2 != NULL) ? v2[n] : 0;
      v[3 * n + 2] = (v3 != NULL) ? v3[n] : 0;
    }

    fprintf(file, "          ");
    int rval =
        bfam_vtk_write_binary_data(writeCompressed, file, (char *)v, vSize);
    fprintf(file, "\n");
    if (rval)
      BFAM_WARNING("Error encoding %s", name);

    bfam_free_aligned(v);
  }
  else
  {
    for (bfam_locidx_t n = 0; n < Ntotal; ++n)
    {
      bfam_real_t a = (v1 != NULL) ? v1[n] : 0;
      bfam_real_t b = (v2 != NULL) ? v2[n] : 0;
      bfam_real_t c = (v3 != NULL) ? v3[n] : 0;

      fprintf(file, "         %" BFAM_REAL_FMTe " %" BFAM_REAL_FMTe
                    " %" BFAM_REAL_FMTe "\n",
              a, b, c);
    }
  }

  fprintf(file, "        </DataArray>\n");
}

void bfam_vtk_write_struc_file(bfam_domain_t *domain, bfam_domain_match_t match,
                               const char **tags, const char *prefix,
                               const char **scalars, const char **vectors,
                               const char **components, int binary,
                               int compress)
{
  const int endian = bfam_endian();

  bfam_locidx_t numElements = domain->numSubdomains;

  int rank, size;
  BFAM_MPI_CHECK(MPI_Comm_rank(domain->comm, &rank));
  BFAM_MPI_CHECK(MPI_Comm_size(domain->comm, &size));

  bfam_subdomain_t **subdomains =
      bfam_malloc(numElements * sizeof(bfam_subdomain_t *));
  bfam_locidx_t numSubdomains = 0;

  bfam_domain_get_subdomains(domain, match, tags, numElements, subdomains,
                             &numSubdomains);

  for (bfam_locidx_t s = 0; s < numSubdomains; s++)
  {
    char filename[BFAM_BUFSIZ];
    char suffix[BFAM_BUFSIZ];
    bfam_subdomain_t *subdomain = subdomains[s];

    if (subdomain->vtk_write_suffix)
    {
      subdomain->vtk_write_suffix(subdomain, suffix, BFAM_BUFSIZ);
      snprintf(filename, BFAM_BUFSIZ, BFAM_VTK_VTS_FORMAT, prefix,
               subdomain->name, suffix);
    }
    else
    {
      snprintf(suffix, BFAM_BUFSIZ, "%d", rank);
      snprintf(filename, BFAM_BUFSIZ, BFAM_VTK_VTS_FORMAT, prefix,
               subdomain->name, suffix);
    }

    BFAM_VERBOSE("Writing file: '%s'", filename);
    FILE *file = fopen(filename, "w");
    if (file == NULL)
    {
      BFAM_LERROR("Could not open %s for output!\n", filename);
      return;
    }
    fprintf(file, "<?xml version=\"1.0\"?>\n");
    fprintf(file, "<VTKFile type=\"StructuredGrid\" version=\"0.1\"");

    if (binary && compress)
      fprintf(file, " compressor=\"vtkZLibDataCompressor\"");

    if (endian == BFAM_BIG_ENDIAN)
      fprintf(file, " byte_order=\"BigEndian\">\n");
    else
      fprintf(file, " byte_order=\"LittleEndian\">\n");

    if (subdomain->vtk_write_vts_piece)
      subdomain->vtk_write_vts_piece(subdomain, file, scalars, vectors,
                                     components, binary, compress, rank);
    else
      BFAM_WARNING("Subdomain: %s does not implement vtk_write_vts_piece",
                   subdomain->name);

    fprintf(file, "</VTKFile>\n");

    if (ferror(file))
    {
      BFAM_LERROR("Error writing to %s\n", filename);
    }

    if (fclose(file))
    {
      BFAM_LERROR("Error closing %s\n", filename);
    }
  }

  bfam_free(subdomains);
}
