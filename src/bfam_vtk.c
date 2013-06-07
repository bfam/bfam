#include <bfam_base.h>
#include <bfam_vtk.h>
#include <bfam_log.h>

#define BFAM_VTK_VTU_FORMAT "%s_%05d.vtu"

static void
bfam_vtk_write_pfile(int size, const char *prefix, const char **scalars,
    const char **vectors, const char **components, int binary,
    int compress)
{
  const int endian = bfam_endian();

  char filename[BFAM_BUFSIZ];
  snprintf(filename, BFAM_BUFSIZ, "%s.pvtu", prefix);

  const char* format;
  if(binary)
    format = "binary";
  else
    format = "ascii";


  BFAM_VERBOSE("Writing file: '%s'", filename);
  FILE *file = fopen(filename, "w");

  if(file == NULL)
  {
    BFAM_LERROR("Could not open %s for output!\n", filename);
    return;
  }

  fprintf(file, "<?xml version=\"1.0\"?>\n");
  fprintf(file, "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\"");

  if(binary && compress)
    fprintf(file, " compressor=\"vtkZLibDataCompressor\"");

  if(endian == BFAM_BIG_ENDIAN)
    fprintf(file, " byte_order=\"BigEndian\">\n");
  else
    fprintf(file, " byte_order=\"LittleEndian\">\n");

  fprintf(file, "  <PUnstructuredGrid GhostLevel=\"0\">\n");
  fprintf(file, "    <PPoints>\n");
  fprintf(file, "      <PDataArray type=\"%s\" Name=\"Position\""
          " NumberOfComponents=\"3\" format=\"%s\"/>\n",
          BFAM_REAL_VTK, format);
  fprintf(file, "    </PPoints>\n");

  fprintf(file, "    <PCellData Scalars=\"mpirank,subdomain_id\">\n");
  fprintf(file, "      <PDataArray type=\"%s\" Name=\"mpirank\""
           " format=\"%s\"/>\n", BFAM_LOCIDX_VTK, format);
  fprintf(file, "      <PDataArray type=\"%s\" Name=\"subdomain_id\""
           " format=\"%s\"/>\n", BFAM_LOCIDX_VTK, format);
  fprintf(file, "    </PCellData>\n");
  for(int s = 0; s < size; ++s)
  {
    fprintf(file, "    <Piece Source=\""BFAM_VTK_VTU_FORMAT"\"/>\n", prefix, s);
  }
  fprintf(file, "  </PUnstructuredGrid>\n");
  fprintf(file, "</VTKFile>\n");

  if(ferror(file))
  {
    BFAM_LERROR("Error writing to %s\n", filename);
  }

  if(fclose(file))
  {
    BFAM_LERROR("Error closing %s\n", filename);
  }

}

void
bfam_vtk_write_file(bfam_domain_t *domain, bfam_domain_match_t match,
    const char **tags, const char *prefix, const char **scalars,
    const char **vectors, const char **components, int binary,
    int compress)
{
  const int endian = bfam_endian();

  bfam_locidx_t numElements = domain->numSubdomains;

  int rank, size;
  BFAM_MPI_CHECK(MPI_Comm_rank(domain->comm->comm, &rank));
  BFAM_MPI_CHECK(MPI_Comm_size(domain->comm->comm, &size));

  if(rank == 0)
    bfam_vtk_write_pfile(size, prefix, scalars, vectors, components, binary,
        compress);

  bfam_subdomain_t **subdomains =
    bfam_malloc(numElements*sizeof(bfam_subdomain_t*));
  bfam_locidx_t numSubdomains;

  bfam_domain_get_subdomains(domain, match, tags,
    numElements, subdomains, &numSubdomains);

  char filename[BFAM_BUFSIZ];
  snprintf(filename, BFAM_BUFSIZ, BFAM_VTK_VTU_FORMAT, prefix, rank);

  BFAM_VERBOSE("Writing file: '%s'", filename);
  FILE *file = fopen(filename, "w");

  if(file == NULL)
  {
    BFAM_LERROR("Could not open %s for output!\n", filename);
    return;
  }

  fprintf(file, "<?xml version=\"1.0\"?>\n");
  fprintf(file, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"");

  if(binary && compress)
    fprintf(file, " compressor=\"vtkZLibDataCompressor\"");

  if(endian == BFAM_BIG_ENDIAN)
    fprintf(file, " byte_order=\"BigEndian\">\n");
  else
    fprintf(file, " byte_order=\"LittleEndian\">\n");

  fprintf(file, "  <UnstructuredGrid>\n");

  for(bfam_locidx_t s = 0; s < numSubdomains; ++s)
  {
    bfam_subdomain_t *subdomain = subdomains[s];

    if(subdomain->vtk_write_vtu_piece)
    {
      subdomain->vtk_write_vtu_piece(subdomains[s], file, scalars, vectors,
          components, binary, compress, rank, s);
    }
    else
    {
      BFAM_WARNING("Subdomain: %s does not implement vtk_write_vtu_piece",
        subdomain->name);
    }
  }

  fprintf(file, "  </UnstructuredGrid>\n");
  fprintf(file, "</VTKFile>\n");

  if(ferror(file))
  {
    BFAM_LERROR("Error writing to %s\n", filename);
  }

  if(fclose(file))
  {
    BFAM_LERROR("Error closing %s\n", filename);
  }

  bfam_free(subdomains);
}
