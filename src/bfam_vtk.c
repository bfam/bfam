#include <bfam_base.h>
#include <bfam_vtk.h>
#include <bfam_log.h>

void
bfam_vtk_write_file(bfam_domain_t *domain, bfam_domain_match_t match,
    const char **tags, const char *prefix, const char **scalars,
    const char **vectors, const char **components, int binary,
    int compress)
{
  const int endian = bfam_endian();

  bfam_locidx_t numElements = domain->numSubdomains;

  int rank;
  BFAM_MPI_CHECK(MPI_Comm_rank(domain->comm->comm, &rank));

  bfam_subdomain_t **subdomains =
    bfam_malloc(numElements*sizeof(bfam_subdomain_t*));
  bfam_locidx_t numSubdomains;

  bfam_domain_get_subdomains(domain, match, tags,
    numElements, subdomains, &numSubdomains);

  char filename[BFAM_BUFSIZ];
  snprintf(filename, BFAM_BUFSIZ, "%s_%05d.vtu", prefix, rank);

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
