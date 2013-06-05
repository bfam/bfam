#include <bfam_base.h>
#include <bfam_vtk.h>
#include <bfam_log.h>

void
bfam_vtk_write_file(bfam_domain_t *domain, bfam_domain_match_t match,
    const char **tags, const char *prefix, const char **scalars,
    const char **vectors, const char **components)
{
  bfam_locidx_t numElements = domain->numSubdomains;

  int rank;
  BFAM_MPI_CHECK(MPI_Comm_rank(domain->comm->comm, &rank));

  bfam_subdomain_t **subdomains =
    bfam_malloc(numElements*sizeof(bfam_subdomain_t*));
  bfam_locidx_t numSubdomains;

  bfam_domain_get_subdomains(domain, match, tags,
    numElements, subdomains, &numSubdomains);

  for(bfam_locidx_t s = 0; s < numSubdomains; ++s)
  {
    bfam_subdomain_t *subdomain = subdomains[s];

    char newPrefix[BFAM_BUFSIZ];
    snprintf(newPrefix, BFAM_BUFSIZ, "%s_%s_%05d", prefix, subdomain->name,
        rank);

    if(subdomain->vtk_write_file)
    {
      subdomain->vtk_write_file(subdomains[s], newPrefix, scalars, vectors,
          components);
    }
    else
    {
      BFAM_WARNING("Subdomain: %s does not implement vtk_write_file",
        subdomain->name);
    }
  }

  bfam_free(subdomains);
}
