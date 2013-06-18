#include <bfam_base.h>
#include <bfam_communicator.h>
#include <bfam_log.h>

bfam_communicator_t*
bfam_communicator_new(bfam_domain_t *domain, bfam_domain_match_t match,
    const char **tags, MPI_Comm comm)
{
  bfam_communicator_t* newCommunicator =
    bfam_malloc(sizeof(bfam_communicator_t));

  bfam_communicator_init(newCommunicator, domain, match, tags, comm);
  return newCommunicator;
}

void
bfam_communicator_init(bfam_communicator_t* communicator,
    bfam_domain_t *domain, bfam_domain_match_t match, const char **tags,
    MPI_Comm comm)
{
  BFAM_VERBOSE("Communicator Init");
  communicator-> comm = comm;

  bfam_subdomain_t *subdomains[domain->numSubdomains];
  bfam_domain_get_subdomains(domain, match, tags,
      domain->numSubdomains,subdomains,&communicator->numSubdomains);
  communicator->subdomains =
    bfam_malloc(communicator->numSubdomains*sizeof(bfam_subdomain_t*));
}

void
bfam_communicator_free(bfam_communicator_t *communicator)
{
  BFAM_VERBOSE("Communicator Free");
  bfam_free(communicator->subdomains);
}

void
bfam_communicator_start(bfam_communicator_t *communicator)
{
  BFAM_VERBOSE("Communicator Start");
}

void
bfam_communicator_finish(bfam_communicator_t *communicator)
{
  BFAM_VERBOSE("Communicator Finish");
}
