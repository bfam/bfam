#include <bfam_base.h>
#include <bfam_communicator.h>
#include <bfam_log.h>

typedef struct bfam_communicator_map_entry
{
  bfam_locidx_t np; /* Neighbor's processor number */
  bfam_locidx_t ns; /* Neighbor's subdomain id */
  bfam_locidx_t ms; /* My subdomain id */
  size_t   send_sz; /* size of data I send */
  size_t   recv_sz; /* size of data I receive */
} bfam_communicator_map_entry_t;


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

  /* get the subdomains */
  bfam_subdomain_t *subdomains[domain->numSubdomains];
  bfam_domain_get_subdomains(domain, match, tags,
      domain->numSubdomains,subdomains,&communicator->numSubdomains);

  /* figure out the info for everyone */
  bfam_communicator_map_entry_t map[communicator->numSubdomains];
  for(int s = 0; s < communicator->numSubdomains;s++)
  {
    map[s].ms = subdomains[s]->id;
    BFAM_ABORT_IF(subdomains[s]->glue_comm_info == NULL,
        "glue_comm_info not initialized for sudbomain %s",
        subdomains[s]->name);
    subdomains[s]->glue_comm_info(subdomains[s],&map[s].np,&map[s].ns,
        &map[s].send_sz,&map[s].recv_sz);
  }
}

void
bfam_communicator_free(bfam_communicator_t *communicator)
{
  BFAM_VERBOSE("Communicator Free");
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
