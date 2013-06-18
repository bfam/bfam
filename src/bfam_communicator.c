#include <bfam_base.h>
#include <bfam_communicator.h>
#include <bfam_log.h>

typedef struct bfam_communicator_map_entry
{
  bfam_locidx_t np; /* Neighbor's processor number */
  bfam_locidx_t ns; /* Neighbor's subdomain id */
  bfam_locidx_t ms; /* My subdomain id */
  bfam_locidx_t rank; /*my rank*/
  size_t   send_sz; /* size of data I send */
  size_t   recv_sz; /* size of data I receive */
  bfam_subdomain_t* subdomain; /*pointer to subdomain*/
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

static int
bfam_communicator_send_compare(const void *a, const void *b)
{
  const bfam_communicator_map_entry_t *mapA =
    (const bfam_communicator_map_entry_t *)a;
  const bfam_communicator_map_entry_t *mapB =
    (const bfam_communicator_map_entry_t *)b;

  int rval = 0;

  if(mapA->np < mapB->np)
    rval = -1;
  else if(mapA->np > mapB->np)
    rval =  1;
  else if(mapA->ns > mapB->ns)
    rval = -1;
  else if(mapA->ns < mapB->ns)
    rval =  1;
  else if(mapA->ms > mapB->ms)
    rval = -1;
  else if(mapA->ms < mapB->ms)
    rval =  1;
  else if(mapA->subdomain->id > mapB->subdomain->id)
  {
    BFAM_ABORT_IF_NOT(mapA->np == mapA->rank,
        "if comparing sub id must be local comm");
    rval = -1;
  }
  else if(mapA->subdomain->id < mapB->subdomain->id)
  {
    BFAM_ABORT_IF_NOT(mapA->np == mapA->rank,
        "if comparing sub id must be local comm");
    rval =  1;
  }
  BFAM_ABORT_IF_NOT(rval,"Should not be same map!");

  return rval;
}

static int
bfam_communicator_recv_compare(const void *a, const void *b)
{
  const bfam_communicator_map_entry_t *mapA =
    (const bfam_communicator_map_entry_t *)a;
  const bfam_communicator_map_entry_t *mapB =
    (const bfam_communicator_map_entry_t *)b;

  int rval = 0;

  if(mapA->np < mapB->np)
    rval = -1;
  else if(mapA->np > mapB->np)
    rval =  1;
  else if(mapA->ms > mapB->ms)
    rval = -1;
  else if(mapA->ms < mapB->ms)
    rval =  1;
  else if(mapA->ns > mapB->ns)
    rval = -1;
  else if(mapA->ns < mapB->ns)
    rval =  1;

  return rval;
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

  size_t send_sz = 0;
  size_t recv_sz = 0;

  /* figure out the info for everyone */
  bfam_communicator_map_entry_t map[communicator->numSubdomains];
  for(int s = 0; s < communicator->numSubdomains;s++)
  {
    map[s].subdomain = subdomains[s];
    BFAM_ABORT_IF(subdomains[s]->glue_comm_info == NULL,
        "glue_comm_info not initialized for sudbomain %s",
        subdomains[s]->name);
    subdomains[s]->glue_comm_info(subdomains[s],&map[s].np,
        &map[s].ms,&map[s].ns,&map[s].send_sz,&map[s].recv_sz);
    BFAM_MPI_CHECK(MPI_Comm_rank(comm, &map[s].rank));
    send_sz += map[s].send_sz;
    recv_sz += map[s].recv_sz;
  }

  communicator->send_buf = bfam_malloc(send_sz);
  communicator->recv_buf = bfam_malloc(recv_sz);

  /* sort for send */
  qsort((void*) map, communicator->numSubdomains,
      sizeof(bfam_communicator_map_entry_t), bfam_communicator_send_compare);
  for(int s = 0; s < communicator->numSubdomains;s++)
  {
  }



  /* sort for recv */
  qsort((void*) map, communicator->numSubdomains,
      sizeof(bfam_communicator_map_entry_t), bfam_communicator_recv_compare);
}

void
bfam_communicator_free(bfam_communicator_t *communicator)
{
  BFAM_VERBOSE("Communicator Free");
  bfam_free(communicator->send_buf);
  bfam_free(communicator->recv_buf);
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
