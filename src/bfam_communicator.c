#include <bfam_base.h>
#include <bfam_communicator.h>
#include <bfam_log.h>

#define BFAM_COMM_NUM_SRT (6)
typedef struct bfam_communicator_map_entry
{
  bfam_locidx_t np; /* Neighbor's processor number */
  bfam_locidx_t s[BFAM_COMM_NUM_SRT];
  /* sort param 0: typically Neighbor's subdomain id */
  /* sort param 1: typically My subdomain id */
  /* sort param 2: typically local ID (sbp) */
  /* sort param 3: typically fce ID (sbp) */
  bfam_locidx_t rank;          /*my rank*/
  bfam_subdomain_t *subdomain; /*pointer to subdomain*/

  bfam_locidx_t orig_order; /* original order */
} bfam_communicator_map_entry_t;

bfam_communicator_t *bfam_communicator_new(bfam_domain_t *domain,
                                           bfam_domain_match_t match,
                                           const char **tags, MPI_Comm comm,
                                           int tag, void *user_args)
{
  bfam_communicator_t *newCommunicator =
      bfam_malloc(sizeof(bfam_communicator_t));

  bfam_communicator_init(newCommunicator, domain, match, tags, comm, tag,
                         user_args);
  return newCommunicator;
}

static int bfam_communicator_send_compare(const void *a, const void *b)
{
  const bfam_communicator_map_entry_t *mapA =
      (const bfam_communicator_map_entry_t *)a;
  const bfam_communicator_map_entry_t *mapB =
      (const bfam_communicator_map_entry_t *)b;

  if (mapA->np < mapB->np)
    return -1;
  if (mapA->np > mapB->np)
    return 1;
  for (int i = 0; i < BFAM_COMM_NUM_SRT; i++)
  {
    if (mapA->s[i] < mapB->s[i])
      return -1;
    if (mapA->s[i] > mapB->s[i])
      return 1;
  }
  BFAM_ABORT_IF(0 == 0, "Should not be same map!");

  return 0;
}

static int bfam_communicator_recv_compare(const void *a, const void *b)
{
  const bfam_communicator_map_entry_t *mapA =
      (const bfam_communicator_map_entry_t *)a;
  const bfam_communicator_map_entry_t *mapB =
      (const bfam_communicator_map_entry_t *)b;

  if (mapA->np < mapB->np)
    return -1;
  if (mapA->np > mapB->np)
    return 1;
  for (int i = 0; i < BFAM_COMM_NUM_SRT / 2; i++)
  {
    if (mapA->s[2 * i + 1] < mapB->s[2 * i + 1])
      return -1;
    if (mapA->s[2 * i + 1] > mapB->s[2 * i + 1])
      return 1;
    if (mapA->s[2 * i] < mapB->s[2 * i])
      return -1;
    if (mapA->s[2 * i] > mapB->s[2 * i])
      return 1;
  }
  BFAM_ABORT_IF(0 == 0, "Should not be same map!");

  return 0;
}

void bfam_communicator_init(bfam_communicator_t *communicator,
                            bfam_domain_t *domain, bfam_domain_match_t match,
                            const char **tags, MPI_Comm comm, int tag,
                            void *user_args)
{
  BFAM_LDEBUG("Communicator Init");
  communicator->comm = comm;
  communicator->tag = tag;

  /* get the subdomains */
  bfam_subdomain_t *subdomains[domain->numSubdomains];
  bfam_domain_get_subdomains(domain, match, tags, domain->numSubdomains,
                             subdomains, &communicator->num_subs);

  communicator->sub_data =
      bfam_malloc(communicator->num_subs * sizeof(bfam_comm_subdata_t));

  size_t send_sz = 0;
  size_t recv_sz = 0;
  communicator->num_procs = 0;
  communicator->user_args = user_args;

  /* figure out the info for everyone */
  bfam_communicator_map_entry_t map[communicator->num_subs];
  bfam_critbit0_tree_t procs = {0};
  char procStr[BFAM_BUFSIZ];
  for (int s = 0; s < communicator->num_subs; s++)
  {
    map[s].subdomain = subdomains[s];

    BFAM_ABORT_IF(subdomains[s]->glue_comm_info == NULL,
                  "glue_comm_info not initialized for sudbomain %s",
                  subdomains[s]->name);

    for (int i = 0; i < BFAM_COMM_NUM_SRT; i++)
      map[s].s[i] = 0;
    subdomains[s]->glue_comm_info(
        subdomains[s], &map[s].np, map[s].s, BFAM_COMM_NUM_SRT,
        &communicator->sub_data[s].send_sz, &communicator->sub_data[s].recv_sz,
        communicator->user_args);

    BFAM_MPI_CHECK(MPI_Comm_rank(comm, &map[s].rank));

    send_sz += communicator->sub_data[s].send_sz;
    recv_sz += communicator->sub_data[s].recv_sz;

    communicator->sub_data[s].subdomain = subdomains[s];

    communicator->sub_data[s].send_buf = NULL;
    communicator->sub_data[s].recv_buf = NULL;

    map[s].orig_order = s;

    snprintf(procStr, BFAM_BUFSIZ, "%jd", (intmax_t)map[s].np);
    if (!bfam_critbit0_contains(&procs, procStr))
    {
      communicator->num_procs++;
      bfam_critbit0_insert(&procs, procStr);
    }
  }
  bfam_critbit0_clear(&procs);

  /* allocate everything now */
  communicator->send_sz = send_sz;
  communicator->send_buf = bfam_malloc(communicator->send_sz);

  communicator->recv_sz = recv_sz;
  communicator->recv_buf = bfam_malloc(communicator->recv_sz);

  communicator->proc_data =
      bfam_malloc(communicator->num_procs * sizeof(bfam_comm_procdata_t));

  communicator->send_request =
      bfam_malloc(2 * communicator->num_procs * sizeof(MPI_Request));
  communicator->recv_request =
      communicator->send_request + communicator->num_procs;

  communicator->send_status =
      bfam_malloc(2 * communicator->num_procs * sizeof(MPI_Status));
  communicator->recv_status =
      communicator->send_status + communicator->num_procs;

  for (int i = 0; communicator->num_procs > i; i++)
  {
    communicator->proc_data[i].send_sz = 0;
    communicator->proc_data[i].send_buf = NULL;

    communicator->proc_data[i].recv_sz = 0;
    communicator->proc_data[i].recv_buf = NULL;

    communicator->send_request[i] = MPI_REQUEST_NULL;
    communicator->recv_request[i] = MPI_REQUEST_NULL;
  }

  /* sort for send  and fill struct */
  qsort((void *)map, communicator->num_subs,
        sizeof(bfam_communicator_map_entry_t), bfam_communicator_send_compare);
  char *send_buf_ptr = communicator->send_buf;
  bfam_locidx_t np = -1;   /* global proc ID */
  bfam_locidx_t proc = -1; /* local storage proc ID */
  for (int s = 0; s < communicator->num_subs; s++)
  {
    bfam_locidx_t t = map[s].orig_order;
    communicator->sub_data[t].send_buf = send_buf_ptr;

    if (map[s].np != np)
    {
      np = map[s].np;
      proc++;
      communicator->proc_data[proc].send_buf = send_buf_ptr;
      communicator->proc_data[proc].rank = np;
    }
    BFAM_ABORT_IF_NOT(communicator->proc_data[proc].rank == np,
                      "problem with local proc ID");
    communicator->proc_data[proc].send_sz += communicator->sub_data[t].send_sz;

    send_buf_ptr += communicator->sub_data[t].send_sz;
  }

  /* sort for recv and fill struct*/
  qsort((void *)map, communicator->num_subs,
        sizeof(bfam_communicator_map_entry_t), bfam_communicator_recv_compare);
  char *recv_buf_ptr = communicator->recv_buf;
  np = -1;   /* global proc ID */
  proc = -1; /* local storage proc ID */
  for (int s = 0; s < communicator->num_subs; s++)
  {
    bfam_locidx_t t = map[s].orig_order;
    communicator->sub_data[t].recv_buf = recv_buf_ptr;
    if (map[s].np != np)
    {
      np = map[s].np;
      proc++;
      communicator->proc_data[proc].recv_buf = recv_buf_ptr;
    }
    BFAM_ABORT_IF_NOT(communicator->proc_data[proc].rank == np,
                      "problem with local proc ID");
    communicator->proc_data[proc].recv_sz += communicator->sub_data[t].recv_sz;

    recv_buf_ptr += communicator->sub_data[t].recv_sz;
  }
}

void bfam_communicator_free(bfam_communicator_t *communicator)
{
  BFAM_LDEBUG("Communicator Free");

  /* Just make sure there are no pending requests */
  BFAM_MPI_CHECK(MPI_Waitall(2 * communicator->num_procs,
                             communicator->send_request,
                             communicator->send_status));

  bfam_free(communicator->send_buf);
  bfam_free(communicator->recv_buf);
  bfam_free(communicator->sub_data);
  bfam_free(communicator->proc_data);
  bfam_free(communicator->send_request);
  bfam_free(communicator->send_status);
}

void bfam_communicator_start(bfam_communicator_t *comm)
{
  BFAM_LDEBUG("Communicator Start");
  BFAM_MPI_CHECK(
      MPI_Waitall(2 * comm->num_procs, comm->send_request, comm->send_status));

  /* post recvs */
  for (int p = 0; p < comm->num_procs; p++)
  {
    BFAM_MPI_CHECK(MPI_Irecv(comm->proc_data[p].recv_buf,
                             comm->proc_data[p].recv_sz, MPI_BYTE,
                             comm->proc_data[p].rank, comm->tag, comm->comm,
                             &comm->recv_request[p]));
  }

  /* fill the data */
  for (int s = 0; s < comm->num_subs; s++)
  {
    bfam_subdomain_t *sub = comm->sub_data[s].subdomain;
    BFAM_ABORT_IF(sub->glue_put_send_buffer == NULL,
                  "glue_put_send_buffer is NULL for %s", sub->name);
    sub->glue_put_send_buffer(sub, comm->sub_data[s].send_buf,
                              comm->sub_data[s].send_sz, comm->user_args);
  }

  /* post sends */
  for (int p = 0; p < comm->num_procs; p++)
  {
    BFAM_MPI_CHECK(MPI_Isend(comm->proc_data[p].send_buf,
                             comm->proc_data[p].send_sz, MPI_BYTE,
                             comm->proc_data[p].rank, comm->tag, comm->comm,
                             &comm->send_request[p]));
  }
}

void bfam_communicator_finish(bfam_communicator_t *comm)
{
  BFAM_LDEBUG("Communicator Finish");
  BFAM_MPI_CHECK(
      MPI_Waitall(comm->num_procs, comm->recv_request, comm->recv_status));

  /* get the data */
  for (int s = 0; s < comm->num_subs; s++)
  {
    bfam_subdomain_t *sub = comm->sub_data[s].subdomain;
    BFAM_ABORT_IF(sub->glue_get_recv_buffer == NULL,
                  "glue_get_recv_buffer is NULL for %s", sub->name);
    sub->glue_get_recv_buffer(sub, comm->sub_data[s].recv_buf,
                              comm->sub_data[s].recv_sz, comm->user_args);
  }
}
