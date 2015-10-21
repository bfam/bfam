#include <bfam_domain.h>
#include <bfam_base.h>
#include <bfam_log.h>

bfam_domain_t *bfam_domain_new(MPI_Comm domComm)
{
  bfam_domain_t *newDomain = bfam_malloc(sizeof(bfam_domain_t));
  bfam_domain_init(newDomain, domComm);
  return newDomain;
}

void bfam_domain_init(bfam_domain_t *thisDomain, MPI_Comm domComm)
{
  const bfam_locidx_t sizeSubdomains = 16;
  thisDomain->comm = domComm; // Perhaps we should duplicate it?
  thisDomain->num_subdomains = 0;
  thisDomain->sizeSubdomains = sizeSubdomains;
  thisDomain->subdomains =
      bfam_malloc(sizeSubdomains * sizeof(bfam_subdomain_t *));
  bfam_dictionary_init(&(thisDomain->name2num));
}

void bfam_domain_free(bfam_domain_t *thisDomain)
{
  for (bfam_locidx_t i = 0; i < thisDomain->num_subdomains; i++)
  {
    thisDomain->subdomains[i]->free(thisDomain->subdomains[i]);
    bfam_free(thisDomain->subdomains[i]);
  }
  thisDomain->comm = MPI_COMM_NULL;
  thisDomain->num_subdomains = 0;
  thisDomain->sizeSubdomains = 0;
  bfam_free(thisDomain->subdomains);
  thisDomain->subdomains = NULL;
  bfam_dictionary_clear(&thisDomain->name2num);
}

static int bfam_domain_compare_subdomain_by_id(const void *a, const void *b)
{
  const bfam_subdomain_t *subA = *(bfam_subdomain_t * const *)a;
  const bfam_subdomain_t *subB = *(bfam_subdomain_t * const *)b;

  int rval;

  if (subA->id < subB->id)
    rval = -1;
  else if (subA->id > subB->id)
    rval = 1;
  else
    rval = 0;

  return rval;
}

void bfam_domain_get_subdomains(bfam_domain_t *thisDomain,
                                bfam_domain_match_t matchType,
                                const char **tags, bfam_locidx_t numEntries,
                                bfam_subdomain_t **subdomains,
                                bfam_locidx_t *num_subdomains)
{
  BFAM_ASSERT(subdomains != NULL);
  BFAM_ASSERT(num_subdomains != NULL);

  if (numEntries <= 0)
    return;

  *num_subdomains = 0;
  for (bfam_locidx_t d = 0; d < thisDomain->num_subdomains; ++d)
  {
    bfam_subdomain_t *subdomain = thisDomain->subdomains[d];
    int matched = 0;
    switch (matchType)
    {
    case BFAM_DOMAIN_OR:
      matched = 0;
      for (size_t t = 0; !matched && tags[t]; ++t)
      {
        int hasTag = bfam_subdomain_has_tag(subdomain, tags[t]);
        matched = hasTag || matched;
      }
      break;
    case BFAM_DOMAIN_AND:
      matched = 1;
      for (size_t t = 0; matched && tags[t]; ++t)
        matched = matched && bfam_subdomain_has_tag(subdomain, tags[t]);
      break;
    default:
      BFAM_ABORT("Unsupported Match Type");
    }

    if (matched)
    {
      subdomains[*num_subdomains] = subdomain;
      ++(*num_subdomains);
    }

    if (*num_subdomains == numEntries)
      break;
  }

  /*
   * Sort subdomains by id number
   */
  qsort(subdomains, *num_subdomains, sizeof(bfam_subdomain_t *),
        bfam_domain_compare_subdomain_by_id);

  return;
}

static int bfam_domain_get_subdomains_critbit_or(const char *elem, void *arg)
{
  int *matched = (int *)((void **)arg)[0];
  bfam_subdomain_t *subdomain = (bfam_subdomain_t *)((void **)arg)[1];

  (*matched) = (*matched) || bfam_subdomain_has_tag(subdomain, elem);
  return !(*matched);
}

static int bfam_domain_get_subdomains_critbit_and(const char *elem, void *arg)
{
  int *matched = (int *)((void **)arg)[0];
  bfam_subdomain_t *subdomain = (bfam_subdomain_t *)((void **)arg)[1];

  (*matched) = (*matched) && bfam_subdomain_has_tag(subdomain, elem);
  return 1;
}

void bfam_domain_get_subdomains_critbit(bfam_domain_t *thisDomain,
                                        bfam_domain_match_t matchType,
                                        bfam_critbit0_tree_t *tags,
                                        bfam_locidx_t numEntries,
                                        bfam_subdomain_t **subdomains,
                                        bfam_locidx_t *num_subdomains)
{
  BFAM_ASSERT(subdomains != NULL);
  BFAM_ASSERT(num_subdomains != NULL);

  if (numEntries <= 0)
    return;

  *num_subdomains = 0;
  for (bfam_locidx_t d = 0; d < thisDomain->num_subdomains; ++d)
  {
    bfam_subdomain_t *subdomain = thisDomain->subdomains[d];
    int matched = 0;
    void *arg[2];
    arg[0] = &matched;
    arg[1] = subdomain;
    switch (matchType)
    {
    case BFAM_DOMAIN_OR:
      matched = 0;
      bfam_critbit0_allprefixed(tags, "",
                                &bfam_domain_get_subdomains_critbit_or, arg);
      break;
    case BFAM_DOMAIN_AND:
      matched = 1;
      bfam_critbit0_allprefixed(tags, "",
                                &bfam_domain_get_subdomains_critbit_and, arg);
      break;
    default:
      BFAM_ABORT("Unsupported Match Type");
    }

    if (matched)
    {
      subdomains[*num_subdomains] = subdomain;
      ++(*num_subdomains);
    }

    if (*num_subdomains == numEntries)
      break;
  }

  /*
   * Sort subdomains by id number
   */
  qsort(subdomains, *num_subdomains, sizeof(bfam_subdomain_t *),
        bfam_domain_compare_subdomain_by_id);

  return;
}

bfam_locidx_t bfam_domain_add_subdomain(bfam_domain_t *thisDomain,
                                        bfam_subdomain_t *newSubdomain)
{
  bfam_locidx_t sub_id;
  // double size
  if (thisDomain->num_subdomains == thisDomain->sizeSubdomains)
  {
    BFAM_VERBOSE("Doubling domain size");
    thisDomain->sizeSubdomains = 2 * thisDomain->sizeSubdomains;
    thisDomain->subdomains =
        bfam_realloc(thisDomain->subdomains,
                     thisDomain->sizeSubdomains * sizeof(bfam_subdomain_t *));
  }

  // create the key value pair
  BFAM_VERBOSE("adding subdomain %3" BFAM_LOCIDX_PRId " with name %s",
               thisDomain->num_subdomains, newSubdomain->name);

  // check if it's already there
  if (1 == bfam_dictionary_insert_locidx(&(thisDomain->name2num),
                                         newSubdomain->name,
                                         thisDomain->num_subdomains))
  {
    BFAM_ABORT("domain already contains subdomain \"%s\"", newSubdomain->name);
  }

  // add block
  sub_id = thisDomain->num_subdomains;
  thisDomain->subdomains[sub_id] = newSubdomain;
  thisDomain->num_subdomains++;
  return sub_id;
}

bfam_subdomain_t *bfam_domain_get_subdomain_by_num(bfam_domain_t *thisDomain,
                                                   bfam_locidx_t id)
{
  BFAM_ABORT_IF(id >= thisDomain->num_subdomains || id < 0,
                "Bad subdomain id: %jd", (intmax_t)id);
  return thisDomain->subdomains[id];
}

void bfam_domain_add_tag(bfam_domain_t *thisDomain, bfam_domain_match_t match,
                         const char **mtags, const char *tag)
{
  const char *tags[] = {tag, NULL};
  bfam_domain_add_tags(thisDomain, match, mtags, tags);
}

void bfam_domain_add_tags(bfam_domain_t *thisDomain, bfam_domain_match_t match,
                          const char **mtags, const char **tags)
{
  bfam_subdomain_t **subdomains =
      bfam_malloc(thisDomain->num_subdomains * sizeof(bfam_subdomain_t **));

  bfam_locidx_t num_subdomains = 0;

  bfam_domain_get_subdomains(thisDomain, match, mtags,
                             thisDomain->num_subdomains, subdomains,
                             &num_subdomains);

  for (bfam_locidx_t s = 0; s < num_subdomains; ++s)
    for (size_t t = 0; tags[t]; ++t)
      bfam_subdomain_add_tag(subdomains[s], tags[t]);

  bfam_free(subdomains);
}

void bfam_domain_add_field(bfam_domain_t *thisDomain, bfam_domain_match_t match,
                           const char **tags, const char *field)
{
  const char *fields[] = {field, NULL};
  bfam_domain_add_fields(thisDomain, match, tags, fields);
}

void bfam_domain_add_field_critbit(bfam_domain_t *thisDomain,
                                   bfam_domain_match_t match,
                                   bfam_critbit0_tree_t *tags,
                                   const char *field)
{
  const char *fields[] = {field, NULL};
  bfam_domain_add_fields_critbit(thisDomain, match, tags, fields);
}

void bfam_domain_add_fields(bfam_domain_t *thisDomain,
                            bfam_domain_match_t match, const char **tags,
                            const char **fields)
{
  bfam_subdomain_t **subdomains =
      bfam_malloc(thisDomain->num_subdomains * sizeof(bfam_subdomain_t **));

  bfam_locidx_t num_subdomains = 0;

  bfam_domain_get_subdomains(thisDomain, match, tags,
                             thisDomain->num_subdomains, subdomains,
                             &num_subdomains);

  for (bfam_locidx_t s = 0; s < num_subdomains; ++s)
    for (size_t f = 0; fields[f]; ++f)
      bfam_subdomain_field_add(subdomains[s], fields[f]);

  bfam_free(subdomains);
}

void bfam_domain_add_fields_critbit(bfam_domain_t *thisDomain,
                                    bfam_domain_match_t match,
                                    bfam_critbit0_tree_t *tags,
                                    const char **fields)
{
  bfam_subdomain_t **subdomains =
      bfam_malloc(thisDomain->num_subdomains * sizeof(bfam_subdomain_t **));

  bfam_locidx_t num_subdomains = 0;

  bfam_domain_get_subdomains_critbit(thisDomain, match, tags,
                                     thisDomain->num_subdomains, subdomains,
                                     &num_subdomains);

  for (bfam_locidx_t s = 0; s < num_subdomains; ++s)
    for (size_t f = 0; fields[f]; ++f)
      bfam_subdomain_field_add(subdomains[s], fields[f]);

  bfam_free(subdomains);
}

void bfam_domain_add_minus_field(bfam_domain_t *thisDomain,
                                 bfam_domain_match_t match, const char **tags,
                                 const char *field)
{
  const char *fields[] = {field, NULL};
  bfam_domain_add_minus_fields(thisDomain, match, tags, fields);
}

void bfam_domain_add_minus_field_critbit(bfam_domain_t *thisDomain,
                                         bfam_domain_match_t match,
                                         bfam_critbit0_tree_t *tags,
                                         const char *field)
{
  const char *fields[] = {field, NULL};
  bfam_domain_add_minus_fields_critbit(thisDomain, match, tags, fields);
}

void bfam_domain_add_minus_fields(bfam_domain_t *thisDomain,
                                  bfam_domain_match_t match, const char **tags,
                                  const char **fields)
{
  bfam_subdomain_t **subdomains =
      bfam_malloc(thisDomain->num_subdomains * sizeof(bfam_subdomain_t **));

  bfam_locidx_t num_subdomains = 0;

  bfam_domain_get_subdomains(thisDomain, match, tags,
                             thisDomain->num_subdomains, subdomains,
                             &num_subdomains);

  for (bfam_locidx_t s = 0; s < num_subdomains; ++s)
    for (size_t f = 0; fields[f]; ++f)
      bfam_subdomain_field_minus_add(subdomains[s], fields[f]);

  bfam_free(subdomains);
}

void bfam_domain_add_minus_fields_critbit(bfam_domain_t *thisDomain,
                                          bfam_domain_match_t match,
                                          bfam_critbit0_tree_t *tags,
                                          const char **fields)
{
  bfam_subdomain_t **subdomains =
      bfam_malloc(thisDomain->num_subdomains * sizeof(bfam_subdomain_t **));

  bfam_locidx_t num_subdomains = 0;

  bfam_domain_get_subdomains_critbit(thisDomain, match, tags,
                                     thisDomain->num_subdomains, subdomains,
                                     &num_subdomains);

  for (bfam_locidx_t s = 0; s < num_subdomains; ++s)
    for (size_t f = 0; fields[f]; ++f)
      bfam_subdomain_field_minus_add(subdomains[s], fields[f]);

  bfam_free(subdomains);
}

void bfam_domain_add_plus_field(bfam_domain_t *thisDomain,
                                bfam_domain_match_t match, const char **tags,
                                const char *field)
{
  const char *fields[] = {field, NULL};
  bfam_domain_add_plus_fields(thisDomain, match, tags, fields);
}

void bfam_domain_add_plus_field_critbit(bfam_domain_t *thisDomain,
                                        bfam_domain_match_t match,
                                        bfam_critbit0_tree_t *tags,
                                        const char *field)
{
  const char *fields[] = {field, NULL};
  bfam_domain_add_plus_fields_critbit(thisDomain, match, tags, fields);
}

void bfam_domain_add_plus_fields(bfam_domain_t *thisDomain,
                                 bfam_domain_match_t match, const char **tags,
                                 const char **fields)
{
  bfam_subdomain_t **subdomains =
      bfam_malloc(thisDomain->num_subdomains * sizeof(bfam_subdomain_t **));

  bfam_locidx_t num_subdomains = 0;

  bfam_domain_get_subdomains(thisDomain, match, tags,
                             thisDomain->num_subdomains, subdomains,
                             &num_subdomains);

  for (bfam_locidx_t s = 0; s < num_subdomains; ++s)
    for (size_t f = 0; fields[f]; ++f)
      bfam_subdomain_field_plus_add(subdomains[s], fields[f]);

  bfam_free(subdomains);
}

void bfam_domain_add_plus_fields_critbit(bfam_domain_t *thisDomain,
                                         bfam_domain_match_t match,
                                         bfam_critbit0_tree_t *tags,
                                         const char **fields)
{
  bfam_subdomain_t **subdomains =
      bfam_malloc(thisDomain->num_subdomains * sizeof(bfam_subdomain_t **));

  bfam_locidx_t num_subdomains = 0;

  bfam_domain_get_subdomains_critbit(thisDomain, match, tags,
                                     thisDomain->num_subdomains, subdomains,
                                     &num_subdomains);

  for (bfam_locidx_t s = 0; s < num_subdomains; ++s)
    for (size_t f = 0; fields[f]; ++f)
      bfam_subdomain_field_plus_add(subdomains[s], fields[f]);

  bfam_free(subdomains);
}

void bfam_domain_init_field(bfam_domain_t *thisDomain,
                            bfam_domain_match_t match, const char **tags,
                            const char *field, bfam_real_t time,
                            bfam_subdomain_init_field_t init_field, void *arg)
{
  bfam_subdomain_t **subdomains =
      bfam_malloc(thisDomain->num_subdomains * sizeof(bfam_subdomain_t **));

  bfam_locidx_t num_subdomains = 0;

  bfam_domain_get_subdomains(thisDomain, match, tags,
                             thisDomain->num_subdomains, subdomains,
                             &num_subdomains);

  for (bfam_locidx_t s = 0; s < num_subdomains; ++s)
    bfam_subdomain_field_init(subdomains[s], field, time, init_field, arg);

  bfam_free(subdomains);
}

void bfam_domain_init_fields_critbit(bfam_domain_t *thisDomain,
                                     bfam_domain_match_t match,
                                     bfam_critbit0_tree_t *tags,
                                     const char *field, bfam_real_t time,
                                     bfam_subdomain_init_field_t init_field,
                                     void *arg)
{
  bfam_subdomain_t **subdomains =
      bfam_malloc(thisDomain->num_subdomains * sizeof(bfam_subdomain_t **));

  bfam_locidx_t num_subdomains = 0;

  bfam_domain_get_subdomains_critbit(thisDomain, match, tags,
                                     thisDomain->num_subdomains, subdomains,
                                     &num_subdomains);

  for (bfam_locidx_t s = 0; s < num_subdomains; ++s)
    bfam_subdomain_field_init(subdomains[s], field, time, init_field, arg);

  bfam_free(subdomains);
}
