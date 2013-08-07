#include <bfam.h>

/**
 * TEST SUBDOMAIN TYPE
 */
typedef struct bfam_subdomain_lsrk_test
{
  bfam_subdomain_t base;
  bfam_locidx_t  rank; /* my rank */
  int               N; /* order of scheme to test */
} bfam_subdomain_lsrk_test_t;

static int
bfam_subdomain_lsrk_test_field_add(bfam_subdomain_t *subdomain, const char *name)
{
  bfam_subdomain_lsrk_test_t *s = (bfam_subdomain_lsrk_test_t*) subdomain;

  if(bfam_dictionary_get_value_ptr(&s->base.fields,name))
    return 1;

  bfam_real_t *field = bfam_malloc_aligned(sizeof(bfam_real_t));

  int rval = bfam_dictionary_insert_ptr(&s->base.fields, name, field);

  BFAM_ASSERT(rval != 1);

  if(rval == 0)
    bfam_free_aligned(field);

  return rval;
}

static int
bfam_subdomain_lsrk_test_free_fields(const char * key, void *val,
    void *arg)
{
  bfam_free_aligned(val);

  return 1;
}

void
bfam_subdomain_lsrk_test_free(bfam_subdomain_t *thisSubdomain)
{
  bfam_subdomain_lsrk_test_t *sub = (bfam_subdomain_lsrk_test_t*) thisSubdomain;

  bfam_dictionary_allprefixed_ptr(&sub->base.fields,"",
      &bfam_subdomain_lsrk_test_free_fields,NULL);
  bfam_dictionary_allprefixed_ptr(&sub->base.fields_p,"",
      &bfam_subdomain_lsrk_test_free_fields,NULL);
  bfam_dictionary_allprefixed_ptr(&sub->base.fields_m,"",
      &bfam_subdomain_lsrk_test_free_fields,NULL);
  bfam_dictionary_allprefixed_ptr(&sub->base.fields_face,"",
      &bfam_subdomain_lsrk_test_free_fields,NULL);
  bfam_subdomain_free(thisSubdomain);
}

bfam_long_real_t
bfam_subdomain_lsrk_test_exact(bfam_subdomain_lsrk_test_t *subdom, bfam_long_real_t t)
{
  return pow(t,(bfam_long_real_t) subdom->N);
}

void
bfam_subdomain_lsrk_test_init(bfam_subdomain_lsrk_test_t       *subdomain,
                             const bfam_locidx_t        id,
                             const char                *name,
                             const int                  N)
{
  bfam_subdomain_init(&subdomain->base, id, name);
  bfam_subdomain_add_tag(&subdomain->base, "_subdomain_lsrk_test");
  subdomain->base.field_add = bfam_subdomain_lsrk_test_field_add;
  subdomain->base.free = bfam_subdomain_lsrk_test_free;
  subdomain->N = N;
}


bfam_subdomain_lsrk_test_t*
bfam_subdomain_lsrk_test_new(const bfam_locidx_t     id,
                         const char             *name,
                         const int               N)
{
  bfam_subdomain_lsrk_test_t* newSubdomain =
    bfam_malloc(sizeof(bfam_subdomain_lsrk_test_t));

  bfam_subdomain_lsrk_test_init(newSubdomain, id, name, N);

  return newSubdomain;
}


/**
 * routine to create the rates
 */
static void
test_lsrk_aux_rates(bfam_subdomain_t *sub, const char *prefix)
{
  char f1[BFAM_BUFSIZ];
  snprintf(f1,BFAM_BUFSIZ,"%s_f1",prefix);
  int rval = sub->field_add(sub,f1);
  BFAM_ASSERT(rval != 1);

  char f3[BFAM_BUFSIZ];
  snprintf(f3,BFAM_BUFSIZ,"%s_f3",prefix);
  rval = sub->field_add(sub,f3);
  BFAM_ASSERT(rval != 1);
}

void
test_lsrk(bfam_ts_lsrk_method_t type,int N,MPI_Comm mpicomm, int commtag)
{
  bfam_ts_lsrk_t* ts;

  bfam_domain_t* dom = bfam_domain_new(mpicomm);
  bfam_subdomain_lsrk_test_t* subDom = NULL;

  subDom = bfam_subdomain_lsrk_test_new(1,"a",N);
  bfam_subdomain_add_tag((bfam_subdomain_t*)subDom,"yes");
  bfam_subdomain_add_tag((bfam_subdomain_t*)subDom,"other");
  bfam_subdomain_add_tag((bfam_subdomain_t*)subDom,"lsrk_test");
  bfam_domain_add_subdomain(dom, (bfam_subdomain_t*) subDom);

  subDom = bfam_subdomain_lsrk_test_new(2,"b",N);
  bfam_subdomain_add_tag((bfam_subdomain_t*)subDom,"no");
  bfam_subdomain_add_tag((bfam_subdomain_t*)subDom,"other");
  bfam_domain_add_subdomain(dom, (bfam_subdomain_t*) subDom);

  subDom = bfam_subdomain_lsrk_test_new(3,"c",N);
  bfam_subdomain_add_tag((bfam_subdomain_t*)subDom,"lsrk_test");
  bfam_subdomain_add_tag((bfam_subdomain_t*)subDom,"yes");
  bfam_subdomain_add_tag((bfam_subdomain_t*)subDom,"random");
  bfam_domain_add_subdomain(dom, (bfam_subdomain_t*) subDom);

  subDom = bfam_subdomain_lsrk_test_new(4,"d",N);
  bfam_subdomain_add_tag((bfam_subdomain_t*)subDom,"other");
  bfam_subdomain_add_tag((bfam_subdomain_t*)subDom,"yes");
  bfam_subdomain_add_tag((bfam_subdomain_t*)subDom,"random");
  bfam_domain_add_subdomain(dom, (bfam_subdomain_t*) subDom);


  const char *ts_tags[] = {NULL};
  const char *comm_tags[] = {"NOOP",NULL};

  bfam_domain_add_field(dom, BFAM_DOMAIN_AND, ts_tags, "f1");
  bfam_domain_add_field(dom, BFAM_DOMAIN_AND, ts_tags, "f2");
  bfam_domain_add_field(dom, BFAM_DOMAIN_AND, ts_tags, "f3");

  ts = bfam_ts_lsrk_new(dom,type,
      BFAM_DOMAIN_AND,ts_tags,
      BFAM_DOMAIN_AND,comm_tags,
      mpicomm,commtag,
      &test_lsrk_aux_rates,
      NULL,NULL,NULL,NULL);

  // free
  bfam_ts_lsrk_free(ts);
  bfam_free(ts);

  bfam_domain_free(dom);
  bfam_free(dom);
}

int
main (int argc, char *argv[])
{
  BFAM_MPI_CHECK(MPI_Init(&argc,&argv));
  test_lsrk(BFAM_TS_LSRK_KC54,4,MPI_COMM_WORLD,10);
  // test_lsrk(BFAM_TS_LSRK_W33 ,3,MPI_COMM_WORLD,10);
  // test_lsrk(BFAM_TS_LSRK_HEUN,2,MPI_COMM_WORLD,10);
  // test_lsrk(BFAM_TS_LSRK_FE  ,1,MPI_COMM_WORLD,10);
  BFAM_MPI_CHECK(MPI_Finalize());

  return EXIT_SUCCESS;
}
