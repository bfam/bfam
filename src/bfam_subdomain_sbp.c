#include <bfam_subdomain_sbp.h>
#include <bfam_log.h>
#include <bfam_vtk.h>

bfam_subdomain_sbp_t*
bfam_subdomain_sbp_new(const bfam_locidx_t     id,
                            const char             *name,
                            const int               dim,
                            const bfam_gloidx_t    *N,
                            const bfam_locidx_t    *Nl,
                            const bfam_locidx_t    *Nb,
                            const bfam_gloidx_t    *gx)
{
  bfam_subdomain_sbp_t *newSub = bfam_malloc(sizeof(bfam_subdomain_sbp_t));
  bfam_subdomain_sbp_init(newSub,id,name,dim,N,Nl,Nb,gx);
  return newSub;
}

static int
bfam_subdomain_sbp_field_add(bfam_subdomain_t *subdomain, const char *name)
{
  bfam_subdomain_sbp_t *s = (bfam_subdomain_sbp_t*) subdomain;

  if(bfam_dictionary_get_value_ptr(&s->base.fields,name))
    return 1;

  size_t fieldSize = 0;
  for(int d = 0; d < s->dim; d++)
    fieldSize *= (s->Nl[d]+1+s->Nb[2*d]+s->Nb[2*d+1]);

  bfam_real_t *field = bfam_malloc_aligned(fieldSize);

  int rval = bfam_dictionary_insert_ptr(&s->base.fields, name, field);

  BFAM_ASSERT(rval != 1);

  if(rval == 0)
    bfam_free_aligned(field);

  return rval;
}

static void
bfam_subdomain_sbp_field_init(bfam_subdomain_t *subdomain,
    const char *name, bfam_real_t time, bfam_subdomain_init_field_t init_field,
    void *arg)
{
  bfam_subdomain_sbp_t *s = (bfam_subdomain_sbp_t*) subdomain;

  bfam_real_t *field = bfam_dictionary_get_value_ptr(&s->base.fields,name);

  BFAM_ABORT_IF(field==NULL, "Init: Field %s not found in subdomain %s",
      name, subdomain->name);

  // size_t fieldLength = s->Np*s->K;

  // init_field(fieldLength, time, s->x, s->y, s->z, subdomain, arg, field);
}

void
bfam_subdomain_sbp_init(bfam_subdomain_sbp_t *subdomain,
                            const bfam_locidx_t     id,
                            const char             *name,
                            const int               dim,
                            const bfam_gloidx_t    *N,
                            const bfam_locidx_t    *Nl,
                            const bfam_locidx_t    *Nb,
                            const bfam_gloidx_t    *gx)
{
  bfam_subdomain_init(&subdomain->base, id, name);
  bfam_subdomain_add_tag(&subdomain->base, "_subdomain_sbp");
  subdomain->base.free = bfam_subdomain_sbp_free;
  subdomain->base.field_add = bfam_subdomain_sbp_field_add;
  subdomain->base.field_init = bfam_subdomain_sbp_field_init;

  subdomain->dim = dim;

  subdomain->N  = bfam_malloc(  dim*sizeof(bfam_gloidx_t));
  subdomain->Nl = bfam_malloc(  dim*sizeof(bfam_locidx_t));
  subdomain->Nb = bfam_malloc(2*dim*sizeof(bfam_locidx_t));
  subdomain->gx = bfam_malloc(  dim*sizeof(bfam_gloidx_t));

  for(int d = 0; d < dim; d++)
  {
    subdomain->N [  d  ] = N [  d  ];
    subdomain->Nl[  d  ] = Nl[  d  ];
    subdomain->Nb[2*d  ] = Nb[2*d  ];
    subdomain->Nb[2*d+1] = Nb[2*d+1];
    subdomain->gx[  d  ] = gx[  d  ];

    BFAM_ABORT_IF(N [d] < 0 || Nl[d] < 0 || Nb[2*d] < 0 || Nb[2*d+1] < 0 ||
                  gx[d] < 0 || gx[d] > N[d] || gx[d]+Nl[d] > N[d],
                  "%s: problem with a dimension %d", d, subdomain->base.name);
  }
}

static int
bfam_subdomain_sbp_free_fields(const char * key, void *val,
    void *arg)
{
  bfam_free_aligned(val);

  return 1;
}

void
bfam_subdomain_sbp_free(bfam_subdomain_t *subdomain)
{
  bfam_subdomain_sbp_t *sub = (bfam_subdomain_sbp_t*) subdomain;

  bfam_dictionary_allprefixed_ptr(&sub->base.fields,"",
      &bfam_subdomain_sbp_free_fields,NULL);

  bfam_subdomain_free(subdomain);

  bfam_free(sub->N );
  bfam_free(sub->Nl);
  bfam_free(sub->Nb);
  bfam_free(sub->gx);
}
