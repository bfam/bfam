#include <bfam.h>
#include <bfam_dictionary.h>

void bfam_dictionary_init(bfam_dictionary_t *d)
{
  d->t.root = NULL;
}

void bfam_dictionary_clear(bfam_dictionary_t *d)
{
  bfam_critbit0_clear(&(d->t));
}

int bfam_dictionary_contains_check(const char * value, void * arg)
{
  (*(int*)arg)++;
  return 1;
}

int bfam_dictionary_contains(bfam_dictionary_t *d, const char *key)
{
  const size_t keylen = strlen(key);

  char *u = bfam_malloc(sizeof(char)*(keylen+2));

  memcpy(u, key, keylen);
  u[keylen] = BFAM_KEYVALUE_SPLIT;
  u[keylen+1] = '\0';
  int found = 0;
  int rval = bfam_critbit0_allprefixed(&(d->t),u,
      &bfam_dictionary_contains_check,&found);
  bfam_free(u);
  return found;
}

int bfam_dictionary_insert(bfam_dictionary_t *d, const char *key,
                                  const char *val)
{
  const size_t keylen = strlen(key);
  const size_t vallen = strlen(val);

  char *keyval = bfam_malloc(sizeof(char)*(keylen+vallen+2));

  memcpy(keyval, key, keylen);

  if(bfam_dictionary_contains(d,key))
  {
    free(keyval);
    return 1;
  }

  keyval[keylen]   = BFAM_KEYVALUE_SPLIT;
  memcpy(&keyval[keylen+1], val, vallen);
  keyval[keylen+vallen+1] = '\0';

  int rval = bfam_critbit0_insert(&(d->t), keyval);

  bfam_free(keyval);
  return rval;
}

int bfam_dictionary_get_value_handle(const char * keyval, void * arg)
{
  char * key = (char  *)((void **)arg)[0];
  char** val = (char **)((void **)arg)[1];
  int keylen = strlen(key);

  *val = &keyval[keylen];

  return 1;
}


char* bfam_dictionary_get_value(bfam_dictionary_t *d, const char *key)
{
  if(!bfam_dictionary_contains(d,key))
    return NULL;

  const size_t keylen = strlen(key);

  char *u = bfam_malloc(sizeof(char)*(keylen+2));

  memcpy(u, key, keylen);
  u[keylen] = BFAM_KEYVALUE_SPLIT;
  u[keylen+1] = '\0';

  char *value = NULL;
  void *arg[2];
  arg[0] = u;
  arg[1] = &value;

  int rval = bfam_critbit0_allprefixed(&(d->t),u,
      &bfam_dictionary_get_value_handle,arg);

  bfam_free(u);

  return value;
}
