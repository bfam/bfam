#include <bfam.h>

int
prefix(const char * key, const char *val, void *arg)
{
  //printf("key: %s \nval: %s \n",key,val);
  const char **check = (const char **)arg;

  for (unsigned i = 0; check[i]; ++i)
    if(strcmp(val,check[i])) return 1;

  BFAM_ABORT("all prefix fail");

  return 0;
}

static void
test_contains()
{
  bfam_dictionary_t dict = {{0}};

  static const char *keys[] =
      {"a", "aa", "b", "bb", "ab", "ba", "aba", "bab", NULL};

  static const char *notkeys[] =
      {"_a", "_aa", "_b", "_bb", "_ab", "_ba", "_aba", "_bab", NULL};

  static const char *values[] =
      {"1", "22", "333", "4444", "55555", "666666", "7777777",
        "88888888", NULL};

  for (unsigned i = 0; keys[i]; ++i)
    bfam_dictionary_insert(&dict, keys[i], values[i]);


  for (unsigned i = 0; keys[i]; ++i)
  {
    if(1 != bfam_dictionary_insert(&dict, keys[i], values[i]))
      BFAM_ABORT("double insert fail");
  }

  for (unsigned i = 0; keys[i]; ++i)
  {
    if (!bfam_dictionary_contains(&dict, keys[i]))
      BFAM_ABORT("Contains fail");
  }

  for (unsigned i = 0; keys[i]; ++i)
  {
    char *val = bfam_dictionary_get_value(&dict, keys[i]);
    if(strcmp(values[i],val))
      BFAM_ABORT("Return is key fail");
  }

  for (unsigned i = 0; notkeys[i]; ++i)
  {
    if (NULL != bfam_dictionary_get_value(&dict, notkeys[i]))
      BFAM_ABORT("Return not key fail");
  }

  static const char *vals[] =
      {"666666", "88888888", NULL};
  bfam_dictionary_allprefixed(&dict,"ba",&prefix,&vals);


  bfam_dictionary_clear(&dict);
}

int
main (int argc, char *argv[])
{
  test_contains();

  return EXIT_SUCCESS;
}
