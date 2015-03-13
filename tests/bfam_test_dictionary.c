#include <bfam.h>

int prefix(const char *key, const char *val, void *arg)
{
  // printf("key: %s \nval: %s \n",key,val);
  const char **check = (const char **)arg;

  for (unsigned i = 0; check[i]; ++i)
    if (strcmp(val, check[i]))
      return 1;

  BFAM_ABORT("all prefix fail");

  return 0;
}

int prefix_ptr(const char *key, void *val, void *arg)
{
  const void **check = (const void **)arg;

  // printf("%s: %p\n",key,val);
  for (unsigned i = 0; check[i]; ++i)
    if (val == check[i])
      return 1;

  BFAM_ABORT("ptr all prefix fail");

  return 0;
}

static void test_contains()
{
  bfam_dictionary_t dict;
  bfam_dictionary_init(&dict);

  static const char *keys[] = {"a",  "aa",  "b",   "bb", "ab",
                               "ba", "aba", "bab", NULL};

  static const char *notkeys[] = {"_a",  "_aa",  "_b",   "_bb", "_ab",
                                  "_ba", "_aba", "_bab", NULL};

  static const char *values[] = {
      "1", "22", "333", "4444", "55555", "666666", "7777777", "88888888", NULL};

  for (unsigned i = 0; keys[i]; ++i)
    bfam_dictionary_insert(&dict, keys[i], values[i]);

  BFAM_ABORT_IF(dict.num_entries != 8, "Not all key,value pairs added");

  for (unsigned i = 0; keys[i]; ++i)
  {
    if (1 != bfam_dictionary_insert(&dict, keys[i], values[i]))
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
    if (strcmp(values[i], val))
      BFAM_ABORT("Return is key fail");
  }

  for (unsigned i = 0; notkeys[i]; ++i)
  {
    if (NULL != bfam_dictionary_get_value(&dict, notkeys[i]))
      BFAM_ABORT("Return not key fail");
  }

  static const char *vals[] = {"666666", "88888888", NULL};
  bfam_dictionary_allprefixed(&dict, "ba", &prefix, &vals);

  bfam_dictionary_clear(&dict);
}

static void test_contains_ptr()
{
  bfam_dictionary_t dict;
  bfam_dictionary_init(&dict);

  static const char *keys[] = {"@^1", "@2", "^3", "@^44", NULL};

  static const char *values[] = {"1_", "2_", "3_", "44_", NULL};

  for (unsigned i = 0; keys[i]; ++i)
    bfam_dictionary_insert_ptr(&dict, keys[i], &values[i]);

  for (unsigned i = 0; keys[i]; ++i)
  {
    if (1 != bfam_dictionary_insert_ptr(&dict, keys[i], &values[i]))
      BFAM_ABORT("double insert fail");
  }

  for (unsigned i = 0; keys[i]; ++i)
  {
    if (!bfam_dictionary_contains(&dict, keys[i]))
      BFAM_ABORT("Contains fail");
  }

  for (unsigned i = 0; keys[i]; ++i)
  {
    void *val = bfam_dictionary_get_value_ptr(&dict, keys[i]);
    if (val != &values[i])
      BFAM_ABORT("Return is key fail");
  }

  static const void *vals[] = {&values[0], &values[3]};
  bfam_dictionary_allprefixed_ptr(&dict, "@^", &prefix_ptr, &vals);

  bfam_dictionary_clear(&dict);
}

static void test_contains_locidx()
{
  bfam_dictionary_t dict;
  bfam_dictionary_init(&dict);

  static const char *keys[] = {"@^1", "@2", "^3", "@^44", NULL};

  static const bfam_locidx_t values[] = {1, 2, 3, 44};

  for (unsigned i = 0; keys[i]; ++i)
    bfam_dictionary_insert_locidx(&dict, keys[i], values[i]);

  for (unsigned i = 0; keys[i]; ++i)
  {
    if (1 != bfam_dictionary_insert_locidx(&dict, keys[i], values[i]))
      BFAM_ABORT("double insert fail");
  }

  for (unsigned i = 0; keys[i]; ++i)
  {
    if (!bfam_dictionary_contains(&dict, keys[i]))
      BFAM_ABORT("Contains fail");
  }

  for (unsigned i = 0; keys[i]; ++i)
  {
    bfam_locidx_t val;
    int rval = bfam_dictionary_get_value_locidx(&dict, keys[i], &val);
    if (rval != 1 || val != values[i])
      BFAM_ABORT("Return is key fail %d %d %d", rval, val, values[i]);
  }

  bfam_dictionary_clear(&dict);
}

static void test_contains_delete()
{
  bfam_dictionary_t dict;
  bfam_dictionary_init(&dict);

  static const char *keys[] = {"@^1", "@2", "^3", "@^44", NULL};

  static const int values[] = {1, 2, 3, 44};

  for (unsigned i = 0; keys[i]; ++i)
    bfam_dictionary_insert_int(&dict, keys[i], values[i]);

  for (unsigned i = 0; keys[i]; ++i)
  {
    if (1 != bfam_dictionary_delete(&dict, keys[i]))
      BFAM_ABORT("delete fail");
  }

  if (0 != bfam_dictionary_delete(&dict, "crap"))
    BFAM_ABORT("non delete fail");

  bfam_dictionary_clear(&dict);
}

int main(int argc, char *argv[])
{
  test_contains();
  test_contains_ptr();
  test_contains_locidx();

  return EXIT_SUCCESS;
}
