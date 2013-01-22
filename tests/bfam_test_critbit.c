#include <bfam.h>

static void
test_contains()
{
  bfam_critbit0_tree_t tree = {0};

  static const char *elems[] =
      {"a", "aa", "b", "bb", "ab", "ba", "aba", "bab", NULL};

  for (unsigned i = 0; elems[i]; ++i)
    bfam_critbit0_insert(&tree, elems[i]);

  for (unsigned i = 0; elems[i]; ++i)
  {
    if (!bfam_critbit0_contains(&tree, elems[i]))
      BFAM_ABORT("Contains fail");
  }

  bfam_critbit0_clear(&tree);
}

static void
test_delete()
{
  bfam_critbit0_tree_t tree = {0};

  static const char *elems[] =
      {"a", "aa", "b", "bb", "ab", "ba", "aba", "bab", NULL};

  for (unsigned i = 1; elems[i]; ++i)
  {
    bfam_critbit0_clear(&tree);

    for (unsigned j = 0; j < i; ++j)
      bfam_critbit0_insert(&tree, elems[j]);

    for (unsigned j = 0; j < i; ++j)
    {
      if (!bfam_critbit0_contains(&tree, elems[j]))
        abort();
    }
    for (unsigned j = 0; j < i; ++j)
    {
      if (1 != bfam_critbit0_delete(&tree, elems[j]))
        abort();
    }
    for (unsigned j = 0; j < i; ++j)
    {
      if (bfam_critbit0_contains(&tree, elems[j]))
        abort();
    }
  }

  bfam_critbit0_clear(&tree);
}

int
main (int argc, char *argv[])
{
  test_contains();
  test_delete();

  return EXIT_SUCCESS;
}
