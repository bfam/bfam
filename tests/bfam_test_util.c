#include <bfam.h>
#include <bfam_util.h>


int
main (int argc, char *argv[])
{
  int failures = 0;

  const char *strings[] = {"test1", "test2", "test3", NULL};
  const char *list = "test1,test2,test3";
  char join[BFAM_BUFSIZ];
  bfam_util_strcsl(join, strings);

  failures += strcmp(list, join);

  if(failures)
    return EXIT_FAILURE;
  else
    return EXIT_SUCCESS;
}
