#include <bfam_base.h>
#include <bfam_util.h>

void
bfam_util_strcsl(char *str, const char **list)
{
  str[0] = '\0';

  if(list != NULL)
  {
    for(size_t l = 0; list[l]; ++l)
    {
      strcat(str, list[l]);
      if(list[l+1])
        strcat(str, ",");
    }
  }
}

