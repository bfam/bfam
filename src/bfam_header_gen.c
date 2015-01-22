#include <assert.h>
#include <ctype.h>
#include <libgen.h>
#include <stdio.h>
#include <string.h>

int main(int argc, char *argv[])
{
  int i;
  FILE *fin;
  FILE *fout;
  char *base;
  char *dot;
  char BASE[BUFSIZ];
  int c;

  if (argc < 3)
  {
    fprintf(stderr, "Usage: %s FILE HEADER_FILE\n", argv[0]);
    return 1;
  }

  fin = fopen(argv[1], "r");
  if (!fin)
  {
    fprintf(stderr, "Error: Cannot open %s for reading.\n", argv[1]);
    return 1;
  }

  fout = fopen(argv[2], "w");
  if (!fout)
  {
    fprintf(stderr, "Error: Cannot open %s for writing.\n", argv[2]);
    return 1;
  }

  base = basename(argv[2]);

  for (i = 0; base[i]; ++i)
  {
    BASE[i] = (char)toupper(base[i]);
    if (BASE[i] == '.')
      BASE[i] = '_';
  }
  BASE[i] = '\0';

  dot = strrchr(base, '.');
  if (dot != NULL)
    *dot = '\0';

  fprintf(fout, "#ifndef %s\n", BASE);
  fprintf(fout, "#define %s\n", BASE);

  fprintf(fout, "static const char %s[] = {\n", base);
  while ((c = getc(fin)) != EOF)
  {
    fprintf(fout, "0x%02x, ", (unsigned int)c);
  }
  fprintf(fout, "0x%02x\n};\n", '\0');

  fprintf(fout, "#endif\n");

  fclose(fout);
  fclose(fin);
  return 0;
}
