#include <bfam_log.h>

static const char *const bfam_log_level_names[] = {
    "ALWAYS", "TRACE", "DEBUG", " VERB", " INFO", " WARN", "ERROR", "SILENT",
};

static const int bfam_log_root_rank = 0;
static int bfam_log_rank = 0;
static FILE *bfam_log_stream = NULL;
static int bfam_log_threshold = BFAM_LL_ALWAYS;

void bfam_log_init(int rank, FILE *stream, int threshold)
{

  bfam_log_rank = rank;
  if (stream == NULL)
  {
    bfam_log_stream = stdout;
  }
  else
  {
    bfam_log_stream = stream;
  }

  if (threshold == BFAM_LL_DEFAULT)
  {
    bfam_log_threshold = BFAM_LL_INFO;
  }
  else
  {
    BFAM_ABORT_IF_NOT(threshold <= BFAM_LL_SILENT &&
                          threshold >= BFAM_LL_ALWAYS,
                      "Invalid logging threshold");
    bfam_log_threshold = threshold;
  }
}

void bfam_log_printf(const char *file, int line, int category, int level, ...)
{
  FILE *stream = (bfam_log_stream) ? bfam_log_stream : stdout;
  va_list ap;
  const char *fmt;

  BFAM_ASSERT(level <= BFAM_LL_SILENT && level >= BFAM_LL_ALWAYS);

  if (bfam_log_threshold > level)
    return;

  if (category == BFAM_LC_ROOT && bfam_log_rank != bfam_log_root_rank)
    return;

  if (category == BFAM_LC_ALL)
  {
    fprintf(stream, "[%3d] ", bfam_log_rank);
  }
  else
  {
    fprintf(stream, "[   ] ");
  }
  fprintf(stream, "%s: ", bfam_log_level_names[level]);

  if (level == BFAM_LL_TRACE)
  {

    fprintf(stream, "%s:%3d ", file, line);
  }

  va_start(ap, level);
  fmt = va_arg(ap, const char *);
  vfprintf(stream, fmt, ap);
  fprintf(stream, "\n");
  va_end(ap);
}
