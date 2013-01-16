#ifndef BFAM_H
#define BFAM_H

#include <bfam_config.h>

#include <stdio.h>
#include <stdlib.h>

#if defined (BFAM_HAVE_SYSEXITS_H)
#include <sysexits.h>
#elif defined (BFAM_HAVE_SYS_SYSEXITS_H)
#include <sys/sysexits.h>
#else
#define EX_OSERR EXIT_FAILURE
#define EX_USAGE EXIT_FAILURE
#endif

#include <mpi.h>

#include <bfam_version.h>

#endif
