## Makefile.am in the src subdirectory
##
## This file is ignored when building with cmake.
## It is used only when building with the autotools.
##
## It is included by the toplevel Makefile.am.
## We do not use recursive make.

libbfam_generated_files = src/bfam_config.h
libbfam_source_files = \
src/bfam_jacobi.c \
src/bfam_util.c
libbfam_include_files = \
src/bfam_base.h \
src/bfam_clock.h \
src/bfam_communicator.h \
src/bfam_critbit.h \
src/bfam_dictionary.h \
src/bfam_domain.h \
src/bfam_domain_pxest_2.h \
src/bfam_domain_pxest_3.h \
src/bfam_gopt.h \
src/bfam.h \
src/bfam_jacobi.h \
src/bfam_kron.h \
src/bfam_log.h \
src/bfam_opencl.h \
src/bfam_subdomain_dgx.h \
src/bfam_subdomain.h \
src/bfam_timestep_adams.h \
src/bfam_timestep.h \
src/bfam_timestep_lsrk.h \
src/bfam_util.h \
src/bfam_version.h \
src/bfam_vtk.h

DISTCLEANFILES += $(libbfam_generated_files)
include_HEADERS += $(libbfam_include_files)
nodist_include_HEADERS += $(libbfam_generated_files)

lib_LTLIBRARIES += libbfam.la
libbfam_la_SOURCES = $(libbfam_source_files) $(libbfam_include_files)
libbfam_la_CPPFLAGS = -I @top_srcdir@/src
