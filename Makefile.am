## Makefile.am in the toplevel directory
##
## This file is ignored when building with cmake.
## It is used only when building with the autotools.

ACLOCAL_AMFLAGS = -I autotools

DISTCLEANFILES =
EXTRA_DIST =
lib_LTLIBRARIES =
include_HEADERS =
nodist_include_HEADERS =

EXTRA_DIST += \
	autotools/bootstrap \
	autotools/git-version-gen autotools/git2cl

include src/Makefile.am
