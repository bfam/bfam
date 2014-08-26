dnl
dnl This file is part of the bfam code.
dnl

AC_INIT([bfam],
        [m4_esyscmd([autotools/git-version-gen .tarball-version])],
        [bfam@librelist.com])
AC_PREREQ(2.61)
AC_CONFIG_HEADERS([src/pre_config.h])
AC_CONFIG_SRCDIR([src/bfam.h])
AC_CONFIG_AUX_DIR([autotools])
AC_CONFIG_MACRO_DIR([autotools])
AC_PREFIX_DEFAULT([$PWD/local])
AX_PREFIX_CONFIG_H([src/bfam_config.h], [BFAM])
AM_INIT_AUTOMAKE([parallel-tests subdir-objects])
AM_SILENT_RULES
SC_VERSION([BFAM])

echo "o---------------------------------------"
echo "| Checking options"
echo "o---------------------------------------"

AC_DEFINE([AUTOTOOLS], 1, [bfam is configured by autotools])
BFAM_ARG_ENABLE([debug], [enable debug mode (assertions and extra checks)],
                [DEBUG])

echo "o---------------------------------------"
echo "| Checking MPI and related programs"
echo "o---------------------------------------"

dnl A nonempty second/third argument causes to enable F77+FC/CXX, respectively.
SC_MPI_CONFIG([BFAM])
SC_MPI_ENGAGE([BFAM])
# This is needed for compatibility with automake >= 1.12
m4_ifdef([AM_PROG_AR],[AM_PROG_AR])
LT_INIT

echo "o---------------------------------------"
echo "| Checking libraries"
echo "o---------------------------------------"

BFAM_CHECK_LIBRARIES([BFAM])

echo "o---------------------------------------"
echo "| Checking headers"
echo "o---------------------------------------"

AC_CHECK_HEADERS([execinfo.h signal.h sys/time.h sys/types.h time.h])
AC_CHECK_HEADERS([lua.h lua5.1/lua.h lua5.2/lua.h])

echo "o---------------------------------------"
echo "| Checking functions"
echo "o---------------------------------------"

dnl AC_CHECK_FUNCS([backtrace backtrace_symbols strtol strtoll])

echo "o---------------------------------------"
echo "| Checking keywords and types"
echo "o---------------------------------------"

AC_C_BIGENDIAN
AC_C_CONST
AC_C_INLINE
AC_C_RESTRICT
AC_CHECK_SIZEOF([long])
AC_CHECK_SIZEOF([long long])
AC_CHECK_SIZEOF([unsigned long])
AC_CHECK_SIZEOF([unsigned long long])
AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T

# Print summary.

dnl AC_DEFINE_UNQUOTED(F77,         ["${F77}"],         [F77 compiler])
dnl AC_DEFINE_UNQUOTED(FFLAGS,      ["${FFLAGS}"],      [F77 compiler flags])
dnl AC_DEFINE_UNQUOTED(FC,          ["${FC}"],          [FC compiler])
dnl AC_DEFINE_UNQUOTED(FCFLAGS,     ["${FCFLAGS}"],     [FC compiler flags])
AC_DEFINE_UNQUOTED(CPP,         ["${CPP}"],         [C preprocessor])
AC_DEFINE_UNQUOTED(CPPFLAGS,    ["${CPPFLAGS}"],    [C preprocessor flags])
AC_DEFINE_UNQUOTED(CC,          ["${CC}"],          [C compiler])
AC_DEFINE_UNQUOTED(CFLAGS,      ["${CFLAGS}"],      [C compiler flags])
dnl AC_DEFINE_UNQUOTED(CXX,         ["${CXX}"],         [CXX compiler])
dnl AC_DEFINE_UNQUOTED(CXXFLAGS,    ["${CXXFLAGS}"],    [CXX compiler flags])
AC_DEFINE_UNQUOTED(LDFLAGS,     ["${LDFLAGS}"],     [Linker flags])
AC_DEFINE_UNQUOTED(LIBS,        ["${LIBS}"],        [Libraries])
dnl AC_DEFINE_UNQUOTED(BLAS_LIBS,   ["${BLAS_LIBS}"],   [BLAS libraries])
dnl AC_DEFINE_UNQUOTED(LAPACK_LIBS, ["${LAPACK_LIBS}"], [LAPACK libraries])
dnl AC_DEFINE_UNQUOTED(FLIBS,       ["${FLIBS}"],       [Fortran libraries])

echo "o----------------------------------"
echo "| Results for bfam are"
echo "o----------------------------------"
dnl echo "| F77:         $F77"
dnl echo "| FFLAGS:      $FFLAGS"
dnl echo "| FC:          $FC"
dnl echo "| FCFLAGS:     $FCFLAGS"
echo "| CPP:         $CPP"
echo "| CPPFLAGS:    $CPPFLAGS"
echo "| CC:          $CC"
echo "| CFLAGS:      $CFLAGS"
dnl echo "| CXX:         $CXX"
dnl echo "| CXXFLAGS:    $CXXFLAGS"
echo "| LDFLAGS:     $LDFLAGS"
echo "| LIBS:        $LIBS"
dnl echo "| FLIBS:       $FLIBS"
dnl echo "| BLAS_LIBS:   $BLAS_LIBS"
dnl echo "| LAPACK_LIBS: $LAPACK_LIBS"
echo "o----------------------------------"

# Create output files.
dnl AC_CONFIG_FILES([Makefile Makefile.bfam.pre Doxyfile])
AC_CONFIG_FILES([Makefile])
AC_OUTPUT

# Final notifications.
BFAM_FINAL_MESSAGES([BFAM])
