dnl
dnl bfam_include.m4 - custom macros
dnl

dnl Documentation for macro names: brackets indicate optional arguments

dnl bfam_ARG_ENABLE(NAME, COMMENT, TOKEN)
dnl Check for --enable/disable-NAME using shell variable BFAM_ENABLE_TOKEN
dnl If shell variable is set beforehand it overrides the option
dnl If enabled, define TOKEN to 1 and set conditional BFAM_TOKEN
dnl Default is disabled
dnl
AC_DEFUN([BFAM_ARG_ENABLE],
         [SC_ARG_ENABLE_PREFIX([$1], [$2], [$3], [BFAM])])

dnl BFAM_ARG_DISABLE(NAME, COMMENT, TOKEN)
dnl Check for --enable/disable-NAME using shell variable BFAM_ENABLE_TOKEN
dnl If shell variable is set beforehand it overrides the option
dnl If enabled, define TOKEN to 1 and set conditional BFAM_TOKEN
dnl Default is enabled
dnl
AC_DEFUN([BFAM_ARG_DISABLE],
         [SC_ARG_DISABLE_PREFIX([$1], [$2], [$3], [BFAM])])

dnl BFAM_ARG_WITH(NAME, COMMENT, TOKEN)
dnl Check for --with/without-NAME using shell variable BFAM_WITH_TOKEN
dnl If shell variable is set beforehand it overrides the option
dnl If with, define TOKEN to 1 and set conditional BFAM_TOKEN
dnl Default is without
dnl
AC_DEFUN([BFAM_ARG_WITH],
         [SC_ARG_WITH_PREFIX([$1], [$2], [$3], [BFAM])])

dnl BFAM_ARG_WITHOUT(NAME, COMMENT, TOKEN)
dnl Check for --with/without-NAME using shell variable BFAM_WITH_TOKEN
dnl If shell variable is set beforehand it overrides the option
dnl If with, define TOKEN to 1 and set conditional BFAM_TOKEN
dnl Default is with
dnl
AC_DEFUN([BFAM_ARG_WITHOUT],
         [SC_ARG_WITHOUT_PREFIX([$1], [$2], [$3], [BFAM])])

dnl BFAM_CHECK_LIBRARIES(PREFIX)
dnl This macro bundles the checks for all libraries and link tests
dnl that are required by bfam.  It can be used by other packages that
dnl link to bfam to add appropriate options to LIBS.
dnl
AC_DEFUN([BFAM_CHECK_LIBRARIES],
[
SC_REQUIRE_LIB([m], [fabs])
SC_CHECK_LIB([z], [adler32_combine], [ZLIB], [$1])
SC_CHECK_LIB([lua52 lua5.2 lua51 lua5.1 lua lua5], [lua_createtable],
	     [LUA], [$1])
])

dnl BFAM_AS_SUBPACKAGE(PREFIX)
dnl Call from a package that is using bfam as a subpackage.
dnl Sets PREFIX_DIST_DENY=yes if bfam is make install'd.
dnl
AC_DEFUN([BFAM_AS_SUBPACKAGE],
         [SC_ME_AS_SUBPACKAGE([$1], [m4_tolower([$1])], [BFAM], [bfam])])

dnl BFAM_FINAL_MESSAGES(PREFIX)
dnl This macro prints messages at the end of the configure run.
dnl
AC_DEFUN([BFAM_FINAL_MESSAGES],
[
AC_MSG_NOTICE([- $1 -------------------------------------------------
The current version of bfam's autotools support just builds a library
from selected sources.  As of yet, there are no examples or tests.
])
if test "x$$1_HAVE_LUA" = x ; then
AC_MSG_ERROR([- $1 -------------------------------------------------
We did not find a recent lua containing the function lua_createtable.
You can fix this by installing the lua-dev package for your software
environment or by compiling a working lua and pointing LIBS to it.
])
fi
])
