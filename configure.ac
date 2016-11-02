AC_PREREQ(2.62)
AC_INIT([relectro],[0.01])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/dat_files.c])
AC_PROG_CC
AC_PROG_CXX
AM_PROG_CC_STDC
AM_PROG_CC_C_O
AC_SEARCH_LIBS([pthread_create],[pthread])
AC_SEARCH_LIBS([clock_gettime],[rt])
AC_SEARCH_LIBS([sqrt],[m])
PKG_CHECK_MODULES([FFTW3], [fftw3])
AC_FUNC_MALLOC
AC_CHECK_FUNCS([clock_gettime])
AC_CHECK_FUNCS([memset])
AC_CHECK_HEADERS([fcntl.h])
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
