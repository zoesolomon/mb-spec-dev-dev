dnl (c) countless GNU projects  

AC_PREREQ(2.54)

dnl define the version of this package
m4_define([version_major], [0])
m4_define([version_minor], [0])
m4_define([version_micro], [0])

m4_define([the_version],[version_major.version_minor.version_micro])

dnl initialize autoconf
AC_INIT([se-param], [the_version], [gmedders@ucsd.edu])
AM_INIT_AUTOMAKE([foreign 1.9 no-define dist-bzip2 no-dist-gzip])

AC_CONFIG_SRCDIR([src/tools/io-tools.cpp])
AC_CONFIG_HEADERS(config.h)

dnl support silent build rules
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

dnl need this one
AC_PROG_CC
AC_PROG_CXX
AC_LANG([C++])

AM_MAINTAINER_MODE

AC_ARG_ENABLE(ansi,
    AC_HELP_STRING([--enable-ansi=@<:@no/yes@:>@],
                   [pass -ansi -pedantic to gcc @<:@default=no@:>@]),,
                   enable_ansi=no)

dnl build time sanity check
AM_SANITY_CHECK

dnl checks for programs
AC_PROG_RANLIB
AC_PROG_INSTALL
AC_PROG_MAKE_SET

changequote(,)dnl
if test "x$GXX" = "xyes"; then
    case " $CXXFLAGS " in
        *[\ \	]-Wall[\ \	]*) ;;
        *) CXXFLAGS="$CXXFLAGS -Wall" ;;
    esac

    if test "x$enable_ansi" = "xyes"; then
        case " $CXXFLAGS " in
            *[\ \	]-ansi[\ \	]*) ;;
            *) CXXFLAGS="$CXXFLAGS -ansi" ;;
        esac

        case " $CXXFLAGS " in
            *[\ \	]-pedantic[\ \	]*) ;;
            *) CXXFLAGS="$CXXFLAGS -pedantic" ;;
        esac
    fi
fi
changequote([,])dnl

dnl GSL
saved_LIBS="$LIBS"
AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm],[],
[AC_MSG_ERROR([
************************************************************************
**  cannot continue without GSL (http://www.google.com/search?q=GSL)  **
************************************************************************])])
LIBS="$saved_LIBS"

GSL_LIBS="-lgsl -lgslcblas"
AC_SUBST(GSL_LIBS)

dnl MKL
AC_ARG_ENABLE([mkl],
    AC_HELP_STRING([--enable-mkl=@<:@yes/no@:>@],
                   [use MKL for fitting @<:@default=yes@:>@]),,
                   enable_mkl=yes)

AC_ARG_WITH([mkl-cflags],
    AC_HELP_STRING([--with-mkl-cflags=FLAGS],
                   [MKL compile flags]),
                   [MKL_CFLAGS=$with_mkl_cflags],
                   [MKL_CFLAGS='-m64 -I/opt/intel/mkl/include'])

AC_ARG_WITH([mkl-libs],
    AC_HELP_STRING([--with-mkl-libs=FLAGS],
                   [MKL linking flags]),
                   [MKL_LIBS=$with_mkl_libs],
                   [MKL_LIBS='-L/opt/intel/lib/intel64 -L/opt/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread'])

if test "x$enable_mkl" = "xyes"; then
saved_LIBS=$LIBS
saved_CXXFLAGS=$CXXFLAGS
LIBS="$LIBS $MKL_LIBS"
CXXFLAGS="$CXXFLAGS $MKL_CFLAGS"
AC_MSG_CHECKING([MKL usability])
    AC_TRY_LINK([#include <mkl_lapack.h>],
                [dgelsd(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);],
                [AC_MSG_RESULT([ok])],
                [AC_MSG_RESULT([not really]); enable_mkl=no])
CXXFLAGS=$saved_CXXFLAGS
LIBS=$saved_LIBS
fi

if test "x$enable_mkl" != "xyes"; then
    MKL_CFLAGS=''
    MKL_LIBS=''
fi

if test "x$enable_mkl" = "xyes"; then
    AC_DEFINE([USE_MKL], [1], [use MKL for fitting])
    AC_SUBST(MKL_CFLAGS)
    AC_SUBST(MKL_LIBS)
fi

dnl fftw3 needed for pimd
dnl PKG_CHECK_MODULES([FFTW], [fftw3])
dnl

AC_CONFIG_FILES([
  Makefile
  src/Makefile
  src/tools/Makefile
  src/infrared/Makefile
  src/raman/Makefile
  src/sfg/Makefile
  src/power/Makefile
])

AC_OUTPUT

dnl
dnl summary
dnl

echo "<=============================================================>
 Configuration summary :
<------------------------------------------------------------->

       CC  = ${CC}
    CFLAGS = ${CFLAGS}
       CXX = ${CXX}
  CXXFLAGS = ${CXXFLAGS}
  CPPFLAGS = ${CPPFLAGS}
   LDFLAGS = ${LDFLAGS}
  IFC_LIBS = ${IFC_LIBS}

    PREFIX = ${prefix}

<=============================================================>"
