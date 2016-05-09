dnl Process this file with autoconf to produce a configure script.

AC_COPYRIGHT([
Copyright 2014, 2015 Jean-Pierre Flori, ANSSI.

This file is part of the expsum Program.

The expsum Program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published
by the Free Software Foundation; either version 3 of the License, or (at
your option) any later version.

The expsum Program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU General Public License
along with the expsum Program; see the file COPYING.  If not, see
http://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.
])

dnl Directories
AC_INIT([expsum], [1.0], [jean-pierre.flori@ssi.gouv.fr])
AC_CONFIG_HEADERS([src/config.h])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([m4])
test_CFLAGS=${CFLAGS+set}

dnl Automake
AM_INIT_AUTOMAKE()
AM_MAINTAINER_MODE([disable])

dnl Libtool
LT_INIT()

dnl Configuration

AC_ARG_WITH([pari],
   [AC_HELP_STRING([--with-pari=DIR], [PARI install directory])],
   [CPPFLAGS="$CPPFLAGS -I$withval/include"
    LDFLAGS="$LDFLAGS -L$withval/lib"],
   [AC_MSG_NOTICE([assuming that PARI is installed in a standard path])])

AC_ARG_WITH([m2],
   [AC_HELP_STRING([--with-m2=<odd integer>], [Odd extension degree])],
   [M2=$withval],
   [M2=9])

AC_ARG_WITH([modulus],
   [AC_HELP_STRING([--with-modulus=<1:minimal_weight,2:conway,other>],
   [Type of modulus to use])],
   [MODULUS=$withval],
   [MODULUS=1])

AC_ARG_WITH([threads],
   [AC_HELP_STRING([--with-threads=<integer>], [Number of threads to use])],
   [USE_THREADS=1
    NB_THREADS=$withval
    MAX_THREADS=$withval
    LIBS="-lpthread $LIBS"],
   [USE_THREADS=0
    NB_THREADS=1
    MAX_THREADS=1])

AC_ARG_WITH([assembly],
   [AS_HELP_STRING([--with-assembly=<yes,more,no>],
   [Turn on assembly code (default is yes)])],
   [],
   [with_assembly=yes])

AC_ARG_WITH([sage],
   [AC_HELP_STRING([--with-sage=PATH],[Sage script path])],
   [SAGE_CMD=$withval],
   [SAGE_CMD=sage])

AC_ARG_ENABLE([pclmul],
   [AS_HELP_STRING([--enable-pclmul],
   [Turn on pclmul code (default is yes)])],
   [],
   [enable_pclmul=yes])

AC_ARG_ENABLE([popcnt],
   [AS_HELP_STRING([--enable-popcnt],
   [Turn on popcnt code (default is yes)])],
   [],
   [enable_popcnt=yes])

AC_ARG_ENABLE([debug],
   [AS_HELP_STRING([--enable-debug],
   [Disable optimisations (default is no)])],
   [],
   [enable_debug=no])

dnl Information on host
AC_CANONICAL_HOST

dnl Exercise compiler
AC_LANG([C])
AC_PROG_CC
AC_PROG_CPP

dnl Libraries
AC_SEARCH_LIBS([pari_init], [pari], [], [
  AC_MSG_ERROR([unable to find the PARI library])
])

dnl Assembly
TUNE_ASSEMBLY([])
GENERATE_ASSEMBLY([])

dnl Preprocessor stuff
AC_DEFINE_UNQUOTED([USE_THREADS], [$USE_THREADS],
  [Define to use threads])
AC_DEFINE_UNQUOTED([NB_THREADS], [$NB_THREADS],
  [Number of threads to use])
AC_DEFINE_UNQUOTED([USE_ASSEMBLY], [$USE_ASSEMBLY],
  [Define to use assembly])
AC_DEFINE_UNQUOTED([USE_MORE_ASSEMBLY], [$USE_MORE_ASSEMBLY],
  [Define to use more assembly])

dnl Header generation
GENERATE_HEADER([])

dnl Tuning
TUNE_CFLAGS([])

dnl Output
AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
