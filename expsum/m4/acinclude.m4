AC_DEFUN([GENERATE_HEADER], [
AC_CONFIG_COMMANDS([src/expsum.h], [
AS_IF([
$ECHO "$MAX_THREADS
$MODULUS
$M2
" | "$SAGE_CMD" ./scripts/expsum.py > /dev/null], [
], [
AC_MSG_ERROR([failed to generate header!])
])
], [
  MAX_THREADS=$MAX_THREADS
  MODULUS=$MODULUS
  M2=$M2
  SAGE_CMD="$SAGE_CMD"
])
  ])

AC_DEFUN([TUNE_ASSEMBLY], [
AS_CASE([$with_assembly],
   [more], [
    USE_ASSEMBLY=1
    USE_MORE_ASSEMBLY=1
   ],
   [no], [
    USE_ASSEMBLY=0
    USE_MORE_ASSEMBLY=0
   ],
   [
    USE_ASSEMBLY=1
    USE_MORE_ASSEMBLY=0
 ])
])

dnl currently handled by the Sage script
AC_DEFUN([GENERATE_ASSEMBLY], [
])

AC_DEFUN([TUNE_CFLAGS], [
gcc_cflags="-mpclmul -msse4.2 -Wall -pedantic"
gcc_cflags_nodebug="-O3 -funroll-loops -fomit-frame-pointer"
gcc_cflags_debug="-g -O0 -fno-inline-functions"
dnl If compiler is gcc, then use some specific flags.
dnl But don't touch user other flags.
AS_IF([test "$test_CFLAGS" != set && test "$GCC" = yes],
  [AS_IF([test "$enable_debug" = yes], [
    CFLAGS="$gcc_cflags_debug $gcc_cflags $test_CFLAGS"
   ], [
    CFLAGS="$gcc_cflags_nodebug $gcc_cflags $test_CFLAGS"
   ])
])
AC_MSG_NOTICE([using CFLAGS="$CFLAGS"])
])
