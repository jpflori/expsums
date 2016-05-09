AC_DEFUN([GENERATE_HEADER],[
echo "$MODULUS\n$M2\n$" | $SAGE ./scripts/expsum.py
  ])


AC_DEFUN([TUNE_CFLAGS], [
gcc_cflags="-Wall -pedantic"
gcc_cflags_nodebug=-march=native -O3 -funroll-loops -fomit-frame-pointer
gcc_cflags_debug=-g -O0 -fno-inline-functions
dnl If compiler is gcc, then use some specific flags.
dnl But don't touch user other flags.
if test "$test_CFLAGS" != set && test -n "$GCC"; then
  if test "$enable_debug" = yes; then
    CFLAGS="$gcc_cflags_debug $gcc_cflags $CFLAGS"
  else
    CFLAGS="$gcc_cflags_nodebug $gcc_cflags $CFLAGS"
  fi
fi
  ])


AC_DEFUN([PCLMUL_EXAMPLE],[AC_LANG_SOURCE([
#include <wmmintrin.h>
int main() {
return _mm_cvtsi128_si64(_mm_clmulepi64_si128(_mm_cvtsi64_si128(17), _mm_cvtsi64_si128(42), 0));
}
])])

# Check whether we need some flag such as -mpclmul in order to enable pclmulqdq
# support
AC_DEFUN([CHECK_PCLMUL_SUPPORT],[
 ac_save_CFLAGS="$CFLAGS"
 AC_CACHE_CHECK([whether $CC can compile pclmulqdq and if it is supported by the hardware], [gf2x_cv_cc_supports_pclmul],[
  gf2x_cv_cc_supports_pclmul=no
  if test "x${enable_pclmul}" = xno ; then
   echo $ECHO_N " disabled, "
  else
   AC_RUN_IFELSE([PCLMUL_EXAMPLE()],[
    gf2x_cv_cc_supports_pclmul=yes
   ],[
    CFLAGS="$ac_save_CFLAGS -mpclmul"
    AC_RUN_IFELSE([PCLMUL_EXAMPLE()],[
     gf2x_cv_cc_supports_pclmul="requires -mpclmul"
    ],[
     gf2x_cv_cc_supports_pclmul=no
    ])
   ],[
   echo $ECHO_N " cross-compiling, "
   gf2x_cv_cc_supports_pclmul=no
   ])
  fi
 ])
 ac_save_CPPFLAGS=$CPPFLAGS
 if test "$gf2x_cv_cc_supports_pclmul" = "requires -mpclmul" ;then
  # Tweaking CFLAGS is often not enough.
  AC_CACHE_CHECK([whether -mpclmul is also needed by the preprocessor],
   [gf2x_cv_cpp_requires_mpclmul_flag],[
   AC_PREPROC_IFELSE([PCLMUL_EXAMPLE()],[
    gf2x_cv_cpp_requires_mpclmul_flag=no
   ],[
    CPPFLAGS="$ac_save_CPPFLAGS -mpclmul"
    AC_PREPROC_IFELSE([PCLMUL_EXAMPLE()],[
    gf2x_cv_cpp_requires_mpclmul_flag=yes
    ],[
     AC_MSG_ERROR([Sorry, the preprocessor can't parse pclmul !])
    ])
   ])
  ])
 fi
 CFLAGS="$ac_save_CFLAGS"
 CPPFLAGS="$ac_save_CPPFLAGS"
 if test "$gf2x_cv_cc_supports_pclmul" = "requires -mpclmul" ;then
  CFLAGS="$CFLAGS -mpclmul"
 fi
 if test "$gf2x_cv_cpp_requires_mpclmul_flag" = "yes" ; then
  CPPFLAGS="$CPPFLAGS -mpclmul"
 fi
 if test "$gf2x_cv_cc_supports_pclmul" != "no" ;then
  AC_DEFINE([GF2X_HAVE_PCLMUL_SUPPORT],[1],[Define if pclmul as present in the source tree is supported by the compiler and hardware])
 fi
])# CHECK_PCLMUL_SUPPORT


AC_DEFUN([POPCNT_EXAMPLE],[AC_LANG_SOURCE([
#include <wmmintrin.h>
#include <assert.h>
int main() {
assert(sizeof(unsigned long) == 8); /* assume 64-bit */
#if defined(__GNUC__) && __GNUC__ == 4 &&__GNUC_MINOR__ == 1
#define _gf2x_mm_cvtsi64_m64(u) _mm_cvtsi64x_m64((u))
#else
#define _gf2x_mm_cvtsi64_m64(u) _mm_cvtsi64_m64((u))
#endif
/* _m128i from 2 int64_t's */
#define _gf2x_mm_setr_epi64(lo, hi)                      		\
    _mm_setr_epi64(                                      		\
            _gf2x_mm_cvtsi64_m64((int64_t) (lo)),       		\
            _gf2x_mm_cvtsi64_m64((int64_t) (hi))        		\
        )
/* _m128i from 1 int64_t's */
#define _gf2x_mm_set1_epi64(u) _mm_set1_epi64( _gf2x_mm_cvtsi64_m64((int64_t) (u)))
__m128i a = _gf2x_mm_set1_epi64(17);
__m128i b = _gf2x_mm_set1_epi64(42);
union { __m128i s; unsigned long x[[2]]; } proxy;
proxy.s = _mm_clmulepi64_si128(a, b, 0);
return proxy.x[[0]] - 650;
}
])])

# Check whether we need some flag such as -mpopcnt in order to enable popcnt
# support
AC_DEFUN([CHECK_POPCNT_SUPPORT],[
 ac_save_CFLAGS="$CFLAGS"
 AC_CACHE_CHECK([whether $CC can compile popcnt and if it is supported by the hardware], [gf2x_cv_cc_supports_popcnt],[
  gf2x_cv_cc_supports_popcnt=no
  if test "x${enable_popcnt}" = xno ; then
   echo $ECHO_N " disabled, "
  else
   AC_RUN_IFELSE([POPCNT_EXAMPLE()],[
    gf2x_cv_cc_supports_popcnt=yes
   ],[
    CFLAGS="$ac_save_CFLAGS -mpopcnt"
    AC_RUN_IFELSE([POPCNT_EXAMPLE()],[
     gf2x_cv_cc_supports_popcnt="requires -mpopcnt"
    ],[
     gf2x_cv_cc_supports_popcnt=no
    ])
   ],[
   echo $ECHO_N " cross-compiling, "
   gf2x_cv_cc_supports_popcnt=no
   ])
  fi
 ])
 ac_save_CPPFLAGS=$CPPFLAGS
 if test "$gf2x_cv_cc_supports_popcnt" = "requires -mpopcnt" ;then
  # Tweaking CFLAGS is often not enough.
  AC_CACHE_CHECK([whether -mpopcnt is also needed by the preprocessor],
   [gf2x_cv_cpp_requires_mpopcnt_flag],[
   AC_PREPROC_IFELSE([POPCNT_EXAMPLE()],[
    gf2x_cv_cpp_requires_mpopcnt_flag=no
   ],[
    CPPFLAGS="$ac_save_CPPFLAGS -mpopcnt"
    AC_PREPROC_IFELSE([POPCNT_EXAMPLE()],[
    gf2x_cv_cpp_requires_mpopcnt_flag=yes
    ],[
     AC_MSG_ERROR([Sorry, the preprocessor can't parse popcnt !])
    ])
   ])
  ])
 fi
 CFLAGS="$ac_save_CFLAGS"
 CPPFLAGS="$ac_save_CPPFLAGS"
 if test "$gf2x_cv_cc_supports_popcnt" = "requires -mpopcnt" ;then
  CFLAGS="$CFLAGS -mpopcnt"
 fi
 if test "$gf2x_cv_cpp_requires_mpopcnt_flag" = "yes" ; then
  CPPFLAGS="$CPPFLAGS -mpopcnt"
 fi
 if test "$gf2x_cv_cc_supports_popcnt" != "no" ;then
  AC_DEFINE([GF2X_HAVE_POPCNT_SUPPORT],[1],[Define if popcnt as present in the source tree is supported by the compiler and hardware])
 fi
])# CHECK_POPCNT_SUPPORT
