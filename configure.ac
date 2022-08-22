AC_PREREQ([2.60])
AC_INIT([libsha2],[0.1])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([build-aux/m4])

AC_CANONICAL_HOST

AH_TOP([#ifndef LIBSHA2__CONFIG_H])
AH_TOP([#define LIBSHA2__CONFIG_H 1])
AH_BOTTOM([#endif /*LIBSHA2__CONFIG_H*/])

AM_INIT_AUTOMAKE([foreign subdir-objects])
AM_PROG_AR
LT_INIT

dnl make the compilation flags quiet unless V=1 is used
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

PKG_PROG_PKG_CONFIG

AC_PATH_TOOL(AR, ar)
AC_PATH_TOOL(RANLIB, ranlib)
AC_PATH_TOOL(STRIP, strip)

if test "x$CFLAGS" = "x"; then
    CFLAGS="-O3 -g"
fi

AC_PROG_CC
if test x"$ac_cv_prog_cc_c89" = x"no"; then
    AC_MSG_ERROR([c89 compiler support required])
fi
AM_PROG_CC_C_O

AC_PROG_CXX

CFLAGS="$CFLAGS -W"

warn_CFLAGS="-std=c89 -pedantic -Wall -Wextra -Wcast-align -Wnested-externs -Wshadow -Wstrict-prototypes -Wno-unused-function -Wno-long-long -Wno-overlength-strings"
saved_CFLAGS="$CFLAGS"
CFLAGS="$CFLAGS $warn_CFLAGS"
AC_MSG_CHECKING([if ${CC} supports ${warn_CFLAGS}])
AC_COMPILE_IFELSE([AC_LANG_SOURCE([[char foo;]])],
    [ AC_MSG_RESULT([yes]) ],
    [ AC_MSG_RESULT([no])
      CFLAGS="$saved_CFLAGS"
    ])

saved_CFLAGS="$CFLAGS"
CFLAGS="$CFLAGS -fvisibility=hidden"
AC_MSG_CHECKING([if ${CC} supports -fvisibility=hidden])
AC_COMPILE_IFELSE([AC_LANG_SOURCE([[char foo;]])],
    [ AC_MSG_RESULT([yes]) ],
    [ AC_MSG_RESULT([no])
      CFLAGS="$saved_CFLAGS"
    ])

AC_C_INLINE
AC_C_RESTRICT

AX_CHECK_COMPILE_FLAG([-Werror], [CFLAG_WERROR="-Werror"], [CFLAG_WERROR=""])

dnl x86_64

AX_CHECK_COMPILE_FLAG([-msse4.1], [SSE41_CFLAGS="-msse4.1"], [], [$CFLAG_WERROR])
TEMP_CFLAGS="$CFLAGS"
CFLAGS="$CFLAGS $SSE41_CFLAGS"
AC_MSG_CHECKING([for SSE4.1 intrinsics])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
    #include <stdint.h>
    #include <immintrin.h>
  ]],[[
    __m128i l = _mm_set1_epi32(0);
    return _mm_extract_epi32(l, 3);
  ]])],
 [ AC_MSG_RESULT([yes]); enable_sse41=yes; AC_DEFINE([ENABLE_SSE41], [1], [Define this symbol to
 build code that uses SSE4.1 intrinsics]) ],
 [ AC_MSG_RESULT([no])]
)
CFLAGS="$TEMP_CFLAGS"
AC_SUBST(SSE41_CFLAGS)

AX_CHECK_COMPILE_FLAG([-mavx -mavx2], [AVX2_CFLAGS="-mavx -mavx2"], [], [$CFLAG_WERROR])
TEMP_CFLAGS="$CFLAGS"
CFLAGS="$CFLAGS $AVX2_CFLAGS"
AC_MSG_CHECKING([for AVX2 intrinsics])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
    #include <stdint.h>
    #include <immintrin.h>
  ]],[[
    __m256i l = _mm256_set1_epi32(0);
    return _mm256_extract_epi32(l, 7);
  ]])],
 [ AC_MSG_RESULT([yes]); enable_avx2=yes; AC_DEFINE([ENABLE_AVX2], [1], [Define this symbol to build code that uses AVX2 intrinsics]) ],
 [ AC_MSG_RESULT([no])]
)
CFLAGS="$TEMP_CFLAGS"
AC_SUBST(AVX2_CFLAGS)

AX_CHECK_COMPILE_FLAG([-msse4 -msha], [X86_SHANI_CFLAGS="-msse4 -msha"], [], [$CFLAG_WERROR])
TEMP_CFLAGS="$CFLAGS"
CFLAGS="$CFLAGS $X86_SHANI_CFLAGS"
AC_MSG_CHECKING([for SHA-NI intrinsics])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
    #include <stdint.h>
    #include <immintrin.h>
  ]],[[
    __m128i i = _mm_set1_epi32(0);
    __m128i j = _mm_set1_epi32(1);
    __m128i k = _mm_set1_epi32(2);
    return _mm_extract_epi32(_mm_sha256rnds2_epu32(i, i, k), 0);
  ]])],
 [ AC_MSG_RESULT([yes]); enable_shani=yes; AC_DEFINE([ENABLE_SHANI], [1], [Define this symbol to build code that uses SHA-NI intrinsics]) ],
 [ AC_MSG_RESULT([no])]
)
CFLAGS="$TEMP_CFLAGS"
AC_SUBST(X86_SHANI_CFLAGS)

dnl ARMv8

AX_CHECK_COMPILE_FLAG([-march=armv8-a+crc+crypto], [ARM_SHANI_CFLAGS="-march=armv8-a+crc+crypto"], [], [$CFLAG_WERROR])
TEMP_CFLAGS="$CFLAGS"
CFLAGS="$CFLAGS $ARM_SHANI_CFLAGS"
AC_MSG_CHECKING([for ARMv8 SHA-NI intrinsics])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
    #include <arm_acle.h>
    #include <arm_neon.h>
  ]],[[
    uint32x4_t a, b, c;
    vsha256h2q_u32(a, b, c);
    vsha256hq_u32(a, b, c);
    vsha256su0q_u32(a, b);
    vsha256su1q_u32(a, b, c);
  ]])],
 [ AC_MSG_RESULT([yes]); enable_arm_shani=yes; AC_DEFINE([ENABLE_ARM_SHANI], [1], [Define this symbol to build code that uses ARMv8 SHA-NI intrinsics]) ],
 [ AC_MSG_RESULT([no])]
)
CFLAGS="$TEMP_CFLAGS"
AC_SUBST(ARM_SHANI_CFLAGS)

AC_CONFIG_HEADERS([lib/config/libsha2-config.h])
AC_CONFIG_FILES([Makefile lib/Makefile lib/libsha2.pc])

dnl make sure nothing new is exported so that we don't break the cache
PKGCONFIG_PATH_TEMP="$PKG_CONFIG_PATH"
unset PKG_CONFIG_PATH
PKG_CONFIG_PATH="$PKGCONFIG_PATH_TEMP"

AC_OUTPUT
