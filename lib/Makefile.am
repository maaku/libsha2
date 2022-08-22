lib_LTLIBRARIES = libsha2.la
sha2includedir = $(includedir)/sha2
sha2include_HEADERS = $(top_srcdir)/include/sha2/sha256.h
noinst_HEADERS  = common.h
noinst_HEADERS += compat/byteswap.h
noinst_HEADERS += compat/cpuid.h
noinst_HEADERS += compat/endian.h

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = libsha2.pc

libsha2_la_CPPFLAGS = -I$(top_srcdir)/include
libsha2_la_CFLAGS  = $(SSE41_CFLAGS)
libsha2_la_CFLAGS += $(AVX2_CFLAGS)
libsha2_la_CFLAGS += $(X86_SHANI_CFLAGS)
libsha2_la_CFLAGS += $(ARM_SHANI_CFLAGS)
libsha2_la_SOURCES  = common.c
libsha2_la_SOURCES += compat/byteswap.c
libsha2_la_SOURCES += sha256.c
libsha2_la_SOURCES += sha256_armv8.c
libsha2_la_SOURCES += sha256_avx2.c
libsha2_la_SOURCES += sha256_shani.c
libsha2_la_SOURCES += sha256_sse4.c
libsha2_la_SOURCES += sha256_sse41.c
