ACLOCAL_AMFLAGS = -I build-aux/m4

lib_LTLIBRARIES = libsha2.la
sha2includedir = $(includedir)/sha2
sha2include_HEADERS = include/sha2/sha256.h
noinst_HEADERS  = lib/common.h
noinst_HEADERS += lib/compat/byteswap.h
noinst_HEADERS += lib/compat/cpuid.h
noinst_HEADERS += lib/compat/endian.h

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = libsha2.pc

libsha2_la_CPPFLAGS = -I$(srcdir)/include
libsha2_la_SOURCES  = lib/common.c
libsha2_la_SOURCES += lib/compat/byteswap.c
libsha2_la_SOURCES += lib/sha256.c
libsha2_la_SOURCES += lib/sha256_armv8.c
libsha2_la_SOURCES += lib/sha256_avx2.c
libsha2_la_SOURCES += lib/sha256_shani.c
libsha2_la_SOURCES += lib/sha256_sse4.c
libsha2_la_SOURCES += lib/sha256_sse41.c
