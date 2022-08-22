AM_CXXFLAGS = -std=c++17

check_LTLIBRARIES = libgtest.la
libgtest_la_SOURCES = $(top_srcdir)/googletest/googletest/src/gtest-all.cc
libgtest_la_CPPFLAGS = -I$(top_srcdir)/googletest/googletest/include -I$(top_srcdir)/googletest/googletest
libgtest_la_LDFLAGS = -pthread

check_PROGRAMS = sha2
sha2_SOURCES = sha2.cc
sha2_LDADD = libgtest.la
sha2_LDFLAGS = -pthread
sha2_CPPFLAGS = -I$(top_srcdir)/googletest/googletest/include -I$(top_srcdir)/googletest/googletest -pthread

TESTS = sha2
