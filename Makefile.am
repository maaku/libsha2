ACLOCAL_AMFLAGS = -I build-aux/m4

AM_CXXFLAGS = -std=c++17

check_LTLIBRARIES = libgtest.la
libgtest_la_SOURCES = googletest/googletest/src/gtest-all.cc
libgtest_la_CPPFLAGS = -I$(top_srcdir)/googletest/googletest/include -I$(top_srcdir)/googletest/googletest
libgtest_la_LDFLAGS = -pthread

check_PROGRAMS = gtest
gtest_SOURCES = test/gtest.cc
gtest_LDADD = libgtest.la
gtest_LDFLAGS = -pthread
gtest_CPPFLAGS = -I$(top_srcdir)/googletest/googletest/include -I$(top_srcdir)/googletest/googletest -pthread

TESTS = gtest

SUBDIRS = lib
