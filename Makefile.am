# global compile and link options
AM_CPPFLAGS = $(BOOST_CPPFLAGS)
AM_LDFLAGS = -llikely

# targets to build and install
lib_LTLIBRARIES = libbaofit.la
bin_PROGRAMS = baofit

# extra targets that should not be installed
## noinst_PROGRAMS =

# instructions for building the library
libbaofit_la_SOURCES = \
	baofit/LyaBaoModel.cc

# library headers to install (nobase prefix preserves any subdirectories)
# Anything that includes config.h should *not* be listed here.
# and is not part of the package public API.
nobase_include_HEADERS = \
	baofit/baofit.h \
	baofit/types.h \
	baofit/RuntimeError.h \
	baofit/LyaBaoModel.h

# instructions for building each program

baofit_SOURCES = src/baofit.cc
baofit_DEPENDENCIES = $(lib_LIBRARIES)
baofit_LDADD = -lboost_program_options -L. -lcosmo -lMinuit2 -lblas
