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
	baofit/AbsCorrelationModel.cc \
	baofit/BaoCorrelationModel.cc \
	baofit/BroadbandModel.cc \
	baofit/XiCorrelationModel.cc \
	baofit/PkCorrelationModel.cc \
	baofit/AbsCorrelationData.cc \
	baofit/QuasarCorrelationData.cc \
	baofit/ComovingCorrelationData.cc \
	baofit/MultipoleCorrelationData.cc \
	baofit/CorrelationFitter.cc \
	baofit/CorrelationAnalyzer.cc \
	baofit/boss.cc

# library headers to install (nobase prefix preserves any subdirectories)
# Anything that includes config.h should *not* be listed here.
# and is not part of the package public API.
nobase_include_HEADERS = \
	baofit/baofit.h \
	baofit/types.h \
	baofit/RuntimeError.h \
	baofit/AbsCorrelationModel.h \
	baofit/BaoCorrelationModel.h \
	baofit/BroadbandModel.h \
	baofit/XiCorrelationModel.h \
	baofit/PkCorrelationModel.h \
	baofit/AbsCorrelationData.h \
	baofit/QuasarCorrelationData.h \
	baofit/ComovingCorrelationData.h \
	baofit/MultipoleCorrelationData.h \
	baofit/CorrelationFitter.h \
	baofit/CorrelationAnalyzer.h \
	baofit/boss.h

# instructions for building each program

baofit_SOURCES = src/baofit.cc
baofit_DEPENDENCIES = $(lib_LIBRARIES)
baofit_LDADD = -lboost_program_options -L. -lbaofit -lcosmo -lMinuit2 -lblas
