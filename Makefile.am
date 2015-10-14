# add path to boost.m4
ACLOCAL_AMFLAGS = -I .

# global compile and link options
AM_CPPFLAGS = $(BOOST_CPPFLAGS)

# targets to build and install
lib_LTLIBRARIES = libbaofit.la
bin_PROGRAMS = baofit

# extra targets that should not be installed
#noinst_PROGRAMS =

# targets that contain unit tests
#check_PROGRAMS = 

# add our pkgconfig file to the install target
pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = baofit.pc

# any library dependencies not already added by configure can be added here
#libbaofit_la_LIBADD = 

# instructions for building the library
libbaofit_la_SOURCES = \
	baofit/AbsCorrelationModel.cc \
	baofit/BaoCorrelationModel.cc \
	baofit/BaoKSpaceCorrelationModel.cc \
	baofit/BaoKSpaceFftCorrelationModel.cc \
	baofit/BroadbandModel.cc \
	baofit/MetalCorrelationModel.cc \
	baofit/NonLinearCorrectionModel.cc \
	baofit/XiCorrelationModel.cc \
	baofit/PkCorrelationModel.cc \
	baofit/AbsCorrelationData.cc \
	baofit/QuasarCorrelationData.cc \
	baofit/ComovingCorrelationData.cc \
	baofit/CorrelationFitter.cc \
	baofit/CorrelationAnalyzer.cc \
	baofit/RadiationModel.cc \
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
	baofit/BaoKSpaceCorrelationModel.h \
	baofit/BaoKSpaceFftCorrelationModel.h \
	baofit/BroadbandModel.h \
	baofit/MetalCorrelationModel.h \
	baofit/NonLinearCorrectionModel.h \
	baofit/XiCorrelationModel.h \
	baofit/PkCorrelationModel.h \
	baofit/AbsCorrelationData.h \
	baofit/QuasarCorrelationData.h \
	baofit/ComovingCorrelationData.h \
	baofit/CorrelationFitter.h \
	baofit/CorrelationAnalyzer.h \
	baofit/RadiationModel.h \
	baofit/boss.h

# instructions for building each program

baofit_SOURCES = src/baofit.cc
baofit_DEPENDENCIES = $(lib_LIBRARIES)
baofit_LDADD = libbaofit.la $(BOOST_PROGRAM_OPTIONS_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIBS)
