AM_CXXFLAGS = $(CXX_FLAGS) -std=c++11

SUBDIRS = .

#noinst_LIBRARIES = ../tools/libtools.a

AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/tools $(MKL_CFLAGS)

bin_PROGRAMS = dipole \
	       response \
               spectrum

dipole_SOURCES = dipole.cpp 
dipole_LDADD = ../tools/libtools.a $(MKL_LIBS)
dipole_DEPS = ../tools/libtools.a

response_SOURCES = response.cpp 
response_LDADD = ../tools/libtools.a $(MKL_LIBS)
response_DEPS = ../tools/libtools.a

spectrum_SOURCES = spectrum.cpp 
spectrum_LDADD = ../tools/libtools.a $(MKL_LIBS)
spectrum_DEPS = ../tools/libtools.a

