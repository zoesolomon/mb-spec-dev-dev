SUBDIRS = .

AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/tools $(MKL_CFLAGS)
AM_CXXFLAGS = $(CXX_FLAGS) -std=c++11

bin_PROGRAMS = polarizability \
	       response_isotropic response_depolarized \
	       spectrum

polarizability_SOURCES = polarizability.cpp 
polarizability_LDADD = ../tools/libtools.a
polarizability_DEPS = ../tools/libtools.a

response_isotropic_SOURCES = response_isotropic.cpp 
response_isotropic_LDADD = ../tools/libtools.a
response_isotropic_DEPS = ../tools/libtools.a

response_depolarized_SOURCES = response_depolarized.cpp
response_depolarized_LDADD = ../tools/libtools.a
response_depolarized_DEPS = ../tools/libtools.a

spectrum_SOURCES = spectrum.cpp 
spectrum_LDADD = 
spectrum_DEPS = 

