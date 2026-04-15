SUBDIRS = .

AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/tools $(MKL_CFLAGS)

bin_PROGRAMS = response \
	       spectrum

response_SOURCES = response.cpp 
response_LDADD = ../tools/libtools.a
response_DEPS = ../tools/libtools.a

spectrum_SOURCES = spectrum.cpp 
spectrum_LDADD = 
spectrum_DEPS = 

