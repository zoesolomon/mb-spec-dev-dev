SUBDIRS = .

noinst_LIBRARIES = libtools.a

AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/tools \
	   $(MKL_CFLAGS)
AM_CXXFLAGS = $(CXX_FLAGS) -std=c++11

libtools_a_SOURCES = constants.h \
    orient.cpp orient.h pbc-tools.cpp pbc-tools.h \
    mat_tools.cpp mat_tools.h \
    io-tools.cpp io-tools.h  \
    read-lammpstrj.cpp read-lammpstrj.h \
    ps.h ps.cpp gammq.h gammq.cpp \
    ltp2011.h ltp2011.cpp \
    ttm4-smear.h ttm4-smear.cpp ttm4-es.h ttm4-es.cpp \
    ttm4-periodic-es.h ttm4-periodic-es.cpp ttm4-dip.h ttm4-dip.cpp \
    avila.cpp avila.h \
    poly2b5.h poly2b5.cpp \
    dip_poly2b5.h dip_poly2b5.cpp \
    lp-site.h lp-site.cpp \
    poly-l-2b3.h poly-l-2b3.cpp \
    alpha-el-poly2b3.h alpha-el-poly2b3.cpp \
    alpha-ttm4.h alpha-ttm4.cpp \
    system.h system.cpp 

