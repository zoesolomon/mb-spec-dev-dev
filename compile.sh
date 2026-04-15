#!/bin/sh


./autogen.sh
./configure CXX=mpicxx CC=mpicc --enable-maintainer-mode --disable-silent-rules PKG_CONFIG_PATH=/usr/lib64/pkgconfig CFLAGS='-g -fopenmp -O0 -Wall -std=c99' CXXFLAGS='-g -fopenmp -O0 -Wall'
./config.status
make -j4
