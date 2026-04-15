#ifndef DIP_POLY2B5_H
#define DIP_POLY2B5_H

#include <cstddef>
#include <iostream>

namespace h2o {

struct dip_poly2b5 {

    dip_poly2b5();

    static const size_t ncoeffs_H = 1826;
    static const size_t ncoeffs_O = 1052;
    static const size_t ncoeffs = ncoeffs_H + ncoeffs_O;

    static void d(const double* w1, const double* w2, 
		  double* poly_x, double* poly_y, double* poly_z);

    static void eval(const double* w1, const double* w2, double* dip);

};

} // namespace h2o

#endif // DIP_POLY2B5_H
