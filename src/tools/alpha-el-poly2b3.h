#ifndef ALP_EL_POLY2B3_H
#define ALP_EL_POLY2B3_H

#include <cstddef>
#include <iostream>

namespace h2o {

struct alpha_el_poly2b3 {

    alpha_el_poly2b3();

    static const size_t ncoeffs_O = 435;
    static const size_t ncoeffs_H = 603;
    static const size_t ncoeffs_L = 588;
    static const size_t ncoeffs = ncoeffs_O + ncoeffs_H + ncoeffs_L;

    static void d(const double* inp_w1, const double* inp_w2, 
	          double* poly_xx, double* poly_yy, double* poly_zz,
		  double* poly_xy, double* poly_xz, double* poly_yz);

    static void eval(const double* inp_w1, const double* inp_w2,
	  	     double* alpha);

    void crd(const double* inp_w1, const double* inp_w2, double* crd);

    static constexpr double r2i = 5.5;
    static constexpr double r2f = 7.5;
};

} // namespace h2o

#endif // ALP_EL_POLY2B3_H
