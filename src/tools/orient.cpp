#include "orient.h"
#include <algorithm>

////////////////////////////////////////////////////////////////////////////////

namespace tools {

    const double mass_O = 15.9949;
    const double mass_H =  1.0079;

////////////////////////////////////////////////////////////////////////////////

void euler_angles_to_matrix
    (const double& t1, const double& t2, const double& t3, double* m)
{
    const double c1 = std::cos(t1);
    const double s1 = std::sin(t1);

    const double c2 = std::cos(t2);
    const double s2 = std::sin(t2);

    const double c3 = std::cos(t3);
    const double s3 = std::sin(t3);

    m[0] = c1*c2*c3 - s1*s3;
    m[1] = -c2*s3*c1 - c3*s1;
    m[2] = c1*s2;

    m[3] = c1*s3 + c3*c2*s1;
    m[4] = c1*c3 - c2*s1*s3;
    m[5] = s1*s2;

    m[6] = -c3*s2;
    m[7] = s3*s2;
    m[8] = c2;
}

////////////////////////////////////////////////////////////////////////////////

void rotate(const double* m, const double* src, double* dst)
{
    for (size_t i = 0; i < 3; ++i) {
        dst[i] = 0.0;
        for (size_t j = 0; j < 3; ++j)
            dst[i] += m[i*3 + j]*src[j];
    }
}

////////////////////////////////////////////////////////////////////////////////

void cpt_Ilab(const size_t& nw, const double* xyz, double* I)
{

    std::fill(I, I + 9, 0.0);

    for (size_t n = 0; n < nw; ++n){
	for (size_t i = 0; i < 3; ++i){
	    const double mass = ( ( i == 0) ? mass_O : mass_H );

	    const double x = xyz[9*n + 3*i    ];
	    const double y = xyz[9*n + 3*i + 1];
	    const double z = xyz[9*n + 3*i + 2];

	    I[0] += mass*(y*y + z*z); // xx
	    I[1] += mass*(x*y      ); // xy
	    I[2] += mass*(x*z      ); // xz

	    I[3] += mass*(x*y      ); // yx
	    I[4] += mass*(x*x + z*z); // yy
	    I[5] += mass*(y*z      ); // yz

	    I[6] += mass*(x*z      ); // zx
	    I[7] += mass*(y*z      ); // zy
	    I[8] += mass*(x*x + y*y); // zz

	}

    }

}

////////////////////////////////////////////////////////////////////////////////

void center(const size_t& nw, const double* inp, double* xyz)
{

    double com[3] = {0.0};
    double tot_mass = 0.0;

    for (size_t n = 0; n < nw; ++n){
	const size_t n9 = n*9;
	for (size_t i = 0; i < 3; ++i){
	    const double mass = ( ( i == 0) ? mass_O : mass_H );
	    tot_mass += mass;

	    for (size_t j = 0; j < 3; ++j)
		com[j] += mass*inp[n9 + 3*i + j];
	}
    }

    for (size_t i = 0; i < 3; ++i)
	com[i] /= tot_mass;

    for (size_t i = 0; i < 3*nw; ++i){
	const size_t i3 = 3*i;
	for (size_t j = 0; j < 3; ++j){
	    xyz[i3 + j] = inp[i3 + j] - com[j];
	}
    }

}


void cross(const double* y, const double* z, double* x)
{
    x[0] = z[2] * y[1] - z[1] * y[2];
    x[1] = z[0] * y[2] - z[2] * y[0];
    x[2] = z[1] * y[0] - z[0] * y[1];

    double mag = norm(x);

    for(size_t i = 0; i < 3; ++i)
	x[i] /= mag;
}

void norm_vector(const double* y, const double* z, double* x)
{
    for(size_t i = 0; i < 3; ++i)
	x[i] = y[i] - z[i];

    double mag = norm(x);

    for(size_t i = 0; i < 3; ++i)
	x[i] /= mag;
}

////////////////////////////////////////////////////////////////////////////////

} // namespace tools

////////////////////////////////////////////////////////////////////////////////
