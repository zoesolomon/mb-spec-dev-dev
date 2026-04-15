#ifndef ORIENT_H
#define ORIENT_H

#include <cmath>
#include <cstddef>

////////////////////////////////////////////////////////////////////////////////

namespace tools {

//----------------------------------------------------------------------------//

inline double distance(const double* p1, const double* p2)
{
    double result(0);

    for (size_t k = 0; k < 3; ++k) {
        const double delta = p1[k] - p2[k];
        result += delta*delta;
    }

    return std::sqrt(result);
}

inline double dot(const double* a, const double* b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

inline double norm(const double* v)
{
    return std::sqrt(dot(v, v));
}

inline double angle_vect(const double* v, const double* u)
{
    return 180.0/M_PI*std::acos(dot(v, u)/norm(u)/norm(v));
}

inline double angle(const double* a, const double* b, const double* c)
{
    double v[3];
    double u[3];
    for(size_t i = 0; i < 3; ++i){
	v[i] = a[i] - b[i];
	u[i] = c[i] - b[i];
    }
    return angle_vect(v, u);
}


inline void convert_2_unit_vector(double* x)
{
    double my_norm = norm(x);

    for(size_t i = 0; i < 3; ++i)
	x[i] /= my_norm;
}

//----------------------------------------------------------------------------//

// ZYZ
void euler_angles_to_matrix
    (const double& t1, // 0 ... 2*\pi
     const double& t2, // 0 ...   \pi
     const double& t3, // 0 ... 2*\pi
     double* row_major_3x3);

void rotate(const double*, const double*, double*);

void cpt_Ilab(const size_t&, const double*, double*);

void center(const size_t&, const double*, double*);

//----------------------------------------------------------------------------//

void cross(const double*, const double*, double*);

void norm_vector(const double*, const double*, double*);

//----------------------------------------------------------------------------//

} // namespace tools

////////////////////////////////////////////////////////////////////////////////

#endif // ORIENT_H
