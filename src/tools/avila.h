#ifndef AVILA_H
#define AVILA_H

#include <cstddef>

////////////////////////////////////////////////////////////////////////////////

// from http://dx.doi.org/10.1063/1.1867437

namespace x2o {

//----------------------------------------------------------------------------//

struct avila {
    static const double r0; // A
    static const double phi0; // rad

    // in A^3

    static double a_xx(const double& dS1, const double& dS2, const double& dS3);
    static double a_yy(const double& dS1, const double& dS2, const double& dS3);
    static double a_zz(const double& dS1, const double& dS2, const double& dS3);

    static double a_xz(const double& dS1, const double& dS2, const double& dS3);

    // cartesian monomer coordinates
    static void cartesian(const double& dS1, const double& dS2,
                          const double& dS3, double xyz[9]);

    static void alpha(const double xyz[9], double a[9]); // row-major
    static void alpha6(const double xyz[9], double a[6]); // xx, yy, zz, xy, xz, yz
};

//----------------------------------------------------------------------------//

} // namespace x2o

////////////////////////////////////////////////////////////////////////////////

#endif // AVILA_H
