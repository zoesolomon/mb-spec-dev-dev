#include "lp-site.h"

namespace mon {

void lp_site(const double* ohh,
	     const double& in_plane_g, const double& out_of_plane_g,
	     double* ohhll)
{
    double oh1[3];
    double oh2[3];

    for (int i = 0; i < 3; ++i) {
	oh1[i] = ohh[i + 3] - ohh[i];
 	oh2[i] = ohh[i + 6] - ohh[i];
    }

    const double v[3] = {
	oh1[1]*oh2[2] - oh1[2]*oh2[1],
	oh1[2]*oh2[0] - oh1[0]*oh2[2],
	oh1[0]*oh2[1] - oh1[1]*oh2[0]
    };

    for (int i = 0; i < 9; ++i) {
	ohhll[i] = ohh[i];
    }

    for (int i = 0; i < 3; ++i) {
	const double in_plane = ohh[i] + 0.5*in_plane_g*(oh1[i] + oh2[i]);
	const double out_of_plane = out_of_plane_g*v[i];

	ohhll[ 9 + i] = in_plane + out_of_plane;
	ohhll[12 + i] = in_plane - out_of_plane;
    }
}

} // namespace mon
