#ifndef PBC_TOOLS_H
#define PBC_TOOLS_H

#include <cstddef>

namespace tools {

    void image(const double*, const double*, const double*, double*);
    void image3(const size_t, const double*, const double*, const double*, double*);

    void center_slab_at_origin(const size_t, double*, double*,
    			       double*, bool shift_half_cell = false);
    void center_slab_at_origin3(const size_t, const size_t, double*, double*,
    			       double*, double*, bool shift_half_cell = false);

    void dump_frame(double*, const int, double*);

    double distance3(size_t, const double*, const double*,
		     const double*, const double*, double* dr = 0);

    void image_monomer(const double*, double*);
    void image_monomer(const double*, const double*, double*);

} // namespace tools

#endif // PBC_TOOLS_H
