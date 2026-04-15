#ifndef TTM4_DIP_H
#define TTM4_DIP_H

#include <cstddef>

namespace ttm {

//
// class incapsulating TTM4 electrostatics
//

struct ttm4_dip {

    // computed for vibrationally averaged (CC-pol) monomer

    static const double molecular_polarizability; // in A^3
    static const double molecular_dipole; // in D

    ttm4_dip();
    ~ttm4_dip();

    void operator()(size_t, const double*, double*);

private:
    size_t  m_nw;
    double* m_mem;

    void allocate(size_t);
};

} // namespace ttm

#endif // TTM4_DIP_H
