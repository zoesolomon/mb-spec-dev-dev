#ifndef TTM4_PERIODIC_ES_H
#define TTM4_PERIODIC_ES_H

#include <cstddef>

namespace ttm {

//
// class incapsulating TTM4 electrostatics
//

struct ttm4_periodic_es {

    // computed for vibrationally averaged (CC-pol) monomer

    static const double molecular_polarizability; // in A^3
    static const double molecular_dipole; // in D

    ttm4_periodic_es();
    ~ttm4_periodic_es();

    double operator()(size_t nw,
		      const double*, // O H H O H H
		      const double*);

    double lapack(size_t nw,
		  const double*, // O H H O H H
	      	  const double*);


    //Calculates Total Dipole
    void dip(const size_t&, double*);

    //Calculates 2B dipole
    void dipind(const size_t&, double*);
    
    //Calculates polarizability
    void ptensor(const size_t&, double*);


private:
    size_t  m_nw;
    double* m_mem;

    void allocate(size_t);
};

} // namespace ttm

#endif // TTM4_PERIODIC_ES_H
