#include <algorithm>

#include "alpha-ttm4.h"
#include "ttm4-es.h"

namespace ttm {

void alpha_tot(size_t nw, double* xyz, double* pt)
{

    ttm::ttm4_es ind_model;

    double Eelec, Eind;
    ind_model(nw, xyz, Eelec, 0, Eind, 0);

    ind_model.ptensor(nw, pt);

}

void alpha_no_1b(size_t nw, double* xyz, double* pt)
{

    // Calculate polarizability of total system
    {

    	ttm::ttm4_es ind_model;

    	double Eelec, Eind;
	ind_model(nw, xyz, Eelec, 0, Eind, 0);

    	ind_model.ptensor(nw, pt);

    }

    // remove 1-body polarizability 
    for(size_t n = 0; n < nw; ++n){

    	ttm::ttm4_es ind_model;

    	double Eelec, Eind;
	ind_model(1, xyz + 9*n, Eelec, 0, Eind, 0);

	double pt_1b[9];
    	ind_model.ptensor(1, pt_1b);

	for(size_t i = 0; i < 9; ++i)
	    pt[i] -= pt_1b[i];

    }

}

} // namespace ttm
