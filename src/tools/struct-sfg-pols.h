#ifndef STRUCT_SFG_POLS_H
#define STRUCT_SFG_POLS_H

#include <vector>
#include <string>
#include <iostream>

namespace sfg {

struct eprop {
    inline eprop();
    inline eprop(const eprop&);
    inline ~eprop();

    size_t nmolec;
    double* dip_k;  // nmolec
    double* dip2_k;  // nmolec
    double* pol_ij;  // nmolec
    double* pol2_ij;  // nmolec
    double* gscale;  // nmolec
};

inline eprop::eprop()
: nmolec(0), dip_k(0), dip2_k(0), pol_ij(0), pol2_ij(0), gscale(0)
{
}

inline eprop::eprop(const eprop& o)
{
    nmolec  = o.nmolec;
    if (nmolec > 0) {
	dip_k = new double[nmolec];
    	std::copy(o.dip_k, o.dip_k + nmolec, dip_k);
	dip2_k = new double[nmolec];
    	std::copy(o.dip2_k, o.dip2_k + nmolec, dip2_k);
	pol_ij = new double[nmolec];
    	std::copy(o.pol_ij, o.pol_ij + nmolec, pol_ij);
	pol2_ij = new double[nmolec];
    	std::copy(o.pol2_ij, o.pol2_ij + nmolec, pol2_ij);
	gscale = new double[nmolec];
    	std::copy(o.gscale, o.gscale + nmolec, gscale);
    }
}

inline eprop::~eprop()
{
    if (nmolec > 0){
	delete[] dip_k;
	delete[] dip2_k;
	delete[] pol_ij;
	delete[] pol2_ij;
	delete[] gscale;
    }
    
}

struct connect {
    inline connect();

    size_t nmolec;
    std::vector< std::vector<int> > partners;
    std::vector<int> central;
};

inline connect::connect()
: nmolec(0), partners(0), central(0)
{
}

} // namespace sfg

#endif // STRUCT_SFG_POLS_H
