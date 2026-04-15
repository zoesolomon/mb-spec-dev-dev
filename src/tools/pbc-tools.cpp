#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>

#include "pbc-tools.h"

////////////////////////////////////////////////////////////////////////////////

namespace tools {

////////////////////////////////////////////////////////////////////////////////

void image(const double* cell, const double* r1, const double* r2,
	   double* crd_out)
{
    double dx = r1[0] - r2[0];
    double dy = r1[1] - r2[1];
    double dz = r1[2] - r2[2];

    double facx = std::floor(dx/cell[0] + 0.5);
    double facy = std::floor(dy/cell[4] + 0.5);
    double facz = std::floor(dz/cell[8] + 0.5);

    crd_out[0] = r2[0] + cell[0]*facx;
    crd_out[1] = r2[1] + cell[4]*facy;
    crd_out[2] = r2[2] + cell[8]*facz;
}

//----------------------------------------------------------------------------//

void image3(const size_t imcon, const double* cell, const double* rcell,
           const double* r1, const double* r2, double* crd_out)
{
//    assert(imcon == 0 || imcon == 1 || imcon == 2 || imcon == 3);

    double dx = r1[0] - r2[0];
    double dy = r1[1] - r2[1];
    double dz = r1[2] - r2[2];

    if (imcon == 1 || imcon == 2) {

       dx -= cell[0]*std::floor(dx/cell[0] + 0.5);
       dy -= cell[4]*std::floor(dy/cell[4] + 0.5);
       dz -= cell[8]*std::floor(dz/cell[8] + 0.5);

       crd_out[0] = r1[0] - dx;
       crd_out[1] = r1[1] - dy;
       crd_out[2] = r1[2] - dz;

    } else if (imcon == 3) {

       const double ssx = rcell[0]*dx + rcell[3]*dy + rcell[6]*dz;
       const double ssy = rcell[1]*dx + rcell[4]*dy + rcell[7]*dz;
       const double ssz = rcell[2]*dx + rcell[5]*dy + rcell[8]*dz;

       const double xss = ssx - std::floor(ssx + 0.5);
       const double yss = ssy - std::floor(ssy + 0.5);
       const double zss = ssz - std::floor(ssz + 0.5);

       dx = cell[0]*xss + cell[3]*yss + cell[6]*zss;
       dy = cell[1]*xss + cell[4]*yss + cell[7]*zss;
       dz = cell[2]*xss + cell[5]*yss + cell[8]*zss;

       crd_out[0] = r1[0] - dx;
       crd_out[1] = r1[1] - dy;
       crd_out[2] = r1[2] - dz;

    }

}

//----------------------------------------------------------------------------//

void center_slab_at_origin(const size_t natm, double* cell, double* mass,
			   double* xyz, bool shift_half_cell)
{
    // First, image all molecules in z-direction relative to z=0
    // (takes care of any 'inverted slab' artifact of PBC)
    double origin[3] = {0.0, 0.0, 0.0};
    if(shift_half_cell){
	for(size_t i = 0; i < natm; ++i)
	    xyz[3*i + 2] += cell[8]/2.0;
    }

    for(size_t i = 0; i < natm; ++i){
	double tmp_atom[3];
	std::copy(xyz + 3*i, xyz + 3*i + 3,
		  tmp_atom);

        image(cell, origin, tmp_atom, xyz + 3*i);
    }

    // set the center of mass of the cell as the origin
    double com[3] = {0.0, 0.0, 0.0};
    double tot_mass(0);
    for(size_t i = 0; i < natm; ++i){
	for(size_t j = 0; j < 3; ++j)
	    com[j] += mass[i]*xyz[i*3 + j];

	tot_mass += mass[i];
    }

    for(size_t j = 0; j < 3; ++j)
	com[j] /= tot_mass;

    for(size_t i = 0; i < natm; ++i)
	for(size_t j = 0; j < 3; ++j)
	    xyz[i*3 + j] -= com[j];

}

//----------------------------------------------------------------------------//

void center_slab_at_origin3(const size_t natm, const size_t imcon, double* cell, 
                double* rcell, double* mass, double* xyz, bool shift_half_cell)
{
    // First, image all molecules in z-direction relative to z=0
    // (takes care of any 'inverted slab' artifact of PBC)
    double origin[3] = {0.0, 0.0, 0.0};
    if(shift_half_cell){
	for(size_t i = 0; i < natm; ++i)
	    xyz[3*i + 2] += cell[8]/2.0;
    }

    for(size_t i = 0; i < natm; ++i){
	double tmp_atom[3];
	std::copy(xyz + 3*i, xyz + 3*i + 3,
		  tmp_atom);

        image3(imcon, cell, rcell, origin, tmp_atom, xyz + 3*i);
    }

    // set the center of mass of the cell as the origin
    double com[3] = {0.0, 0.0, 0.0};
    double tot_mass(0);
    for(size_t i = 0; i < natm; ++i){
	for(size_t j = 0; j < 3; ++j)
	    com[j] += mass[i]*xyz[i*3 + j];

	tot_mass += mass[i];
    }

    for(size_t j = 0; j < 3; ++j)
	com[j] /= tot_mass;

    for(size_t i = 0; i < natm; ++i)
	for(size_t j = 0; j < 3; ++j)
	    xyz[i*3 + j] -= com[j];

}

//----------------------------------------------------------------------------//

void dump_frame(double* cell, const size_t natm, double* xyz)
{

    std::cout << "timestep      1000       930         0         2    0.000000"
	<< std::endl
	<< "   22.00       0.000       0.000\n"
	<< "   0.000       22.00       0.000\n"
	<< "   0.000       0.000       100.0\n";
    for(size_t i = 0; i < natm; ++i){
	if(i%3 == 0)
	    std::cout << "OW";
	else
	    std::cout << "HW";
	std::cout << "             " << i << "    0.000000    0.000000\n";
	for(size_t j = 0; j < 3; ++j)
	    std::cout << std::setw(20) << xyz[3*i + j];
	std::cout << std::endl;
    }

}

//----------------------------------------------------------------------------//

double distance3(size_t imcon, const double* cell, const double* rcell,
		 const double* r1, const double* r2, double* dr)
{
    double dx = r1[0] - r2[0];
    double dy = r1[1] - r2[1];
    double dz = r1[2] - r2[2];

    const double ssx = rcell[0]*dx + rcell[3]*dy + rcell[6]*dz;
    const double ssy = rcell[1]*dx + rcell[4]*dy + rcell[7]*dz;
    const double ssz = rcell[2]*dx + rcell[5]*dy + rcell[8]*dz;

    const double xss = ssx - std::floor(ssx + 0.5);
    const double yss = ssy - std::floor(ssy + 0.5);
    const double zss = ssz - std::floor(ssz + 0.5);

    dx = cell[0]*xss + cell[3]*yss + cell[6]*zss;
    dy = cell[1]*xss + cell[4]*yss + cell[7]*zss;
    dz = cell[2]*xss + cell[5]*yss + cell[8]*zss;

    if(dr){
	dr[0] = dx;
	dr[1] = dy;
	dr[2] = dz;
    }

    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

//----------------------------------------------------------------------------//

void image_monomer(const double* cell,
		   double* mon)
{
    double r1[9];
    std::copy(mon, mon + 9, r1);

    image_monomer(cell, r1, mon);

}
//----------------------------------------------------------------------------//

void image_monomer(const double* cell,
		   const double* r1, double* mon)
{
    std::fill(mon, mon+9, 0.0);

    // take r1 oxygen as central atom, image its two hydrogens

    for(size_t i = 0; i < 3; ++i)
	mon[i] = r1[i];

    for(size_t a = 1; a < 3; ++a)
	image(cell, r1, r1 + 3*a, mon + 3*a);
}

////////////////////////////////////////////////////////////////////////////////

} // namespace tools
