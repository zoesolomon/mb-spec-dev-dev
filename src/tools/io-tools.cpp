#include <cassert>

#include <cstring>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <stdexcept>

#include "io-tools.h"

namespace opt {

//----------------------------------------------------------------------------//

void print_error(size_t lineno, size_t error)
{
    std::ostringstream oss;
    if (error == 0)
	oss << "unexpected text at line " << lineno << " of the XYZ file";
    else if (error == 1)
	oss << "unexpected EOF at line " << lineno << " of the XYZ file";
    else if (error == 2)
	oss << "Matrix is not invertible at line " << lineno << "\n";
    else
	oss << "unexpected error at line " << lineno << " of the XYZ file";

    throw std::runtime_error(oss.str());
}

//----------------------------------------------------------------------------//

void calculate_rcell(size_t lineno, const double* cell, double* rcell)
{

    // Invert unit cell
    std::fill(rcell, rcell + 9, 0.0);

    rcell[0] = cell[4]*cell[8] - cell[5]*cell[7];
    rcell[1] = cell[2]*cell[7] - cell[1]*cell[8];
    rcell[2] = cell[1]*cell[5] - cell[2]*cell[4];
    rcell[3] = cell[5]*cell[6] - cell[3]*cell[8];
    rcell[4] = cell[0]*cell[8] - cell[2]*cell[6];
    rcell[5] = cell[2]*cell[3] - cell[0]*cell[5];
    rcell[6] = cell[3]*cell[7] - cell[4]*cell[6];
    rcell[7] = cell[1]*cell[6] - cell[0]*cell[7];
    rcell[8] = cell[0]*cell[4] - cell[1]*cell[3];

    const double det = cell[0]*rcell[0] + cell[3]*rcell[1] + cell[6]*rcell[2];

    if ( det == 0 )
	print_error(lineno, 2);

    double detinv = 1.0/det;

    for(size_t i = 0; i < 9; ++i)
	rcell[i] *= detinv;

}

//----------------------------------------------------------------------------//


bool match_name(const char* ref_name, const char* test_name)
{
    for(size_t i = 0; i < MAX_CHAR_LENGTH
	&& (ref_name[i] != '\0' || test_name[i] != '\0'); ++i)
	if(ref_name[i] != test_name[i])
	    return false;

    return true;
}

//----------------------------------------------------------------------------//

void print_atom(std::ostream& os, const char* name)
{
    std::string tmpstring(name);
    os << std::setw(8) << tmpstring;
}

////////////////////////////////////////////////////////////////////////////////

} // namespace opt

////////////////////////////////////////////////////////////////////////////////
