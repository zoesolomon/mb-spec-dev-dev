#include <cassert>

#include <cstring>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <stdexcept>

#include "read-lammpstrj.h"
#include "io-tools.h"

////////////////////////////////////////////////////////////////////////////////

namespace opt {
    

//----------------------------------------------------------------------------//

void read_lammpstrj_frame(size_t& lineno, size_t& imcon, std::istream& is, opt::lammpstrj_frame& frame, double t) 
{
    assert(is && !is.eof());
    assert(frame.natm == 0); // empty

    // Creating a junk variable 
    std::string string_junk;
    size_t int_junk;
    double f_junk;
    
    // Reading a line and throwing it
    std::string line;
    std::getline(is, line); 
    if (is.eof()) // check eof before defining/allocating variables
        return;
    ++lineno;

    // Reading a line and saving step number
    std::getline(is, line); 
    if (is.eof()) 
        return;
    ++lineno;
    std::istringstream iss(line);
    iss >> frame.timestep;
    frame.time = t;

    // Reading a line and throwing it
    std::getline(is, line); 
    if (is.eof())
        return;
    ++lineno;

    // Reading a line and saving total atom number
    std::getline(is, line); 
    if (is.eof())
        return;
    ++lineno;
    iss.clear();
    iss.str(line);
    iss >> frame.natm;

    // Reading a line and throwing it
    std::getline(is, line); 
    if (is.eof())
        return;
    ++lineno;

    frame.xyz        = new double[3*frame.natm];
    frame.name       = new char[frame.natm*MAX_CHAR_LENGTH];
    frame.atomic_num = new int[frame.natm];
    frame.mass       = new double[frame.natm];
    frame.box       = new double[6];
    frame.cell       = new double[9];
    frame.rcell       = new double[9];
    frame.dipind       = new double[4*frame.natm];
    frame.dipmol       = new double[4*frame.natm];
    frame.molecid       = new int[frame.natm];

    std::fill(frame.name, frame.name + frame.natm*MAX_CHAR_LENGTH, '0');
    std::fill(frame.cell, frame.cell + 9, 0);

    // read unit cell dimensions (lines 6-8 of frame)
    for(size_t i = 0; i < 3; ++i){

	std::getline(is, line); 
	if (is.eof()) 
	    return;
	++lineno;

        iss.clear();
        iss.str(line);
	iss >> frame.box[2*i] >> frame.box[2*i + 1];
	if (iss.fail()) 
	    opt::print_error(lineno, 0);
	frame.cell[3*i + i] =  frame.box[2*i + 1] - frame.box[2*i];
    }

    if(imcon < 3){		// orthorhombic
	frame.vol = frame.cell[0]*frame.cell[4]*frame.cell[8];

	std::fill(frame.rcell, frame.rcell + 9, 0.0);
	frame.rcell[0] = 1.0/frame.cell[0];
	frame.rcell[4] = 1.0/frame.cell[4];
	frame.rcell[8] = 1.0/frame.cell[8];

    }else{
	std::ostringstream oss;
	oss << "Do not know how to calculate volume for imcon = " 
     	    << imcon << std::endl;
	throw std::runtime_error(oss.str());
    }


    ++lineno;
    std::getline(is, line);

    for (size_t n = 0; n < frame.natm; ++n) {
        if (is.eof()) 
    	    opt::print_error(lineno, 1);

        std::getline(is, line);
        ++lineno;

        iss.clear();
        iss.str(line);
        iss >> int_junk >> frame.molecid[n] >> int_junk >> f_junk 
            >> frame.xyz[3*n + 0] >> frame.xyz[3*n + 1] >> frame.xyz[3*n + 2]
            >> frame.dipmol[4*n + 0] >> frame.dipmol[4*n + 1] >> frame.dipmol[4*n + 2] >> frame.dipmol[4*n + 3]
            >> frame.dipind[4*n + 0] >> frame.dipind[4*n + 1] >> frame.dipind[4*n + 2] >> frame.dipind[4*n + 3];
    
    }
    
    iss.clear();
}

//----------------------------------------------------------------------------//

size_t read_lammpstrj 
	(const char * in, std::vector<opt::lammpstrj_frame>& frames, size_t imcon, double t) {


	std::ifstream ifs(in);
	if (!ifs) {
		std::ostringstream oss;
		oss << "could not open " << in << " for reading\n";
		throw std::runtime_error(oss.str());
	}

	size_t lineno(0);
	size_t frameno(0);

        
	while((!ifs.eof())){
		opt::lammpstrj_frame this_frame;
            
		read_lammpstrj_frame(lineno, imcon, ifs, this_frame, t);
            
		if (ifs.eof() || this_frame.natm == 0) break;
			frames.push_back(this_frame);
			++frameno;
	}

	return frameno;
}

//----------------------------------------------------------------------------//

void read_xyz_frame(size_t& lineno, size_t& imcon, std::istream& is, opt::lammpstrj_frame& frame, double t, std::vector<double> box) 
{
    assert(is && !is.eof());
    assert(frame.natm == 0); // empty

    // Creating a junk variable 
    std::string string_junk;
    size_t int_junk;
    double f_junk;
    
    // Reading a line and saving total atom number
    std::string line;
    std::getline(is, line); 
    if (is.eof())
        return;
    ++lineno;
    std::istringstream iss(line);
    iss.clear();
    iss.str(line);
    iss >> frame.natm;

    // Reading a line and throwing it
    std::getline(is, line); 
    if (is.eof())
        return;
    ++lineno;

    frame.xyz        = new double[3*frame.natm];
    frame.name       = new char[frame.natm*MAX_CHAR_LENGTH];
    frame.atomic_num = new int[frame.natm];
    frame.mass       = new double[frame.natm];
    frame.box       = new double[6];
    frame.cell       = new double[9];
    frame.rcell       = new double[9];
    frame.dipind       = new double[4*frame.natm];
    frame.dipmol       = new double[4*frame.natm];
    frame.molecid       = new int[frame.natm];

    std::fill(frame.name, frame.name + frame.natm*MAX_CHAR_LENGTH, '0');
    std::fill(frame.cell, frame.cell + 9, 0);
    frame.time = t;
    
    for(size_t i = 0; i < 6; ++i){
        frame.box[i] =  box[i];
    }

    for(size_t i = 0; i < 3; ++i){
        frame.cell[3*i + i] =  frame.box[2*i + 1] - frame.box[2*i];
    }

    if(imcon < 3){		// orthorhombic
	frame.vol = frame.cell[0]*frame.cell[4]*frame.cell[8];

	std::fill(frame.rcell, frame.rcell + 9, 0.0);
	frame.rcell[0] = 1.0/frame.cell[0];
	frame.rcell[4] = 1.0/frame.cell[4];
	frame.rcell[8] = 1.0/frame.cell[8];

    }else{
	std::ostringstream oss;
	oss << "Do not know how to calculate volume for imcon = " 
     	    << imcon << std::endl;
	throw std::runtime_error(oss.str());
    }

    for (size_t n = 0; n < frame.natm; ++n) {
        if (is.eof()) 
    	    opt::print_error(lineno, 1);

        std::getline(is, line);
        ++lineno;

        iss.clear();
        iss.str(line);
        iss >> string_junk >> frame.xyz[3*n + 0] >> frame.xyz[3*n + 1] >> frame.xyz[3*n + 2];
        frame.molecid[n] = (n + 1) / 3;
        

    
    }
    
    
    iss.clear();
}


} // namespace

//----------------------------------------------------------------------------//

