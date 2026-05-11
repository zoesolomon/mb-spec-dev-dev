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

    frame.xyz        = new double[3*frame.natm];
    frame.name       = new char[frame.natm*MAX_CHAR_LENGTH];
    frame.atomic_num = new int[frame.natm];
    frame.mass       = new double[frame.natm];
    frame.cell       = new double[9];
    frame.rcell       = new double[9];
    frame.dipind       = new double[4*frame.natm];
    frame.dipmol       = new double[4*frame.natm];
    frame.molecid       = new int[frame.natm];

    std::fill(frame.name, frame.name + frame.natm*MAX_CHAR_LENGTH, '0');
    std::fill(frame.cell, frame.cell + 9, 0);

    // this line should be ITEM: BOX BOUNDS
    // orthorombic - xx yy zz (3 x 2)
    // restricted triclinic (dont need??) - xy xz yz xx yy xx (3 x 3)
    // general triclinic (including monoclinic) - abc origin (3 x 4)
    std::getline(is, line); 
    if (is.eof())
        return;
    ++lineno;

    // read unit cell dimensions (lines 6-8 of frame)

    double* box_temp = new double[12];
    std::fill(box_temp, box_temp + 12, 0);

    for(size_t i = 0; i < 3; ++i){

        std::getline(is, line); 
        if (is.eof()) 
            return;
        ++lineno;

        iss.clear();
        iss.str(line);
        iss >> box_temp[4*i] >> box_temp[4*i + 1] >> box_temp[4*i + 2] >> box_temp[4*i + 3];
        // std::cerr << box_temp[4*i] << " " << box_temp[4*i + 1] << " " << box_temp[4*i + 2] <<  " " << box_temp[4*i + 3] << std::endl;
        
        // if (iss.fail()) 
        //     opt::print_error(lineno, 0);
    }
    

    if (imcon < 3) {

        // check to make sure in lammps orthogonal format
        double temp_sum = box_temp[2] + box_temp[3] + box_temp[6] + box_temp[7] + box_temp[10] + box_temp[11];
        if (temp_sum != 0) {
            throw std::runtime_error("box with imcon = " + std::to_string(imcon) + " must be in orthogonal format");
        }

        frame.cell[0] = box_temp[1] - box_temp[0];
        frame.cell[4] = box_temp[5] - box_temp[4];
        frame.cell[8] = box_temp[9] - box_temp[8];


        frame.vol = frame.cell[0]*frame.cell[4]*frame.cell[8];

        std::fill(frame.rcell, frame.rcell + 9, 0.0);
        frame.rcell[0] = 1.0/frame.cell[0];
        frame.rcell[4] = 1.0/frame.cell[4];
        frame.rcell[8] = 1.0/frame.cell[8];

    } else if (imcon == 3){	// non-orthorhombic

        // check to make sure in lammps triclinic format
        double temp_sum = box_temp[2] + box_temp[3] + box_temp[6] + box_temp[7] + box_temp[10] + box_temp[11];
        double temp2_sum = box_temp[3] + box_temp[7] + box_temp[11];
        if (temp_sum == 0) {
            std::cerr << "box with imcon = " << imcon << "must give tilt factor" << std::endl;
        } else if (temp2_sum != 0) { // restricted triclinic
            double xlo = box_temp[0] - std::min({0.0, box_temp[2], box_temp[6], box_temp[2] + box_temp[6]});
            double xhi = box_temp[1] - std::max({0.0, box_temp[2], box_temp[6], box_temp[2] + box_temp[6]});
            double ylo = box_temp[4] - std::min({0.0, box_temp[10]});
            double yhi = box_temp[5] - std::max({0.0, box_temp[10]});
            
            frame.cell[0] = xhi - xlo;
            frame.cell[3] = box_temp[2];
            frame.cell[4] = yhi - ylo;
            frame.cell[6] = box_temp[6];
            frame.cell[7] = box_temp[10];
            frame.cell[8] = box_temp[9] - box_temp[8];

            frame.vol = frame.cell[0]*frame.cell[4]*frame.cell[8];
        } else {
            for (size_t i = 0; i < 9; i++) {
                if (i == 3 || i == 7 || i == 11) {
                    continue;
                } else {
                    frame.cell[i] = box_temp[i];
                }
            }

            double cross_x = frame.cell[4]*frame.cell[8] - frame.cell[5]*frame.cell[7];
            double cross_y = -frame.cell[1]*frame.cell[8] + frame.cell[7]*frame.cell[2];
            double cross_z = frame.cell[1]*frame.cell[5] - frame.cell[4]*frame.cell[2];

            frame.vol = frame.cell[0]*cross_x + frame.cell[3]*cross_y + frame.cell[6]*cross_z; // [a] . [b] x [c]
        }

	    calculate_rcell(lineno, frame.cell, frame.rcell);

    } else {
	std::ostringstream oss;
	oss << "Do not know how to calculate volume for imcon = " 
     	    << imcon << std::endl;
	throw std::runtime_error(oss.str());
    }
    // std::cerr << frame.cell[0] << " " << frame.cell[1] << " " << frame.cell[2] <<  " " << 
    //         frame.cell[3] << " " << frame.cell[4] << " " << frame.cell[5] << " " <<
    //         frame.cell[6] << " " << frame.cell[7] << " " << frame.cell[8] << " " <<
    //         std::endl;

    // ++lineno;
    // std::getline(is, line);

    // Reading a line and throwing it
    std::getline(is, line); 
    if (is.eof())
        return;
    ++lineno;

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
	(const char * in, std::vector<opt::lammpstrj_frame>& frames, size_t imcon, double dump) {


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
            
		read_lammpstrj_frame(lineno, imcon, ifs, this_frame, dump * frameno);
            
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
    frame.cell       = new double[9];
    frame.rcell       = new double[9];
    frame.dipind       = new double[4*frame.natm];
    frame.dipmol       = new double[4*frame.natm];
    frame.molecid       = new int[frame.natm];

    std::fill(frame.name, frame.name + frame.natm*MAX_CHAR_LENGTH, '0');
    std::fill(frame.cell, frame.cell + 9, 0);
    frame.time = t;
    
    for(size_t i = 0; i < 9; ++i){
        frame.cell[i] =  box[i];
    }

    if (imcon < 3) {

        frame.vol = frame.cell[0]*frame.cell[4]*frame.cell[8];

        std::fill(frame.rcell, frame.rcell + 9, 0.0);
        frame.rcell[0] = 1.0/frame.cell[0];
        frame.rcell[4] = 1.0/frame.cell[4];
        frame.rcell[8] = 1.0/frame.cell[8];

    } else if (imcon == 3){	// non-orthorhombic

        double cross_x = frame.cell[4]*frame.cell[8] - frame.cell[5]*frame.cell[7];
        double cross_y = -frame.cell[1]*frame.cell[8] + frame.cell[7]*frame.cell[2];
        double cross_z = frame.cell[1]*frame.cell[5] - frame.cell[4]*frame.cell[2];

        frame.vol = frame.cell[0]*cross_x + frame.cell[3]*cross_y + frame.cell[6]*cross_z; // [a] . [b] x [c]

	    calculate_rcell(lineno, frame.cell, frame.rcell);

    } else {
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

