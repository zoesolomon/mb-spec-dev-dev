#ifndef STRUCT_LAMMPSTRJ_FRAME_H
#define STRUCT_LAMMPSTRJ_FRAME_H

#include <vector>
//#include <string>
//#include <iostream>

#include "macros.h"

namespace opt {

struct lammpstrj_frame {
    inline lammpstrj_frame();
    inline lammpstrj_frame(const lammpstrj_frame&);
    inline ~lammpstrj_frame();
    inline lammpstrj_frame& operator=(const lammpstrj_frame&);

    size_t natm;
    double* xyz; // 3*natm
    char* name; // natm * MAX_CHAR_LENGTH
    int* atomic_num; // natm;
    double* mass; // natm
    // double* box;
    double* cell;
    double* rcell;
    double* dipind;
    double* dipmol;
    double vol;
    size_t timestep;
    double time;
    size_t imcon;
    int* molecid;
};

inline lammpstrj_frame::lammpstrj_frame()
: natm(0), xyz(0), name(0), atomic_num(0), mass(0), cell(0), rcell(0),
//  natm(0), xyz(0), name(0), atomic_num(0), mass(0), box(0), cell(0), rcell(0),
    dipind(0), dipmol(0), vol(0), 
    timestep(0), time(0), imcon(0), molecid(0)
{
}

inline lammpstrj_frame::lammpstrj_frame(const lammpstrj_frame& o)
{
    natm = o.natm;
    if (natm > 0) {
        xyz = new double[3*natm];
        std::copy(o.xyz, o.xyz + 3*natm, xyz);
	name = new char[natm*MAX_CHAR_LENGTH];
    	std::copy(o.name, o.name + natm*MAX_CHAR_LENGTH, name);
	atomic_num = new int[natm];
    	std::copy(o.atomic_num, o.atomic_num + natm, atomic_num);
    	mass = new double[natm];
    	std::copy(o.mass, o.mass + natm, mass);
    	cell = new double[9];
    	std::copy(o.cell, o.cell + 9, cell);
    	rcell = new double[9];
    	std::copy(o.rcell, o.rcell + 9, rcell);
        // box = new double[9];
        // std::copy(o.box, o.box + 9, box);
        dipind = new double[4*natm];
        std::copy(o.dipind, o.dipind + 4*natm, dipind);
        dipmol = new double[4*natm];
        std::copy(o.dipmol, o.dipmol + 4*natm, dipmol);
	vol = o.vol;
	timestep = o.timestep;
	time = o.time;
	imcon = o.imcon;
	molecid = new int[natm];
	std::copy(o.molecid, o.molecid + natm, molecid);
    }
}

inline lammpstrj_frame::~lammpstrj_frame()
{
    if (natm > 0){
     	delete[] xyz;
	delete[] name;
	delete[] atomic_num;
	delete[] mass;
        // delete[] box;
	delete[] cell;
	delete[] rcell;
	delete[] dipind;
	delete[] dipmol;
        delete[] molecid;
    }
}

inline lammpstrj_frame& lammpstrj_frame::operator=(const lammpstrj_frame& o)
{
    if (this != &o) {
        if (natm > 0){
            delete[] xyz;
	    delete[] name;
            delete[] atomic_num;
            delete[] mass;
            // delete[] box;
	    delete[] cell;
	    delete[] rcell;
	    delete[] dipind;
	    delete[] dipmol;
            delete[] molecid;
	}

        natm = o.natm;
        if (natm > 0) {
            xyz = new double[3*natm];
            std::copy(o.xyz, o.xyz + 3*natm, xyz);
	    name = new char[natm*MAX_CHAR_LENGTH];
	    std::copy(o.name, o.name + natm*MAX_CHAR_LENGTH, name);
            atomic_num = new int[natm];
            std::copy(o.atomic_num, o.atomic_num + natm, atomic_num);
            mass = new double[natm];
            std::copy(o.mass, o.mass + natm, mass);
            // box = new double[9];
            // std::copy(o.box, o.box + 9, box);
	    cell = new double[9];
	    std::copy(o.cell, o.cell + 9, cell);
	    rcell = new double[9];
	    std::copy(o.rcell, o.rcell + 9, rcell);
	    dipind = new double[4*natm];
	    std::copy(o.dipind, o.dipind + 4*natm, dipind);
	    dipmol = new double[4*natm];
	    std::copy(o.dipmol, o.dipmol + 4*natm, dipmol);
	    vol = o.vol;
	    timestep = o.timestep;
	    time = o.time;
	    imcon = o.imcon;
	    molecid = new int[natm];
	    std::copy(o.molecid, o.molecid + natm, molecid);
	}
    }

    return *this;
}

} // namespace opt

#endif // STRUCT_FRAME_H
