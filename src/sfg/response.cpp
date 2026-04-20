//////////////////////////////////////////////////////////////////////////
//  Calculates SFG response function for given polarization condition.  //
//  Options for polarizations combinations are: zzz, zxx, xzx, and xxz. //
//////////////////////////////////////////////////////////////////////////
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <omp.h>

#include <cassert>
#include <cstdlib>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <cmath>
#include <algorithm>

#include "read-lammpstrj.h"
#include "system.h"
#include "pbc-tools.h"
#include "struct-sfg-pols.h"
#define DO_AC yes  	// on to include auto-molecular correlations
#define DO_CC yes       // on to include cross-molecular correlations
//#define SHIFT_HALF_CELL yes
//#define CENTER_CELL yes
//#define BINARY yes

//#define TYPE2 yes
//#define JUST_SWITCH yes

//#define PS_SWITCH yes // on to use PS DMS, off to read from DIPMOL_CMD

#include "ps.h"
#include "avila.h"

////////////////////////////////////////////////////////////////////////////////

namespace {

////////////////////////////////////////////////////////////////////////////////

std::vector<opt::lammpstrj_frame> frames;
static x2o::avila pol_model;

std::vector<sfg::eprop> eprops;
std::vector<sfg::connect> connects;

//const double zc1(1.0);
//const double zc2(2.0);
const double rOO_max(4);

//----------------------------------------------------------------------------//

inline double g_scale(const double* crd, const double zc1, const double zc2)
{
    double z = std::fabs(crd[2]);
    if(z < zc1)
	return 0.0;
    else if(z > zc2)
	return copysign(1.0, crd[2]);
    else{
	const double c = std::cos(M_PI*(z-zc2)/(2*(zc1-zc2)));
	return copysign(c*c, crd[2]); // multiply switch by sign of z (crd[2])
    }
}

//----------------------------------------------------------------------------//

inline void image(const size_t imcon, const double* cell, const double* rcell,
           const double* r1, const double* r2, double* crd_out)
{
    assert(imcon == 0 || imcon == 1 || imcon == 2 || imcon == 3);

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

inline void image_monomer(const size_t imcon, const double* cell, const double* rcell,
       const double* r1, double* mon)
{
    std::fill(mon, mon + 9, 0.0);

    // take r1 oxygen as central atom, image its two hydrogens

    for(size_t i = 0; i < 3; ++i)
        mon[i] = r1[i];

    for(size_t a = 1; a < 3; ++a)
        image(imcon, cell, rcell, r1, r1 + 3*a, mon + 3*a);
}

////////////////////////////////////////////////////////////////////////////////

} // namespace


////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
    if (argc < 3) {
        std::cerr << "usage: " << argv[0] << " <file_name(s)> <mb-spec.json>" << std::endl;
        return EXIT_FAILURE;
    }

    System sys;
    sys.SetUpFromJson("SFG", argv[argc - 1]);

    size_t threads(0);
    {
	std::istringstream iss(argv[1]);
	iss >> threads;
	if (!iss || !iss.eof()) {
	    std::cerr << "could not convert '" << argv[1]
		<< "' to integer" << std::endl;
	    return EXIT_FAILURE;
	}
    }

    double zc1(0);
    {
	zc1 = sys.GetZc1();
	if (!std::is_same<decltype(zc1), double>::value) {
	    std::cerr << "could not convert '" << zc1
		<< "' to real number" << std::endl;
	    return EXIT_FAILURE;
	}
    }
    
    double zc2(0);
    {
	zc2 = sys.GetZc2();
	if (!std::is_same<decltype(zc2), double>::value) {
	    std::cerr << "could not convert '" << zc2
		<< "' to real number" << std::endl;
	    return EXIT_FAILURE;
	}
    }
    
    double int_max(0);
    {
	int_max = sys.GetIntMax();
	if (!std::is_same<decltype(int_max), double>::value) {
	    std::cerr << "could not convert '" << int_max
		<< "' to real number" << std::endl;
	    return EXIT_FAILURE;
	}
    }

    std::string ijk = sys.GetIjk();
    std::transform(ijk.begin(), ijk.end(), ijk.begin(), ::tolower);

    int dipk(0);
    int dip2k(0);
    int polij(0);
    int pol2ij(0);

    std::cerr << "ijk = " << ijk << std::endl;

    // Define possible polarization combinations

    if (ijk.compare("zzz") == 0) {
        dipk = 2;		// z component of dipole vector
        dip2k = 2;		// z component of dipole vector
        polij = 8;		// zz element of polarizabilty tensor
        pol2ij = 8;		// zz element of polarizabilty tensor
    } else if (ijk.compare("zxx") == 0) { //|| (ijk.compare("zyy") == 0) {
        dipk = 0;		// x component of dipole vector
        dip2k = 1;		// y component of dipole vector
        polij = 6;		// zx element of polarizability tensor
        pol2ij = 7;		// zy element of polarizability tensor
    } else if (ijk.compare("xzx") == 0) { //|| (ijk.compare("yzy") == 0) {
        dipk = 0;		// x component of dipole vector
        dip2k = 1;		// y component of dipole vector
        polij = 2;		// xz element of polarizability tensor
        pol2ij = 5;		// yz element of polarizability tensor
    } else if (ijk.compare("xxz") == 0) { //|| (ijk.compare("yyz") == 0) {
        dipk = 2;		// z component of dipole vector
        dip2k = 2;		// z component of dipole vector
        polij = 0;		// xx element of polarizability tensor
        pol2ij = 4;		// yy element of polarizability tensor
    } else {
	std::cerr << "not a valid polarization combination" << std::endl;
	return EXIT_FAILURE;
    }

//    std::cerr << "dip = " << dipk << std::endl;
//    std::cerr << "pol = " << polij << std::endl;
//    std::cerr << "dip2 = " << dip2k << std::endl;
//    std::cerr << "pol2 = " << pol2ij << std::endl;

    std::string nb = sys.GetMb();
    std::transform(nb.begin(), nb.end(), nb.begin(), ::tolower);

    size_t imcon(0);
    {
	imcon = sys.GetImcon();
    if (imcon == 0){
      imcon = 2;
    } else if (imcon < 0 || imcon > 3) {
      std::cerr << "Cannot calculate system of imcon = " << imcon << std::endl;
      return EXIT_FAILURE;
    } else if (imcon > 0 && imcon < 3) {
      if (imcon == 1) {
        std::cerr << "imcon = " << imcon << " (cubic)" << std::endl;
      } else if (imcon == 2) {
        std::cerr << "imcon = " << imcon << " (orthorhombic)" << std::endl;
      } else if (imcon == 3) {
        std::cerr << "imcon = " << imcon << " (monoclinic)" << std::endl;
      }
    } else {
      std::cerr << "imcon value of " << imcon << " not a valid integer" << std::endl;
      return EXIT_FAILURE;
    }
    
    double dump(0);
    {
	dump = sys.GetDumpFrequency();
	if (!std::is_same<decltype(dump), double>::value) {
	    std::cerr << "could not convert '" << dump
		<< "' to real number" << std::endl;
	    return EXIT_FAILURE;
	}
    }

    for(int i = 9; i < argc; ++i){
	try {
	    opt::read_lammpstrj(argv[i], frames, imcon, dump);
	    std::cerr <<" # Frames : " << frames.size() << std::endl;

	} catch (const std::exception& e) {
	    std::cerr << e.what() << std::endl;
	    return EXIT_FAILURE;
	}
    }

    // Open files for Writing

    // First, go over trajectory and pre-compute \mu_k and \alpha_{ij}

    double dt = (frames[1].time - frames[0].time);
    int int_time_ind = (int) (int_max / dt);

    std::cerr << dt << "   " << int_time_ind << std::endl;

    #ifdef CENTER_CELL
    #pragma omp parallel for num_threads(threads)
    for (size_t n = 0; n < frames.size(); ++n) {
        #ifdef SHIFT_HALF_CELL
           if (imcon < 3) {
               tools::center_slab_at_origin(frames[n].natm, 
                  frames[n].cell, frames[n].mass, frames[n].xyz, true);
           } else if (imcon == 3) {
               tools::center_slab_at_origin3(frames[n].natm, imcon, 
                  frames[n].cell, frames[n].rcell, frames[n].mass, frames[n].xyz, true);
           } else {
               std::cerr << "Do not know how to center cell of imcon = " << imcon << std::endl;
           }
        #else
           if (imcon < 3) {
               tools::center_slab_at_origin(frames[n].natm, 
                  frames[n].cell, frames[n].mass, frames[n].xyz, false);
           } else if (imcon == 3) {
               tools::center_slab_at_origin3(frames[n].natm, imcon, 
                  frames[n].cell, frames[n].rcell, frames[n].mass, frames[n].xyz, false);
           } else {
               std::cerr << "Do not know how to center cell of imcon = " << imcon << std::endl;
           }
        #endif
    }
    #endif

    for(size_t n = 0; n < frames.size(); ++n) {
 
        const size_t nmolec = frames[n].natm/3;
        { // calculate electrostatic properties
            sfg::eprop eprop_this_frame;
            eprop_this_frame.nmolec  = nmolec;
            eprop_this_frame.dip_k   = new double[eprop_this_frame.nmolec];
            eprop_this_frame.dip2_k  = new double[eprop_this_frame.nmolec];
            eprop_this_frame.pol_ij  = new double[eprop_this_frame.nmolec];
            eprop_this_frame.pol2_ij = new double[eprop_this_frame.nmolec];
            eprop_this_frame.gscale  = new double[eprop_this_frame.nmolec];
 
            #pragma omp parallel for num_threads(threads)
            for(size_t i = 0; i < nmolec; ++i){
                const size_t i9 = 9*i;
 
                double mon[9];
                image_monomer(imcon, frames[n].cell, frames[n].rcell, frames[n].xyz + i9, mon);
 
            #ifdef JUST_SWITCH
                eprop_this_frame.dip_k[i] = 1.0;
                eprop_this_frame.dip2_k[i] = 1.0;
                eprop_this_frame.pol_ij[i] = 1.0;
                eprop_this_frame.pol2_ij[i] = 1.0;
            #else
 	
                // k-component of the dipole
            #ifdef PS_SWITCH
                double tmp_dip[3];
                ps::dipole(mon, tmp_dip);
                double moldip_k = tmp_dip[dipk];
                double moldip2_k = tmp_dip[dip2k];
            #else
                double moldip_k = frames[n].dipmol[3*4*i+dipk];
                double moldip2_k = frames[n].dipmol[3*4*i+dip2k];
            #endif
 
                if (nb.compare("1b") != 0) {	 // 1B+NB dipoles
                    double inddip_k = frames[n].dipind[3*4*i+dipk];
                    double inddip2_k = frames[n].dipind[3*4*i+dip2k];
                    eprop_this_frame.dip_k[i] = moldip_k + inddip_k;
                    eprop_this_frame.dip2_k[i] = moldip2_k + inddip2_k;
                } else {			         // 1B dipoles
                    eprop_this_frame.dip_k[i] = moldip_k;
                    eprop_this_frame.dip2_k[i] = moldip2_k;
                }
 
                // ij-component of the polarizability
                double pol_mol[9];
      	        std::fill(pol_mol, pol_mol + 9, 0.0);
                pol_model.alpha(mon, pol_mol);
                eprop_this_frame.pol_ij[i] = pol_mol[polij];
                eprop_this_frame.pol2_ij[i] = pol_mol[pol2ij];
     
                #endif
     		
                // scaling that accounts for the double-interface
                double g_i = g_scale(frames[n].xyz + i9, zc1, zc2);
                eprop_this_frame.gscale[i] = g_i;
            }
 
            eprops.push_back(eprop_this_frame);
 
        } // end electrostatics calculation
 
        { // determine connectivity
 
            sfg::connect connect_this_frame;
            connect_this_frame.nmolec = nmolec;
     
            for (size_t i = 0; i < nmolec; ++i) { 
                if (eprops[n].gscale[i] != 0.0) {
                    #ifdef DO_CC //Begin partner determination
                        const size_t i9 = 9*i;
                        std::vector<int> mol;
     
                        #pragma omp parallel for num_threads(threads)
                        for(size_t j = 0; j < nmolec; ++j){
                            double avg_rij(0);
                            size_t count(0);
                            if (j != i){
                            // Assume that int_time_ind is short enough to 
                            // prohibit diffusion between solvation shells
                            // and that the average_rij cannot be < rOO_max
                            // unless rij at the initial time is 
                            // less than 1.5*rOO_max
                            if (tools::distance3(imcon, frames[n].cell,
                                frames[n].rcell,
                                frames[n].xyz + i9,
                                frames[n].xyz + j*9)
                                < 2.0*rOO_max){
     
                                // calculate average rij distance over the next
                                // int_time_ind frames
                                for(size_t a = n;
                                    a < std::min(frames.size(), n + int_time_ind + 1);
                                    a += 250){
         
                                    #pragma omp atomic
                                    avg_rij += tools::distance3(imcon,
                                                 frames[n].cell,
                                                 frames[n].rcell,
                                                 frames[n].xyz + i9,
                                                 frames[a].xyz + j*9);
         
                                    #pragma omp atomic
                                    ++count;
                                }
          
                                if ((avg_rij/count) < rOO_max) {
                                    #pragma omp critical
                                    {
                                    mol.push_back((int)j);
                                    }
                                }
     		            }
     		        }
     	            }
     
                    connect_this_frame.partners.push_back(mol);
                    #endif // end partner determination

                    connect_this_frame.central.push_back(i);
     	        }
            }
     
      	    connects.push_back(connect_this_frame);
        }
    }
 
    if (eprops.size() != frames.size()) {
        std::cerr << "eprops.size() != frames.size()" << std::endl;
        exit(1);
    }
 
    std::cerr << "Deallocating frames" << std::endl;
    //frames.erase(frames.begin(), frames.end());
    std::vector<opt::lammpstrj_frame>().swap(frames);
 
    // With electrostatic properties loaded in eprops
    // Calculated the 1B approximation to the SFG with the 
    // truncated cross-correlation approximation of dx.doi.org/10.1021/jz400683v
     
    double corr[int_time_ind];
    std::fill(corr, corr + int_time_ind, 0.0);
    int ntime[int_time_ind];
    std::fill(ntime, ntime + int_time_ind, 0.0);
 
    std::cerr << "Beginning correlation function" << std::endl;
 
    #pragma omp parallel for num_threads(threads)
    for(size_t m = 0; m < eprops.size(); ++m){ 	// dipole, time = 0
        for(size_t n = m; 						//polarizability, time = 't'
 	      n < std::min( eprops.size(), m + int_time_ind); ++ n){ 
 
            int time = n - m;
 
            double* dip_k_0  = eprops[m].dip_k;
            double* dip2_k_0  = eprops[m].dip2_k;
            double* pol_ij_t = eprops[n].pol_ij;
            double* pol2_ij_t = eprops[n].pol2_ij;

            for(size_t i = 0; i < connects[m].central.size(); ++i){
                const size_t i_ind = connects[m].central[i];
                const double g_i0 = eprops[m].gscale[i_ind];
 
            #ifdef DO_AC // autocorrelation, j == i
                 {
                 #ifdef TYPE2
     	             const double g_it = eprops[n].gscale[i_ind];
                     const double g3 = g_i0*g_it*g_it;
                 #else
 	                 const double g3 = g_i0*g_i0*g_i0;
                 #endif
 
                 {
                 #pragma omp atomic 
                 corr[time] += g3*0.5*(dip_k_0[i_ind]*pol_ij_t[i_ind]+dip2_k_0[i_ind]*pol2_ij_t[i_ind]);
                 }
 
                 }
            #endif
 
            #ifdef DO_CC // truncated cross-correlation, j =/= i
                for(size_t j = 0; j < connects[m].partners[i].size(); ++j){
                    const size_t j_ind = connects[m].partners[i][j];
                    const double g_j0 = eprops[n].gscale[j_ind];
                    {
                    #pragma omp atomic
                    corr[time] += g_i0*g_j0*g_j0*0.5*(dip_k_0[i_ind]*pol_ij_t[j_ind]+dip2_k_0[i_ind]*pol2_ij_t[j_ind]);
                    }
                 }
            #endif
            }
 
 	   #pragma omp atomic
 	   ++ntime[time];
        }
    }
 
    std::cout << std::setprecision(10) << std::scientific;
 
    for(int i = 0; i < int_time_ind; ++i)
        std::cout << std::setw(20) << dt*i 
        << std::setw(20) << corr[i]/((double)ntime[i])
        << std::setw(20) << corr[i]
        << std::setw(20) << ((double)ntime[i])
        << std::endl;
}

