#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <mpi.h>

#include <cmath>
#include <cassert>
#include <cstdlib>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>

#include "read-lammpstrj.h"
#include "constants.h"
#include "system.h"

#include "avila.h"

#include "alpha-el-poly2b3.h"

#include "alpha-ttm4.h"

#include "ttm4-periodic-es.h"

////////////////////////////////////////////////////////////////////////////////

namespace {

////////////////////////////////////////////////////////////////////////////////


static x2o::avila pol_model;

constexpr double r2i = 5.5;
constexpr double r2f = 7.5;

//----------------------------------------------------------------------------//

void print_dimer(double* crd, double time)
{
    std::cerr << "6\n" << time << "\n";
    for(size_t n = 0; n < 2; ++n){
	for(size_t i = 0; i < 3; ++i){
	    if(i%3 == 0)
		std::cerr << "O(frag=" << n << ")  ";
	    else
		std::cerr << "H(frag=" << n << ")  ";
	    for(size_t j = 0; j < 3; ++j)
		std::cerr << std::setw(23) << crd[9*n+3*i+j];
	    std::cerr << std::endl;
	}
    }
}

//----------------------------------------------------------------------------//

inline double distance3(size_t imcon, const double* cell, const double* rcell,
		       const double* r1, const double* r2)
{
    assert(imcon == 0 || imcon == 1 || imcon == 2 || imcon == 3);

    double dx = r1[0] - r2[0];
    double dy = r1[1] - r2[1];
    double dz = r1[2] - r2[2];

    if (imcon == 1 || imcon == 2) {

        dx -= cell[0]*std::floor(dx/cell[0] + 0.5);
        dy -= cell[4]*std::floor(dy/cell[4] + 0.5);
        dz -= cell[8]*std::floor(dz/cell[8] + 0.5);

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

    }

    return std::sqrt(dx*dx + dy*dy + dz*dz);

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

       double facx = std::floor(dx/cell[0] + 0.5);
       double facy = std::floor(dy/cell[4] + 0.5);
       double facz = std::floor(dz/cell[8] + 0.5);

       crd_out[0] = r2[0] + cell[0]*facx;
       crd_out[1] = r2[1] + cell[4]*facy;
       crd_out[2] = r2[2] + cell[8]*facz;

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

inline void image_dimer(const size_t imcon, const double* cell, 
                const double* rcell, const double* r1, const double* r2,
		      	double* dimer)
{
    std::fill(dimer, dimer + 18, 0.0);

    // take r1 oxygen as central atom, image its two hydrogens
     
    for(size_t i = 0; i < 3; ++i)
	    dimer[i] = r1[i];
    
    for(size_t a = 1; a < 3; ++a)
	    image(imcon, cell, rcell, r1, r1 + 3*a, dimer + 3*a);

    // Now image the oxygen of the second molecule

    image(imcon, cell, rcell, r1, r2, dimer + 9);

    // Finally, reattach hydrogens of molec two to oxygen of molec 2

    for(size_t a = 1; a < 3; ++a)
	    image(imcon, cell, rcell, dimer + 9, r2 + 3*a, dimer + 9 + 3*a);
}

//----------------------------------------------------------------------------//

double f_switch2(const double& r)
{
    if (r > r2f)
	return 0.0;
    else if (r > r2i) {
	const double x = (r - r2i)/(r2f - r2i);
	return 1 + x*x*(2*x - 3);
    } else
	return 1.0;
}

////////////////////////////////////////////////////////////////////////////////

  static int world_rank(0);
  static int world_size(0);

} // namespace

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    if (argc < 3) {
        std::cerr << "usage: mpirun -np 4 " << argv[0] << " <file_name> <mb-spec.json>" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    System sys;
    sys.SetUpFromJson("RAMAN", argv[argc - 1]);

    std::vector<double> box(9);

    box = sys.GetBoxVector();

    const double r2i = sys.GetCutoff1();
    const double r2f = sys.GetCutoff2();


    std::string fname = argv[1];
    std::ifstream ifs_cmd(fname);
    size_t pos = fname.rfind('.');
    std::string format  = fname.erase(0, pos + 1);
    std::transform(format.begin(), format.end(), format.begin(), ::tolower);


    // Open files for Writing

    std::cout.setf(std::ios_base::showpoint);
    std::cout.precision(9);

    std::cerr.setf(std::ios_base::scientific);
    std::cerr.precision(12);

    size_t lineno(0);
    size_t frameno(0);
    double* dimer = new double[18];
    size_t num_frame(0);


    size_t imcon(0);
    {
        imcon = sys.GetImcon();
    }

    if (world_rank == 0) {
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
    }

    double dump(0);
    dump = sys.GetDumpFrequency();

    // Toggles names: 1B, NB (1BNB), MB to lowercase
    std::string mb = sys.GetMb();
    std::transform(mb.begin(), mb.end(), mb.begin(), ::tolower);


    while (!ifs_cmd.eof()) {

        std::vector<opt::lammpstrj_frame> frames;

        while ((int) frames.size() < world_size && !ifs_cmd.eof()) {
            opt::lammpstrj_frame this_frame;
            if(format.compare("xyz") == 0) {
                opt::read_xyz_frame(lineno, imcon, ifs_cmd, this_frame, dump * frameno, box);
            } else {
                opt::read_lammpstrj_frame(lineno, imcon, ifs_cmd, this_frame, dump * frameno);
            }
            if (this_frame.natm % 3 != 0 || this_frame.molecid[this_frame.natm - 1] != (this_frame.natm / 3)) {
                if (world_rank == 0) {
                    std::cerr << "MB-α can only be calculated for water" << std::endl;
                }
                return EXIT_FAILURE;
            }
            if(ifs_cmd.eof()) break;
	    if(this_frame.natm == 0)
	        break;
            frames.push_back(this_frame);
        }

        int n_frames = frames.size();
        if (n_frames == 0) break;
    
        opt::lammpstrj_frame frame;
        double total_pol[9];

        if (n_frames > world_rank) {

            frame = frames[world_rank];
            const size_t nmolec = frame.natm/3;

            std::fill(total_pol, total_pol + 9, 0.0);

            double pol1b[9];
            std::fill(pol1b, pol1b + 9, 0.0);

            double pol2b_corr[9];
            std::fill(pol2b_corr, pol2b_corr + 9, 0.0);

            double polMb[9];
            std::fill(polMb, polMb + 9, 0.0);

            for (size_t i = 0; i < nmolec; ++i) {
                size_t i9 = 9*i;

                /////
                // 1B contribution: avila
                /////

                // FOR 1B, 1BNB, MB
                {
      	            double monomer[9];
      	            std::copy(frame.xyz + i9, frame.xyz + i9 + 3, monomer);
      	            for(size_t a = 1; a < 3; ++a)
      	                image(imcon, frame.cell, frame.rcell, frame.xyz + i9,
      		                  frame.xyz + i9 + 3*a, monomer + 3*a);

      	            // Avila 1B polarizability
      	            double pol_mol[9];
      	            std::fill(pol_mol, pol_mol + 9, 0.0);
      	            pol_model.alpha(monomer, pol_mol);
      	            for (size_t a = 0; a < 9; ++a)
      	                pol1b[a] += pol_mol[a];

      	            // TTM4 1B polarizability
      	            if (mb.compare("1b") != 0 ) {
                        double pol_ind_1b[9];
      	                std::fill(pol_ind_1b, pol_ind_1b + 9, 0.0);
      	                ttm::alpha_tot(1, monomer, pol_ind_1b);
      	                for(size_t a = 0; a < 9; ++a)
      	                    polMb[a] -= pol_ind_1b[a];
                    }
                }
                
                /////
                // 2B contribution: s*[ \alpha^{corr} ]
                /////

                // FOR MB
                if (mb.compare("mb") == 0) {
                    for (size_t j = i+1; j < nmolec; ++j) {
      	                const size_t j9 = 9*j;

    	                double rOO = distance3(imcon, frame.cell, frame.rcell, 
      		 	                               frame.xyz + i9, frame.xyz + j9);

      	                if (rOO < r2f) {

      	                    const double s = f_switch2(rOO);
      	                    image_dimer(imcon, frame.cell, frame.rcell, 
      			                frame.xyz + i9, frame.xyz + j9, dimer);

      	                    double corr[9];
      	                    std::fill(corr, corr + 9, 0.0);
      	                    h2o::alpha_el_poly2b3::eval(dimer, dimer + 9, corr);

      	                    for (size_t b = 0; b < 9; ++b) {
      		                pol2b_corr[b] += s*corr[b];
      	                    }
      	                }    
                    }
                }
                ++frameno;

            }

            /////
            // MB contribution: alpha^{ttm4,NB}(1, \dots, N) 
            //                  - \sum_{i  }^N \alpha^{ttm4,1B}(i)
            //
            //    Note: alpha^{ttm4,NB} contains a monomer polarizability
            //          Therefore, the 1B polarizability
            //          must be removed to obtain the many-body polarizability
            //
            //          The 1B polarizability has already been
            //          removed above in the polMb term
            /////
 	
            // FOR 1BNB and MB, NOT 1B
            if (mb.compare("1b") != 0) {

	        double ind[9];
	        std::fill(ind, ind + 9, 0.0);

	        ttm::ttm4_periodic_es ind_model;
	        ind_model(nmolec, frame.xyz, frame.cell);
	        ind_model.ptensor(nmolec, ind);

	        for(size_t a = 0; a < 9; ++a)
		    polMb[a] += ind[a];
            }

	    for(size_t a = 0; a < 9; ++a)
	        total_pol[a] = pol1b[a] + pol2b_corr[a] + polMb[a];

        } else {
            opt::lammpstrj_frame frame = frames[0];
        }

        double all_frame_time[1 * world_size];
        MPI_Gather(&frame.time, 1, MPI_DOUBLE, &all_frame_time, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        double all_frame_pol[9 * world_size];
        MPI_Gather(&total_pol, 9, MPI_DOUBLE, &all_frame_pol, 9, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        if (world_rank == 0) {
            for (int i =0; i < n_frames; ++i)
    	    std::cout << all_frame_time[i] << ' ' 
    	              << std::setw(14) << all_frame_pol[i*9 + 0]   << ' ' 
    	              << std::setw(14) << all_frame_pol[i*9 + 1]   << ' ' 
    	              << std::setw(14) << all_frame_pol[i*9 + 2]   << ' '
                   // << std::setw(14) << all_frame_pol[i*9 + 3]   << ' '
    	              << std::setw(14) << all_frame_pol[i*9 + 4]   << ' '
    	              << std::setw(14) << all_frame_pol[i*9 + 5]   << ' '
                   // << std::setw(14) << all_frame_pol[i*9 + 6]   << ' '
                   // << std::setw(14) << all_frame_pol[i*9 + 7]   << ' '
    	              << std::setw(14) << all_frame_pol[i*9 + 8]   << std::endl;
        } 
    }

    delete[] dimer;
    MPI_Finalize();
}
