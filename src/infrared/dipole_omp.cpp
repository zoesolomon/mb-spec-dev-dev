// #ifdef HAVE_CONFIG_H
// #include "config.h"
// #endif // HAVE_CONFIG_H

// #include <omp.h>

// #include <cmath>
// #include <cassert>
// #include <cstdlib>

// #include <iomanip>
// #include <iostream>
// #include <fstream>
// #include <sstream>
// #include <string>
// #include <utility>
// #include <algorithm>

// #include "read-lammpstrj.h"
// #include "constants.h"

// #include "ttm4-dip.h"
// #include "ltp2011.h"
// #include "ttm4-periodic-es.h"
// #include "system.h"

// #include "dip_poly2b5.h"

// ////////////////////////////////////////////////////////////////////////////////

// namespace {

// ////////////////////////////////////////////////////////////////////////////////


// namespace t = ttm;

// typedef ttm::ttm4_dip ind_model;
// static h2o::ltp2011 dip1b_model;

// //----------------------------------------------------------------------------//

// inline double distance3(size_t imcon, const double* cell, const double* rcell,
// 		       const double* r1, const double* r2)
// {
//     assert(imcon == 0 || imcon == 1 || imcon == 2 || imcon == 3);

//     double dx = r1[0] - r2[0];
//     double dy = r1[1] - r2[1];
//     double dz = r1[2] - r2[2];

//     if (imcon == 1 || imcon == 2) {

//         dx -= cell[0]*std::floor(dx/cell[0] + 0.5);
//         dy -= cell[4]*std::floor(dy/cell[4] + 0.5);
//         dz -= cell[8]*std::floor(dz/cell[8] + 0.5);

//     } else if (imcon == 3) {

//         const double ssx = rcell[0]*dx + rcell[3]*dy + rcell[6]*dz;
//         const double ssy = rcell[1]*dx + rcell[4]*dy + rcell[7]*dz;
//         const double ssz = rcell[2]*dx + rcell[5]*dy + rcell[8]*dz;

//         const double xss = ssx - std::floor(ssx + 0.5);
//         const double yss = ssy - std::floor(ssy + 0.5);
//         const double zss = ssz - std::floor(ssz + 0.5);

//         dx = cell[0]*xss + cell[3]*yss + cell[6]*zss;
//         dy = cell[1]*xss + cell[4]*yss + cell[7]*zss;
//         dz = cell[2]*xss + cell[5]*yss + cell[8]*zss;

//     }

//     return std::sqrt(dx*dx + dy*dy + dz*dz);
// }

// //----------------------------------------------------------------------------//

// inline void image(const size_t imcon, const double* cell, const double* rcell,
//            const double* r1, const double* r2, double* crd_out)
// {
//     assert(imcon == 0 || imcon == 1 || imcon == 2 || imcon == 3);

//     double dx = r1[0] - r2[0];
//     double dy = r1[1] - r2[1];
//     double dz = r1[2] - r2[2];

//     if (imcon == 1 || imcon == 2) {

//        double facx = std::floor(dx/cell[0] + 0.5);
//        double facy = std::floor(dy/cell[4] + 0.5);
//        double facz = std::floor(dz/cell[8] + 0.5);

//        crd_out[0] = r2[0] + cell[0]*facx;
//        crd_out[1] = r2[1] + cell[4]*facy;
//        crd_out[2] = r2[2] + cell[8]*facz;

//     } else if (imcon == 3) {

//        const double ssx = rcell[0]*dx + rcell[3]*dy + rcell[6]*dz;
//        const double ssy = rcell[1]*dx + rcell[4]*dy + rcell[7]*dz;
//        const double ssz = rcell[2]*dx + rcell[5]*dy + rcell[8]*dz;

//        const double xss = ssx - std::floor(ssx + 0.5);
//        const double yss = ssy - std::floor(ssy + 0.5);
//        const double zss = ssz - std::floor(ssz + 0.5);

//        dx = cell[0]*xss + cell[3]*yss + cell[6]*zss;
//        dy = cell[1]*xss + cell[4]*yss + cell[7]*zss;
//        dz = cell[2]*xss + cell[5]*yss + cell[8]*zss;

//        crd_out[0] = r1[0] - dx;
//        crd_out[1] = r1[1] - dy;
//        crd_out[2] = r1[2] - dz;

//     }
// }

// //----------------------------------------------------------------------------//

// inline void image_dimer(const size_t imcon, const double* cell, 
//                 const double* rcell, const double* r1, const double* r2,
// 		      	double* dimer)
// {
//     std::fill(dimer, dimer + 18, 0.0);

//     // take r1 oxygen as central atom, image its two hydrogens
     
//     for(size_t i = 0; i < 3; ++i)
// 	    dimer[i] = r1[i];
    
//     for(size_t a = 1; a < 3; ++a)
// 	    image(imcon, cell, rcell, r1, r1 + 3*a, dimer + 3*a);

//     // Now image the oxygen of the second molecule

//     image(imcon, cell, rcell, r1, r2, dimer + 9);

//     // Finally, reattach hydrogens of molec two to oxygen of molec 2

//     for(size_t a = 1; a < 3; ++a)
// 	    image(imcon, cell, rcell, dimer + 9, r2 + 3*a, dimer + 9 + 3*a);
// }

// //----------------------------------------------------------------------------//

// double f_switch2(const double& r, const double r2i, const double r2f)
// {
//     if (r > r2f)
// 	    return 0.0;
//     else if (r > r2i) {
// 	    const double x = (r - r2i)/(r2f - r2i);
// 	    return 1 + x*x*(2*x - 3);
//     } else
// 	    return 1.0;
// }

// ////////////////////////////////////////////////////////////////////////////////


// } // namespace

// ////////////////////////////////////////////////////////////////////////////////

// int main(int argc, char** argv)
// {

//   if (argc < 3) {
//     std::cerr << "usage: " << argv[0] << " <file_name> <mb-spec.json>" << std::endl;
//     return EXIT_FAILURE;
//   }

//   System sys;
//   sys.SetUpFromJson("INFRARED", argv[argc - 1]);

//   std::vector<double> box(9);

//   box = sys.GetBoxVector();

//   const double r2i = sys.GetCutoff1();
//   const double r2f = sys.GetCutoff2();


//   // Open files for Writing

//   std::cout.setf(std::ios_base::showpoint);
//   std::cout.precision(9);

//   size_t lineno(0);
//   size_t frameno(0);

//   ind_model pot_ind;

//   double* dimer = new double[18]; // FOR MB

//   // begin the parallelization over frames

//   std::string fname = argv[1];
//   std::ifstream ifs_cmd(fname);
//   size_t pos = fname.find('.');
//   std::string format  = fname.erase(0, pos + 1);
//   std::transform(format.begin(), format.end(), format.begin(), ::tolower);

//   size_t imcon(0);
//   {
//     imcon = sys.GetImcon();
//   }

//   if (world_rank == 0) {
//     if (imcon == 0){
//       imcon = 2;
//     } else if (imcon < 0 || imcon > 3) {
//       std::cerr << "Cannot calculate system of imcon = " << imcon << std::endl;
//       return EXIT_FAILURE;
//     } else if (imcon > 0 && imcon < 3) {
//       if (imcon == 1) {
//         std::cerr << "imcon = " << imcon << " (cubic)" << std::endl;
//       } else if (imcon == 2) {
//         std::cerr << "imcon = " << imcon << " (orthorhombic)" << std::endl;
//       } else if (imcon == 3) {
//         std::cerr << "imcon = " << imcon << " (monoclinic)" << std::endl;
//       }
//     } else {
//       std::cerr << "imcon value of " << imcon << " not a valid integer" << std::endl;
//       return EXIT_FAILURE;
//     }
//   }

//   double mu(0);
//   mu = sys.GetMu();

//   double dump(0);
//   dump = sys.GetDumpFrequency();

//   if (world_rank == 0) {
//     if (mu != 0) {
//       std::cerr << "Using LTP 1B DMS" << std::endl;
//     } else {
//       std::cerr << "Reading from " << fname << std::endl;
//     }
//   }


//   // Toggles names: 1B, NB (1BNB), MB to lowercase
//   std::string mb = sys.GetMb();
//   std::transform(mb.begin(), mb.end(), mb.begin(), ::tolower);

//   // duh can't read dipoles if there r no dipoles
//   if (world_rank == 0 && mu == 0 && format.compare("xyz") == 0) {
//     std::cerr << ".lammpstrj format required to read dipoles" << std::endl;
//     return EXIT_FAILURE;
//   }


//   while(!ifs_cmd.eof() ) {
//     std::vector<opt::lammpstrj_frame> frames;

//     // reads in the number of frames equal to num of processors for MPI
//     while((int) frames.size() < world_size && !ifs_cmd.eof()){
//       opt::lammpstrj_frame this_frame;
//       if(format.compare("xyz") == 0) {
//         opt::read_xyz_frame(lineno, imcon, ifs_cmd, this_frame, dump * frameno, box);
//       } else {
//         opt::read_lammpstrj_frame(lineno, imcon, ifs_cmd, this_frame, dump * frameno);
//       }

//       frameno++;

//       if (ifs_cmd.eof()) break;
//       if (this_frame.natm == 0) break;

//       frames.push_back(this_frame);

//     }

//     int n_frames = frames.size();

//     if (n_frames == 0) break;


//     opt::lammpstrj_frame frame;

//     double total_dip[3];
//     double dip_mag = 0.0;


//     if (n_frames > world_rank) {

//       frame = frames[world_rank];


//       const size_t nmolec = frame.natm/3;
//       std::fill (total_dip, total_dip + 3, 0.0);

//       double dip2b_corr[3];
//       double dip2b_ind[3];
//       double dip_ind[3];
//       std::fill(dip_ind, dip_ind + 3, 0.0);

//       // ONLY FOR MB, two body
//       // only for waters
//       //if (mb.compare("1b") != 0) {
//       if (mb.compare("mb") == 0) {
//         if (frame.natm % 3 == 0 && nmolec == frame.molecid[frame.natm - 1]) {
//           std::fill(dip2b_corr, dip2b_corr + 3, 0.0);
//           std::fill(dip2b_ind, dip2b_ind + 3, 0.0);

//           for (size_t i = 0; i < nmolec; ++i) {
//             for (size_t j = i+1; j < nmolec; ++j) {
//               const size_t i9 = 9*i;
//               const size_t j9 = 9*j;
//               double rOO = distance3(imcon, frame.cell, frame.rcell, 
//                                          frame.xyz + i9, frame.xyz + j9);

//               if (rOO < r2f) {
//                 const double s = f_switch2(rOO, r2i, r2f);
//                 image_dimer(imcon, frame.cell, frame.rcell, 
//                                frame.xyz + i9, frame.xyz + j9, dimer);
 
//                 // Get corr two-body dipole
//                 {
//                   double tmp_dip2b[3];
//                   h2o::dip_poly2b5::eval(dimer, dimer + 9, tmp_dip2b);
  
//                   for (size_t b = 0; b < 3; ++b) {
//                     dip2b_corr[b] += s*tmp_dip2b[b];
//                   }
//                 }
  
//                 // Get TTM4 two-body dipole
//                 {
//                   double tmp_dip2b[3];
//                   pot_ind(2, dimer, tmp_dip2b);
 
//                   for(size_t b = 0; b < 3; ++b) {
//                     //dip2b_ind[b] += s*tmp_dip2b[b];
//                     dip2b_ind[b] += (1 - s)*tmp_dip2b[b];
//                     dip_ind[b] += s*tmp_dip2b[b];
//                   }
//                 }
//               }
//             }
//           }
//         } else {
//           if (world_rank == 0) {
//             std::cerr << "'mb' feature can only be used for water" << std::endl;
//           }
//           return EXIT_FAILURE;
//         }
//       }

//       for(size_t i = 0; i < frame.natm; ++i){
//         size_t i3 = 3*i;
//         size_t i4 = 4*i;
//         size_t i9 = 9*i;

//         double dip_mol[3];

//         if (mu == 0) {
//           for (size_t j = 0; j < 3; ++j) {
//             total_dip[j] += frame.dipmol[i4 + j]; // for all
//             if (mb.compare("1b") != 0) {  // for 1bnb + mb
//               total_dip[j] += frame.dipind[i4 + j];
//             }
//           }
//         } else if (mu !=0 || (total_dip[0] == 0 && total_dip[1] == 0 && total_dip[2] == 0)) { // if calc mb-mu or if dipoles read only zeros
//           if (frame.natm % 3 == 0 && nmolec == frame.molecid[frame.natm - 1]) { // check to make sure we have only water
//             if (i % 3 == 0) { // only want to calc for each molecule not each atom
//               double monomer[9];
//               std::copy(frame.xyz + i3, frame.xyz + i3 + 3, monomer);
//               for(size_t a = 1; a < 3; ++a) {
//                 image(imcon, frame.cell, frame.rcell, frame.xyz + i3,
//                 frame.xyz + i3 + 3*a, monomer + 3*a);
//               }    
//               dip1b_model.cartesian(1, monomer, dip_mol);
//               //std::cerr << dip_mol[0] << "\t" << dip_mol[1] << "\t" << dip_mol[2] << std::endl;
//               for (size_t t = 0; t < 3; t++) {
//                 total_dip[t] += dip_mol[t];
//               }
//             }
//           } else {
//             if (world_rank == 0) {
//               std::cerr << "MB-μ can only be calculated for water" << std::endl;
//             }
//             return EXIT_FAILURE;
//           }
//         }
//       }

// // Substitute two-body dipole contribution from induction with corr
// // for MB
//       if (mb.compare("mb") == 0) {
//         for (size_t j = 0; j < 3; ++j) {
//         //  total_dip[j] += dip2b_corr[j];
//          // total_dip[j] += (dip2b_corr[j] + dip2b_ind[j]);
//           total_dip[j] += (dip2b_corr[j] - dip2b_ind[j]);
//         }
//       } //else if (mb.compare("nb") == 0 && mu != 0) {
//         //for (size_t j = 0; j < 3; ++j) {
//         //  total_dip[j] += dip_ind[j];
//         //}
//       //}

//       dip_mag = std::sqrt(total_dip[0]*total_dip[0]
//                   + total_dip[1]*total_dip[1]
//                   + total_dip[2]*total_dip[2]);


//     }

//     double all_frame_time[1 * world_size];
//     MPI_Gather(&frame.time, 1, MPI_DOUBLE, &all_frame_time, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

//     double all_frame_dip[3 * world_size];
//     MPI_Gather(&total_dip, 3, MPI_DOUBLE, &all_frame_dip, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);

//     double all_frame_dip_mag[1 * world_size];
//     MPI_Gather(&dip_mag, 1, MPI_DOUBLE, &all_frame_dip_mag, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);


//     if (world_rank == 0){
//       for (int i = 0; i < n_frames; ++i) 
//         std::cout << all_frame_time[i] << ' '
//                   << std::setw(14) << all_frame_dip[i * 3 + 0] << ' ' 
//                   << std::setw(14) << all_frame_dip[i * 3 + 1] << ' ' 
//                   << std::setw(14) << all_frame_dip[i * 3 + 2] << ' ' 
//                   << std::setw(14) << all_frame_dip_mag[i] << std::endl;
//     }

//   }

//   delete[] dimer;
//   MPI_Finalize();
// }

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <omp.h>

#include <cmath>
#include <cassert>
#include <cstdlib>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <utility>
#include <algorithm>
#include <vector>
#include <array>

#include "read-lammpstrj.h"
#include "constants.h"

#include "ttm4-dip.h"
#include "ltp2011.h"
#include "ttm4-periodic-es.h"
#include "system.h"

#include "dip_poly2b5.h"

////////////////////////////////////////////////////////////////////////////////

namespace {

////////////////////////////////////////////////////////////////////////////////

namespace t = ttm;

typedef ttm::ttm4_dip ind_model;
static h2o::ltp2011 dip1b_model;

//----------------------------------------------------------------------------//

inline double distance3(size_t imcon, const double* cell,
                        const double* rcell,
                        const double* r1,
                        const double* r2)
{
    assert(imcon == 0 || imcon == 1 ||
           imcon == 2 || imcon == 3);

    double dx = r1[0] - r2[0];
    double dy = r1[1] - r2[1];
    double dz = r1[2] - r2[2];

    if (imcon == 1 || imcon == 2) {

        dx -= cell[0]*std::floor(dx/cell[0] + 0.5);
        dy -= cell[4]*std::floor(dy/cell[4] + 0.5);
        dz -= cell[8]*std::floor(dz/cell[8] + 0.5);

    } else if (imcon == 3) {

        const double ssx =
            rcell[0]*dx + rcell[3]*dy + rcell[6]*dz;

        const double ssy =
            rcell[1]*dx + rcell[4]*dy + rcell[7]*dz;

        const double ssz =
            rcell[2]*dx + rcell[5]*dy + rcell[8]*dz;

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

inline void image(const size_t imcon,
                  const double* cell,
                  const double* rcell,
                  const double* r1,
                  const double* r2,
                  double* crd_out)
{
    assert(imcon == 0 || imcon == 1 ||
           imcon == 2 || imcon == 3);

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

       const double ssx =
           rcell[0]*dx + rcell[3]*dy + rcell[6]*dz;

       const double ssy =
           rcell[1]*dx + rcell[4]*dy + rcell[7]*dz;

       const double ssz =
           rcell[2]*dx + rcell[5]*dy + rcell[8]*dz;

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

inline void image_dimer(const size_t imcon,
                        const double* cell,
                        const double* rcell,
                        const double* r1,
                        const double* r2,
                        double* dimer)
{
    std::fill(dimer, dimer + 18, 0.0);

    for(size_t i = 0; i < 3; ++i)
        dimer[i] = r1[i];

    for(size_t a = 1; a < 3; ++a)
        image(imcon, cell, rcell,
              r1,
              r1 + 3*a,
              dimer + 3*a);

    image(imcon, cell, rcell,
          r1,
          r2,
          dimer + 9);

    for(size_t a = 1; a < 3; ++a)
        image(imcon, cell, rcell,
              dimer + 9,
              r2 + 3*a,
              dimer + 9 + 3*a);
}

//----------------------------------------------------------------------------//

double f_switch2(const double& r,
                 const double r2i,
                 const double r2f)
{
    if (r > r2f)
        return 0.0;

    else if (r > r2i) {

        const double x =
            (r - r2i)/(r2f - r2i);

        return 1 + x*x*(2*x - 3);

    } else {
        return 1.0;
    }
}

////////////////////////////////////////////////////////////////////////////////

} // namespace

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
    omp_set_nested(1);

    if (argc < 3) {

        std::cerr
            << "usage: "
            << argv[0]
            << " <file_name> <mb-spec.json>"
            << std::endl;

        return EXIT_FAILURE;
    }

    System sys;

    sys.SetUpFromJson("INFRARED",
                      argv[argc - 1]);

    std::vector<double> box(9);

    box = sys.GetBoxVector();

    const double r2i = sys.GetCutoff1();
    const double r2f = sys.GetCutoff2();

    std::cout.setf(std::ios_base::showpoint);
    std::cout.precision(9);

    size_t lineno(0);
    size_t frameno(0);

    ind_model pot_ind;

    std::string fname = argv[1];

    std::ifstream ifs_cmd(fname);

    size_t pos = fname.find('.');

    std::string format =
        fname.erase(0, pos + 1);

    std::transform(format.begin(),
                   format.end(),
                   format.begin(),
                   ::tolower);

    size_t imcon = sys.GetImcon();

    if (imcon == 0) {

        imcon = 2;

    } else if (imcon < 0 || imcon > 3) {

        std::cerr
            << "Cannot calculate system of imcon = "
            << imcon
            << std::endl;

        return EXIT_FAILURE;

    } else if (imcon == 1) {

        std::cerr
            << "imcon = 1 (cubic)"
            << std::endl;

    } else if (imcon == 2) {

        std::cerr
            << "imcon = 2 (orthorhombic)"
            << std::endl;

    } else if (imcon == 3) {

        std::cerr
            << "imcon = 3 (monoclinic)"
            << std::endl;
    }

    double mu = sys.GetMu();

    double dump = sys.GetDumpFrequency();

    if (mu != 0) {

        std::cerr
            << "Using LTP 1B DMS"
            << std::endl;

    } else {

        std::cerr
            << "Reading from "
            << fname
            << std::endl;
    }

    std::string mb = sys.GetMb();

    std::transform(mb.begin(),
                   mb.end(),
                   mb.begin(),
                   ::tolower);

    if (mu == 0 &&
        format.compare("xyz") == 0) {

        std::cerr
            << ".lammpstrj format required to read dipoles"
            << std::endl;

        return EXIT_FAILURE;
    }

    while(!ifs_cmd.eof()) {

        std::vector<opt::lammpstrj_frame> frames;

        while((int)frames.size() < 32 &&
              !ifs_cmd.eof()) {

            opt::lammpstrj_frame this_frame;

            if(format.compare("xyz") == 0) {

                opt::read_xyz_frame(
                    lineno,
                    imcon,
                    ifs_cmd,
                    this_frame,
                    dump * frameno,
                    box);

            } else {

                opt::read_lammpstrj_frame(
                    lineno,
                    imcon,
                    ifs_cmd,
                    this_frame,
                    dump * frameno);
            }

            frameno++;

            if (ifs_cmd.eof()) break;
            if (this_frame.natm == 0) break;

            frames.push_back(this_frame);
        }

        int n_frames = frames.size();

        if (n_frames == 0)
            break;

        std::vector<std::array<double,3>>
            all_frame_dip(n_frames);

        std::vector<double>
            all_frame_time(n_frames);

        std::vector<double>
            all_frame_dip_mag(n_frames);

        //======================================================
        // FRAME PARALLELISM
        //======================================================

        #pragma omp parallel for schedule(dynamic)
        for (int f = 0; f < n_frames; ++f) {

            opt::lammpstrj_frame& frame =
                frames[f];

            double total_dip[3] =
                {0.0, 0.0, 0.0};

            double dip_mag = 0.0;

            const size_t nmolec =
                frame.natm/3;

            double dip2b_corr[3] =
                {0.0, 0.0, 0.0};

            double dip2b_ind[3] =
                {0.0, 0.0, 0.0};

            double dip_ind[3] =
                {0.0, 0.0, 0.0};

            //==================================================
            // MB TWO-BODY
            //==================================================

            if (mb.compare("mb") == 0) {

                if (frame.natm % 3 == 0 &&
                    nmolec ==
                    frame.molecid[frame.natm - 1]) {

                    //==========================================
                    // PAIR PARALLELISM
                    //==========================================

                    #pragma omp parallel
                    {
                        double local_corr[3] =
                            {0.0, 0.0, 0.0};

                        double local_ind2b[3] =
                            {0.0, 0.0, 0.0};

                        double local_ind[3] =
                            {0.0, 0.0, 0.0};

                        double dimer[18];

                        #pragma omp for schedule(dynamic)
                        for (size_t i = 0;
                             i < nmolec;
                             ++i) {

                            for (size_t j = i + 1;
                                 j < nmolec;
                                 ++j) {

                                const size_t i9 = 9*i;
                                const size_t j9 = 9*j;

                                double rOO =
                                    distance3(
                                        imcon,
                                        frame.cell,
                                        frame.rcell,
                                        frame.xyz + i9,
                                        frame.xyz + j9);

                                if (rOO < r2f) {

                                    const double s =
                                        f_switch2(
                                            rOO,
                                            r2i,
                                            r2f);

                                    image_dimer(
                                        imcon,
                                        frame.cell,
                                        frame.rcell,
                                        frame.xyz + i9,
                                        frame.xyz + j9,
                                        dimer);

                                    // corr dipole
                                    {
                                        double tmp[3];

                                        h2o::dip_poly2b5::eval(
                                            dimer,
                                            dimer + 9,
                                            tmp);

                                        for (size_t b = 0;
                                             b < 3;
                                             ++b) {

                                            local_corr[b] +=
                                                s*tmp[b];
                                        }
                                    }

                                    // induction dipole
                                    {
                                        double tmp[3];

                                        pot_ind(2,
                                                dimer,
                                                tmp);

                                        for(size_t b = 0;
                                            b < 3;
                                            ++b) {

                                            local_ind2b[b] +=
                                                (1 - s)*tmp[b];

                                            local_ind[b] +=
                                                s*tmp[b];
                                        }
                                    }
                                }
                            }
                        }

                        #pragma omp critical
                        {
                            for (size_t b = 0;
                                 b < 3;
                                 ++b) {

                                dip2b_corr[b] +=
                                    local_corr[b];

                                dip2b_ind[b] +=
                                    local_ind2b[b];

                                dip_ind[b] +=
                                    local_ind[b];
                            }
                        }
                    }

                } else {

                    std::cerr
                        << "'mb' feature can only be used for water"
                        << std::endl;

                    continue;
                }
            }

            //==================================================
            // MONOMER / INPUT DIPOLES
            //==================================================

            for(size_t i = 0;
                i < frame.natm;
                ++i){

                size_t i3 = 3*i;
                size_t i4 = 4*i;

                double dip_mol[3];

                if (mu == 0) {

                    for (size_t j = 0;
                         j < 3;
                         ++j) {

                        total_dip[j] +=
                            frame.dipmol[i4 + j];

                        if (mb.compare("1b") != 0) {

                            total_dip[j] +=
                                frame.dipind[i4 + j];
                        }
                    }

                } else {

                    if (frame.natm % 3 == 0 &&
                        nmolec ==
                        frame.molecid[frame.natm - 1]) {

                        if (i % 3 == 0) {

                            double monomer[9];

                            std::copy(
                                frame.xyz + i3,
                                frame.xyz + i3 + 3,
                                monomer);

                            for(size_t a = 1;
                                a < 3;
                                ++a) {

                                image(
                                    imcon,
                                    frame.cell,
                                    frame.rcell,
                                    frame.xyz + i3,
                                    frame.xyz + i3 + 3*a,
                                    monomer + 3*a);
                            }

                            dip1b_model.cartesian(
                                1,
                                monomer,
                                dip_mol);

                            for (size_t t = 0;
                                 t < 3;
                                 ++t) {

                                total_dip[t] +=
                                    dip_mol[t];
                            }
                        }
                    }
                }
            }

            //==================================================
            // MB correction
            //==================================================

            if (mb.compare("mb") == 0) {

                for (size_t j = 0;
                     j < 3;
                     ++j) {

                    total_dip[j] +=
                        (dip2b_corr[j]
                        - dip2b_ind[j]);
                }
            }

            dip_mag =
                std::sqrt(
                    total_dip[0]*total_dip[0]
                  + total_dip[1]*total_dip[1]
                  + total_dip[2]*total_dip[2]);

            all_frame_time[f] =
                frame.time;

            all_frame_dip_mag[f] =
                dip_mag;

            for (size_t j = 0;
                 j < 3;
                 ++j) {

                all_frame_dip[f][j] =
                    total_dip[j];
            }
        }

        //======================================================
        // OUTPUT
        //======================================================

        for (int i = 0;
             i < n_frames;
             ++i) {

            std::cout
                << all_frame_time[i] << ' '
                << std::setw(14)
                << all_frame_dip[i][0] << ' '
                << std::setw(14)
                << all_frame_dip[i][1] << ' '
                << std::setw(14)
                << all_frame_dip[i][2] << ' '
                << std::setw(14)
                << all_frame_dip_mag[i]
                << std::endl;
        }
    }

    return 0;
}