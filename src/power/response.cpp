#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <sstream>

#include "system.h"
#include "read-lammpstrj.h"
#include <omp.h>

namespace {

static size_t ntraj = 0;
static int nsamples = 0;

static double t0;
static double dt;

static std::vector<double> accum;
static std::vector<double> ntrajs;

const double debye = 3.33564e-30; // Coulomb*m


////////////////////////////////////////////////////////////////////////////////

// void do_autocorr(double* times, double* vels, const double int_max, double size, std::string fname)
// {

//     dt = times[1] - times[0];
//     int this_file_nsamples = size;
    
//     int int_time_ind = (int) (int_max / dt);

//     double corr_x[int_time_ind];
//     double corr_y[int_time_ind];
//     double corr_z[int_time_ind];
//     double ntime[int_time_ind];

//     std::fill(corr_x, corr_x + int_time_ind, 0.0);
//     std::fill(corr_y, corr_y + int_time_ind, 0.0);
//     std::fill(corr_z, corr_z + int_time_ind, 0.0);
//     std::fill(ntime, ntime + int_time_ind, 0.0);

//     for(int i = 0; i < this_file_nsamples; ++i){
// 	for(int j = i; j < std::min(this_file_nsamples, i + int_time_ind); ++j){
// 	    int time = (j-i);
//         int i3 = 3*i;
//         int j3 = 3*j;

// 	    corr_x[time] += vels[i3]*vels[j3];
// 	    corr_y[time] += vels[i3 + 1]*vels[j3 + 1];
// 	    corr_z[time] += vels[i3 + 2]*vels[j3 + 2];
// 	    ++ntime[time];
// 	}
//     }

//     for(int i = 0; i < int_time_ind; ++i){
// 	corr_x[i] /= ntime[i];
//     	corr_y[i] /= ntime[i];
//     	corr_z[i] /= ntime[i];
//     }

//     std::cerr << "* loaded " << int_time_ind
//               << " samples from '" << fname << "'"
//               << std::endl;

//     if (nsamples < int_time_ind){

// 	for (int k = 0; k < nsamples; ++k) {
// 	    const double autocorr =
// 		(corr_x[k] + corr_y[k] + corr_z[k]);

// 	    accum[k] += autocorr;
// 	    ++ntrajs[k];
// 	}

// 	for (int k = nsamples; k < int_time_ind; ++k) {
// 	    const double autocorr =
// 		(corr_x[k] + corr_y[k] + corr_z[k]);

// 	    accum.push_back(autocorr);
// 	    ntrajs.push_back(1.0);
// 	}

// 	nsamples = int_time_ind;

//     }else{
// 	for (int k = 0; k < int_time_ind; ++k) {
// 	    const double autocorr =
// 		(corr_x[k] + corr_y[k] + corr_z[k]);
	    
// 	    accum[k] += autocorr;
// 	    ++ntrajs[k];
// 	}

//     }

//     ++ntraj;
// }

// void do_autocorr(std::vector<double> dat, const double int_max, int size, std::string fname)
// void do_autocorr(std::vector<double> dat, const double int_max, std::string fname, double dump)
// {
//     double dt = dump;
//     std::cerr << dt << std::endl;
//     int int_time_ind = (int)(int_max / dt);
//     int size = dat.size();

//     std::vector<double> corr_x(int_time_ind, 0.0);
//     std::vector<double> corr_y(int_time_ind, 0.0);
//     std::vector<double> corr_z(int_time_ind, 0.0);
//     std::vector<double> ntime(int_time_ind, 0.0);


//     for (int i = 0; i < size; ++i) {
//         for (int j = i; j < std::min(size, i + int_time_ind); ++j) {

//             int lag = j - i;

//             int i4 = 4*i;
//             int j4 = 4*j;

//             corr_x[lag] += dat[i4] * dat[j4];
//             corr_y[lag] += dat[i4 + 1] * dat[j4 + 1];
//             corr_z[lag] += dat[i4 + 2] * dat[j4 + 2];

//             ++ntime[lag];
//         }
//     }

//     for (int i = 0; i < int_time_ind; ++i) {
//         corr_x[i] /= ntime[i];
//         corr_y[i] /= ntime[i];
//         corr_z[i] /= ntime[i];
//     }

//     std::cerr << "* loaded " << int_time_ind
//               << " samples from '" << fname << "'"
//               << std::endl;

//     for (int k = 0; k < int_time_ind; ++k) {

//         double autocorr =
//             corr_x[k] + corr_y[k] + corr_z[k];

//         std::cout << k*dt << " " << autocorr << std::endl;
//     }
// }

// void do_autocorr(std::vector<double> dat, const double int_max, std::string fname, double dump)
// {
//   dt = dat[4] - dat[0];
//   int int_time_ind = (int)(int_max / dt);
//   int size = dat.size();

//   int nframes = dat.size() / 4;

//   std::vector<double> corr_vel(3 * int_time_ind, 0.0);
//   std::vector<int> nav(int_time_ind, 0);

//   for (int i = 0; i < nframes; ++i) {
//     for (int j = i; j < std::min(nframes, i + int_time_ind); ++j) {

//       int itime = j - i;
//       nav[itime]++;

//       for (int k = 0; k < 3; k++) {
//         double vi = dat[4*i + 1 + k];
//         double vj = dat[4*j + 1 + k];

//         corr_vel[3*itime + k] += vj * vi;
//       }
//     }
//   }

//   // double  corr_vel[3 * int_time_ind];
//   // int nav[int_time_ind];

//   // std::fill(nav, nav + int_time_ind, 0);
//   // std::fill(corr_vel, corr_vel + 3*int_time_ind, 0.0); 

//   // // autocorr
//   // for (int i = 0; i < size; ++i) {
//   //   for (int j = i; j < std::min(size, i + int_time_ind); ++j) {
//   //     int itime = j - i;

//   //     if (itime > size) continue;

//   //     nav[itime] += 1;

//   //     for (int k = 0; k < 3; k++) {
//   //       corr_vel[3 * i + k] += dat[4 * j + k + 1] * dat[4 * i + k + 1];
//   //     }
//   //   }
//   // }

//   // print tcf
//   for (int l = 0; l < std::min(size, int_time_ind); ++l) {
//     double time = dt * l;
//     double corr_vel_tot = 0;

//     for (int m = 0; m < 3; m++) {
//       corr_vel_tot += corr_vel[3 * l + m] / nav[l];
//     }

//     std::cout << time << " " << corr_vel_tot << std::endl;
//   }
  
  
// }

void do_autocorr(std::vector<double> dat, double int_max, std::string fname, double dump)
{
    const int stride = 4;
    const int nframes = dat.size() / stride;

    if (nframes < 2) {
        std::cerr << "Not enough frames\n";
        return;
    }

    double dt = dat[stride] - dat[0];

    int int_time_ind = static_cast<int>(int_max / dt);
    if (int_time_ind > nframes)
        int_time_ind = nframes;

    std::vector<double> corr_vel(3 * int_time_ind, 0.0);
    std::vector<int> nav(int_time_ind, 0);

    // --- PARALLEL REGION ---
#pragma omp parallel
    {
        std::vector<double> corr_private(3 * int_time_ind, 0.0);
        std::vector<int> nav_private(int_time_ind, 0);

#pragma omp for schedule(dynamic)
for(int i = 0; i < nframes; ++i){
	for(int j = i; j < std::min(nframes, i + int_time_ind); ++j){
                int itime = j - i;
                if (itime >= int_time_ind)
                    break;

                nav_private[itime]++;

                for (int k = 0; k < 3; ++k) {
                    double v1 = dat[stride * i + 1 + k];
                    double v2 = dat[stride * j + 1 + k];

                    corr_private[3 * itime + k] += v2 * v1;
                }
            }
        }

        // --- REDUCTION ---
#pragma omp critical
        {
            for (int i = 0; i < int_time_ind; ++i) {
                nav[i] += nav_private[i];
                for (int k = 0; k < 3; ++k) {
                    corr_vel[3 * i + k] += corr_private[3 * i + k];
                }
            }
        }
    }

    // --- OUTPUT ---
    for (int itime = 0; itime < int_time_ind; ++itime) {

        if (nav[itime] == 0) continue;

        double time = dt * itime;
        double corr_vel_tot = 0.0;

        for (int k = 0; k < 3; ++k) {
            corr_vel_tot += corr_vel[3 * itime + k] / nav[itime];
            // corr_vel_tot += corr_vel[3 * itime + k] / (nav[itime] * (dat.size() / (nframes * stride)));
            
        }

        std::cout << time << " " << corr_vel_tot << "\n";
    }
}

} // namespace

////////////////////////////////////////////////////////////////////////////////


int main(int argc, char **argv) {


  if (argc < 3) {
    std::cerr << argv[0] << " <file_name> <mb-spec.json>" << std::endl;
  }

  System sys;
  sys.SetUpFromJson("POWER", argv[argc - 1]);


  // Open files for Writing

  std::cout.setf(std::ios_base::showpoint);
  std::cout.precision(9);

  size_t lineno(0);
  size_t frameno(0);

  std::string fname = argv[1];
  std::ifstream ifs_cmd(fname);
  size_t pos = fname.find('.');
  std::string format  = fname.erase(0, pos + 1);
  std::transform(format.begin(), format.end(), format.begin(), ::tolower);

  size_t imcon(0);
  {
    imcon = sys.GetImcon();
  }
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
  dump = sys.GetDumpFrequency();
  std::vector<double> box(6);

  std::vector<double> everything;
  std::vector<opt::lammpstrj_frame> frames;


  while(!ifs_cmd.eof() ) {

    opt::lammpstrj_frame this_frame;
    if(format.compare("xyz") == 0) {
        opt::read_xyz_frame(lineno, imcon, ifs_cmd, this_frame, dump * frameno, box);
    } else {
        opt::read_lammpstrj_frame(lineno, imcon, ifs_cmd, this_frame, dump * frameno);
    }

    frameno++;

    if (ifs_cmd.eof()) break;
    if (this_frame.natm == 0) break;

    frames.push_back(this_frame);
  }

    int n_frames = frames.size();

    opt::lammpstrj_frame frame;

    double total_vel[3];
  

    for (int i = 0; i < n_frames; i++) {
        frame = frames[i];
        std::fill (total_vel, total_vel + 3, 0.0);

        for(size_t k = 0; k < frame.natm; ++k){
            size_t k3 = 3*k;
            size_t k4 = 4*k;
            size_t k9 = 9*k;

            for (size_t l = 0; l < 3; ++l) {
                total_vel[l] += frame.xyz[k3 + l]; // for all
                // total_vel[l] += frame.xyz[k3 + l] / 1e35; // for all
            }
        }
        everything.push_back(frame.time);
        for (size_t j = 0; j < 3; j++) {
          everything.push_back(total_vel[j]);
        }

        // everything.push_back(frame.time);
        // for (size_t k = 0; k < frame.natm; ++k) {
        //     everything.push_back(frame.xyz[3*k + 0]);
        //     everything.push_back(frame.xyz[3*k + 1]);
        //     everything.push_back(frame.xyz[3*k + 2]);
        // }
        

        
        // everything.push_back(frame.xyz[0]); 
        // everything.push_back(all_frame_vel[i3 + 1]); 
        // everything.push_back(all_frame_vel[i3 + 2]);
    }
    

    double int_max(0);
    {
    	int_max = sys.GetIntMax();
    	if (!int_max) {
    	    std::cerr << "could not convert '" << int_max << "' to real number" << std::endl;
    	    return EXIT_FAILURE;
    	}
    }

    // do_autocorr(all_frame_time, all_frame_vel, int_max, world_size, argv[1]);

    // if (world_rank == 0) {
    //     do_autocorr(everything, int_max, world_size, argv[1]);
    // }

    do_autocorr(everything, int_max, argv[1], dump);

    //std::cerr << "after autocorr" <<std::endl;

    if (ntraj == 0)
        exit(0);

    std::cout << std::setprecision(8) << std::scientific;

    for (int k = 0; k < nsamples; ++k) {
        std::cout << std::setw(18) << dt*k
                  << std::setw(18) << (accum[k]/ntrajs[k])
                  << '\n';
    }

    return 0;

}
