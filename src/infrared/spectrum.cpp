#include <cmath>
#include <vector>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <random>

#include <cassert>
#include <cstdlib>
#include <sstream>

#include "system.h"

namespace {

bool normalize = false;
size_t kernel_size = 15.0;


const size_t num_nu = 10000; // frequency grid
const double omega_max = 9000.0;
const double ps2s = 1.0e-12; // s -> ps
const double pi2 = 2.0*M_PI;

const double c = 2.99792458e-02; // x cm/ps (speed of light)

const double hplanck = 6.62606876e-34;
const double hbar = hplanck/(2.0*M_PI);
const double kboltz = 1.3806503e-23;
const double epsilon_vacuum = 8.8541878e-12;
const double light_velocity = 2.99792e8; // m/s
const double conv_freq =  2.99792e10;

//const double length = 19.7291394*1e-10; // m
//const double volume = length*length*length; // m^3

//const double T = 298.15; // K

} // namespace

int main(int argc, char** argv)
{

    System sys;
    sys.SetUpFromJson("INFRARED", argv[argc - 1]);
    
    if (argc < 3) {
        std::cerr << "usage: " << argv[0] << " <file_name> <mb-spec.json>" << std::endl;
        return EXIT_FAILURE;
    }

    std::ifstream infile(argv[1]);
    if (!infile) {
        std::cerr << "Error: cannot open file " << argv[1] << "\n";
        return EXIT_FAILURE;
    }


    double T = sys.GetTemperature();
    double volume = sys.GetBoxVolume();
    double tau = sys.GetTau();
    double alpha = sys.GetAlpha();
    double max_int = sys.GetIntMax();

    {
        if(T == 0){
            T = 298.15;
        }
        if (!T) {
            std::cerr << "could not convert '" << T
            << "' to real number" << std::endl;
            return EXIT_FAILURE;
        }
    
        assert(T > 0.0);
    }
    {
        if(volume == 0){
            normalize = true;
            volume = 7679.34930337;
        }

        if (!volume) {
            std::cerr << "could not convert '" << volume
                << "' to real number" << std::endl;
            return EXIT_FAILURE;
        }
    
        assert(volume > 0.0);
	    volume *= 1e-30; // A^3 -> m^3
    }
    {
        if (!std::is_same<decltype(tau), double>::value) {
            std::cerr << "could not convert '" << tau
            << "' to real number" << std::endl;
            return EXIT_FAILURE;
        }
    }
    {
        if (!std::is_same<decltype(alpha), double>::value) {
            std::cerr << "could not convert '" << alpha
            << "' to real number" << std::endl;
            return EXIT_FAILURE;
        }
    }
    {
        if (!std::is_same<decltype(max_int), double>::value) {
            std::cerr << "could not convert '" << max_int
            << "' to real number" << std::endl;
            return EXIT_FAILURE;
        }
    }

    std::vector<double> time;
    std::vector<double> corr;

    size_t Ndata = 0;

    double t(0);
    double c(0);

    while (infile >> t >> c) {
        if (max_int > 0 && t >= max_int) {
            time.push_back(t);
            corr.push_back(0.0);
        } else {
            ++Ndata;
            time.push_back(t);
            corr.push_back(c);
        }
    }

    infile.close();

    std::cerr << "Ndata : " << Ndata << std::endl;
    corr[Ndata] = 0.0;


    double arg = 2.0 * M_PI / (double(2 * Ndata - 1));

    double kT = kboltz * T;
    double beta = 1/kT;
    double norm =(3.0 * hbar * light_velocity * volume * epsilon_vacuum);
    norm = 2.0 * M_PI / norm;

    double dt = (time[1] - time[0]) * ps2s;
    double dnu = 1.0 / (double(2 * Ndata - 1)) / dt;
    double arg1 = beta * hplanck * dnu;

    if(tau != 0)
        for (size_t t = 0; t < Ndata; ++t)
            corr[t] *= std::exp(-(t*dt)/(tau*ps2s));
//    if(alpha != 0){
//        const double alpha2 = alpha*alpha;
//        for (size_t t = 0; t < Ndata; ++t){
//            const double time2 = (t*dt)*(t*dt)/(ps2s*ps2s);
//            corr[t] *= std::exp(-alpha2*time2);
//	}
//    }
//

    std::vector<double> freq;
    std::vector<double> abs;

    for (size_t inu = 0; inu < Ndata; ++inu){

	double nu = dnu * inu;
	double beta_hplanck_nu = arg1 * inu;
	double tmp = beta_hplanck_nu / 2.0;
	double factor_quantum = nu * tmp;

	if (inu == 0) factor_quantum = 1.0;

	double absorbance = 0.0;

	for (size_t itime = Ndata - 1; itime > 0; --itime){

	    double nu_time = (double(itime)) * (double(inu)) * arg;
	    absorbance += std::cos(nu_time) * corr[itime];

	}

	for (size_t itime = 0; itime < Ndata; ++itime){

	    double nu_time = (double(itime)) * (double(inu)) * arg;
	    absorbance += std::cos(nu_time) * corr[itime];

	}

	factor_quantum *= norm * dt;
	absorbance *= factor_quantum/100;
        freq.push_back(nu/conv_freq);
	abs.push_back(absorbance);
	//std::cout << (nu/conv_freq) << '\t'<< absorbance << std::endl;

    }
  

  if (alpha != 0) {
       // gaussian kernel
      //size_t kernel_size = 51;
      //double sigma = 20.0;
      kernel_size *= alpha;
      double sigma = kernel_size / 3;
   
      // create and populate array 0 to n-1
      double center = (kernel_size - 1) / 2.0;
      std::vector<double> kernel(kernel_size);
      std::iota(kernel.begin(), kernel.end(), 0);
   
       // Gaussian function
      std::transform(kernel.begin(), kernel.end(), kernel.begin(), [sigma, center](double k) { 
   	return (1 / (sigma * sqrt(2 * M_PI))) * exp(-0.5 * pow((k - center) / sigma, 2.0)); 
      });
      
      // normalize 
      double sum = std::accumulate(kernel.begin(), kernel.end(), 0.0);
      std::transform(kernel.begin(), kernel.end(), kernel.begin(), [sum](double k) { 
   	return k / sum; 
      });
   
      // convolve 1d
      size_t signal_size = abs.size();
      size_t pad_width = kernel_size / 2;
   
      // padded array
      std::vector<double> padded(signal_size + (2 * pad_width));
      std::vector<double> convolved(signal_size);
   
      for (size_t f = 0; f < pad_width; f++) {
   	padded[f] = abs[pad_width - f];
      } 
      for (size_t m = pad_width; m < (pad_width + signal_size); m++) {
   	padded[m] = abs[m - pad_width];
      } 
      for (size_t b = pad_width + signal_size; b < ((2 * pad_width) + signal_size); b++) {
           padded[b] = abs[signal_size - (b - signal_size)];
      } 
      for (size_t c = 0; c < convolved.size(); c++) {
   	std::vector<double> pdt(kernel_size);
   	for (size_t k = 0; k < kernel_size; k++) {
               pdt[k] = kernel[k] * padded[c + k];
   	}
   	convolved[c] = std::accumulate(pdt.begin(), pdt.end(), 0.0);
      }
   
      for  (size_t n = 0; n < abs.size(); ++n) {
   	std::cout << freq[n] << '\t' << convolved[n] << std::endl;
       }
   }

   else {

      for  (size_t n = 0; n < abs.size(); ++n) {
   	std::cout << freq[n] << '\t' << abs[n] << std::endl;
       }

   }
    return 0;
}    
