#include <cmath>
#include <vector>
#include <iomanip>
#include <iostream>

#include <cassert>
#include <cstdlib>
#include <sstream>
#include <stdexcept>

namespace {

//bool normalize = false;
double T = 298.15; // K
double tau = 0.0;
double alpha = 0.0;
double max_int = 0.0;
// int clump_rate = 0;

const size_t num_nu = 10000; // frequency grid
const double omega_max = 9000.0;
const double pi2 = 2.0*M_PI;

const double c = 2.99792458e-02; // x cm/ps (speed of light)
const double kb = 1.3806503e-23; // Boltzmann const

} // namespace

int main(int argc, char** argv)
{
//    double volume(0.0);

    if (argc != 5) {
	std::cerr << "usage: spectrum-sfg temp tau alpha max < corr_dip.dat > spectrum.dat" << std::endl;
	return EXIT_FAILURE;
    }

    {
        std::istringstream iss(argv[1]);
        iss >> T;
        if (!iss || !iss.eof()) {
    	std::cerr << "could not convert '" << argv[1]
    	    << "' to real number" << std::endl;
    	return EXIT_FAILURE;
        }
    }
#if 0
    {
        std::istringstream iss(argv[2]);
        iss >> volume;
	if(volume == 0){
	    normalize = true;
	    volume = 7679.34930337;
	}
        if (!iss || !iss.eof()) {
    	std::cerr << "could not convert '" << argv[2]
    	    << "' to real number" << std::endl;
    	return EXIT_FAILURE;
        }
    
        assert(volume > 0.0);
	volume *= 1e-30; // A^3 -> m^3
    }
#endif
    {
        std::istringstream iss(argv[2]);
        iss >> tau;
        if (!iss || !iss.eof()) {
    	std::cerr << "could not convert '" << argv[2]
    	    << "' to real number" << std::endl;
    	return EXIT_FAILURE;
        }
    }
    {
        std::istringstream iss(argv[3]);
        iss >> alpha;
        if (!iss || !iss.eof()) {
    	std::cerr << "could not convert '" << argv[3]
    	    << "' to real number" << std::endl;
    	return EXIT_FAILURE;
        }
    }
    {
        std::istringstream iss(argv[4]);
        iss >> max_int;
        if (!iss || !iss.eof()) {
    	std::cerr << "could not convert '" << argv[4]
    	    << "' to real number" << std::endl;
    	return EXIT_FAILURE;
        }
    }

    std::cerr << T << ' ' << tau << ' ' << alpha << ' ' << max_int << ' ' 
	      << std::endl;

    std::vector<double> time;
    std::vector<double> corr;

    double corr_max(0);
    double t_max(0);
    size_t lineno(0);
    while (std::cin) {
    	double t, co;

	std::string line;
	std::getline(std::cin, line);
	if(!std::cin.eof()){
	++lineno;

	std::istringstream iss(line);
	iss >> t >> co;
    	if (iss.fail()) {
	    std::ostringstream oss;
	    oss << "unexpected text at line " << lineno << " of the scan file";
	    throw std::runtime_error(oss.str());
	}

	if(max_int > 0 && t > max_int){
	    time.push_back(t);
	    corr.push_back(0.0);
	}else if (std::cin) {
	    time.push_back(t);
	    corr.push_back(co);
	    t_max = t;
	    corr_max = co;
	}
	}
    }

    std::cerr << corr_max << ' ' << t_max << std::endl;

    double corr_0 = corr[0];

    if(max_int > 0){
	for(size_t i = 0; i < time.size(); ++i){
	   corr[i] -= corr_0;
//	    corr[i] -= corr_max;
//          corr[i] /= corr[0];
    }
    }

    const double ps2wav = 1.0/(c); // nu(ps-1) -> omega (cm-1)

    const double dt = time[1] - time[0]; // in ps

    const double nu_max = omega_max/ps2wav; // ps-1
    const double nu_min = 0.0;

    const size_t nt = time.size()/2;

//    std::cerr << "nt = " << nt << std::endl;

    if(tau != 0)
	for (size_t t = 0; t < nt; ++t)
	    corr[t] *= std::exp(-(t*dt)/tau);
    if(alpha != 0){
	const double alpha2 = alpha*alpha;
	for (size_t t = 0; t < nt; ++t){
	    const double time2 = (t*dt)*(t*dt);
	    corr[t] *= std::exp(-alpha2*time2);
	}
    }

    const double dnu = (nu_max - nu_min)/(num_nu - 1);

    std::vector<double> freq;
    std::vector<double> intensity_im;
    std::vector<double> intensity_re;

    for (size_t k = 0; k < num_nu; ++k) {
	const double nu = nu_min + k*dnu; // ps-1
//	const double nu = nu_min + (k+1)*dnu; // ps-1

	// assume autocorrelation has faded to 0 by upper integration limit
	double accum_im = 0.5*corr[0];
	for (size_t t = 1; t < (nt - 1); ++t)
	    accum_im += std::cos(pi2*nu*t*dt)*corr[t];

	double accum_re = 0.0;
//	double accum_re = 0.5*corr[0];  // works
	for (size_t t = 1; t < (nt - 1); ++t){
	    accum_re += std::sin(pi2*nu*t*dt)*corr[t];
//	    accum_re += std::cos(pi2*nu*t*dt)*corr[t-1];  // works
	}

	accum_im *= dt; //ps
	accum_re *= -1.0*dt; //ps
//	accum_re *= dt; //ps //works

	freq.push_back(ps2wav*nu);
	intensity_im.push_back(accum_im);
	intensity_re.push_back(accum_re);

    }

    // prefactor of \omega_{IR}/{k_B T}
    // spectra reported in arbitrary units, so neglecting 1/{k_B T}
    for (size_t k = 0; k < freq.size(); ++k){
	intensity_im[k] *= freq[k];
	intensity_re[k] *= freq[k];
    }
/*
    // Normalization, just do relative to Im, but scale Re by the same factor
    double max(0);
    for (size_t k = 0; k < freq.size(); ++k)
	if((std::fabs(intensity_im[k]) > max) && freq[k] > 2500.0)
	    max = std::fabs(intensity_re[k]);

    if(normalize){
	std::cerr << "Normalization factor: " << max << std::endl;

	for (size_t k = 0; k < freq.size(); ++k){
	    intensity_im[k] /= max;
	    intensity_re[k] /= max;
	}
    }
*/

#if 0
    for (size_t k = 0; k < freq.size(); k += clump_rate){
	double avg_intensity(0);
	double avg_freq(0);

	if((k + clump_rate) > freq.size())
	    return 0;

	for(size_t i = 0; i < clump_rate; ++i){
	    avg_intensity += intensity[k + i];
	    avg_freq      += freq[k + i];
	}

	avg_intensity /= clump_rate;
	avg_freq      /= clump_rate;

    	std::cout << std::setw(15) << avg_freq 
 	          << std::setw(15) << avg_intensity << std::endl;
    }
#else
    for (size_t k = 0; k < freq.size(); ++k)
    	std::cout << std::setw(15) << freq[k]
 	          << std::setw(15) << intensity_im[k]
 	          << std::setw(15) << intensity_re[k]
		  << std::setw(15) << intensity_re[k]*intensity_re[k] + intensity_im[k]*intensity_im[k] << std::endl;

#endif

}
