#include <cmath>
#include <vector>
#include <iomanip>
#include <iostream>

#include <cassert>
#include <cstdlib>
#include <sstream>
#include "system.h"

#define ANISO yes

#define RAMAN yes
#define DESYM yes
//#define BE_CORR yes

namespace {

bool normalize = false;
double tau = 0.0;
double alpha = 0.0;
double max_int = 0.0;
const size_t num_nu = 10000; // frequency grid
const double omega_max = 9000.0;
const double ps2s = 1.0e-12; // s -> ps
const double pi2 = 2.0*M_PI;

const double c = 2.99792458e-02; // x cm/ps (speed of light)

const double hplanck = 6.62606876e-34; // m^2 kg s^{-1}
const double hbar = hplanck/(2.0*M_PI);
const double kboltz = 1.3806503e-23; // m^2 kg s^{-2} K^{-1}
//const double epsilon_vacuum = 8.8541878e-12;
const double light_velocity = c/ps2s/100.0; // m/s

//const double T = 298.15; // K

} // namespace

int main(int argc, char** argv)
{
    System sys;
    sys.SetUpFromJson("RAMAN", argv[argc - 1]);
    
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

    const double factor = 1.0/(c); // nu(ps-1) -> omega (cm-1)

    const double dt = time[1] - time[0]; // in ps

    const double nu_max = omega_max/factor; // ps-1
    const double nu_min = 0.0;

    const size_t nt = time.size()/2;
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

    if(normalize){
	int ind_max_int = (int) (max_int/dt);
	double area((corr[0] + corr[ind_max_int - 1])/2);
	for (size_t t = 1; t < ind_max_int; ++t)
	    area += corr[t];

	area *= (time[ind_max_int - 1] - time[0])/ind_max_int;

	for (size_t t = 0; t < nt; ++t){
	    corr[t] /= area;
	}
    }

    const double dnu = (nu_max - nu_min)/(num_nu - 1);

    std::vector<double> freq;
    std::vector<double> intensity;

#ifdef ANISO
    for(size_t i = 0; i < nt; ++i)
	corr[i] /= 10.0;
#endif

    for (size_t k = 0; k < num_nu; ++k) {
    	const double nu = nu_min + k*dnu; // ps-1

	// assume autocorrelation has faded to 0 by upper integration limit
	double accum = 0.5*corr[0];
	for (size_t t = 1; t < (nt - 1); ++t)
	    accum += std::cos(pi2*nu*t*dt)*corr[t];

	// changing even integral from \int_{-\infty}^{\infty} 
	// to 2*\int_0^{\infty}
	accum = 2.0*accum*dt; //ps

	freq.push_back(factor*nu);
	intensity.push_back(accum);
	//intensity.push_back(accum*factor_quant/100.0);

    }

#ifdef RAMAN
    for (size_t k = 0; k < freq.size(); ++k)
	//intensity[k] *= 1.0/std::pow((20500.0-freq[k]),4);
	intensity[k] *= freq[k]*light_velocity/std::pow((20500.0-freq[k]),4);
#endif

#ifdef DESYM
    for (size_t k = 0; k < freq.size(); ++k)
	intensity[k] *= std::tanh( 0.5*(1.0/T/kboltz)*hbar
				   * freq[k]*100.0*light_velocity );
				   //* pi2*freq[k]*100.0*light_velocity );
#endif

#ifdef BE_CORR
    for (size_t k = 0; k < freq.size(); ++k)
	intensity[k] /= freq[k];

#endif

    double max(0);
    for (size_t k = 0; k < freq.size(); ++k)
	if((std::fabs(intensity[k]) > max) && freq[k] > 1000.0)
	    max = std::fabs(intensity[k]);

    for (size_t k = 0; k < freq.size(); ++k)
	if(normalize)
	    std::cout << std::setw(15) << freq[k] 
 		      << std::setw(15) << intensity[k]/max << std::endl;
	else
	    std::cout << std::setw(15) << freq[k] 
 		      << std::setw(15) << intensity[k] << std::endl;

}
