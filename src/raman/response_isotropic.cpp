#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <algorithm>

#include "system.h"

namespace {

static size_t ntraj = 0;
static int nsamples = 0;

static double t0;
static double dt;

static std::vector<double> accum;
static std::vector<double> ntrajs;

const double debye = 3.33564e-30; // Coulomb*m

void do_autocorr(const char* filename, const int traj_num, const double int_max)
{
    std::vector<double> aiso;

    std::ifstream ifs(filename);

    if (!ifs) {
        std::cerr << " ** Error ** : could not open '"
                  << filename << "' for reading"
                  << std::endl;
        exit(1);
    }

//    std::string line;
//    getline(ifs,line);
//    std::istringstream iss(line);
//    iss >> vol;

    while (ifs) {
        double t, axx, ayy, azz, tmp;
        ifs >> t >> axx >> tmp >> tmp >> ayy >> tmp >> azz;

        if (aiso.size() == 0 && ntraj == 0)
            t0 = t;

        if (aiso.size() == 1 && ntraj == 0)
            dt = t - t0;

        if (ifs) {
//	    aiso.push_back(1.0/3.0*(axx)*1e-30);
//	    aiso.push_back(1.0/3.0*(axx+ayy+azz)*1e-30); // A^3 -> m^3
	    aiso.push_back(1.0/3.0*(axx+ayy+azz)); // A^3 -> m^3
        }
    }


    int this_file_nsamples = (int) aiso.size();

    int int_time_ind = (int) (int_max / dt);

    double corr_iso[int_time_ind];
    double ntime[int_time_ind];

    std::fill(corr_iso, corr_iso + int_time_ind, 0.0);
    std::fill(ntime, ntime + int_time_ind, 0.0);

    for(int i = 0; i < this_file_nsamples; ++i){
	for(int j = i; j < std::min(this_file_nsamples, i + int_time_ind); ++j){
	    size_t time = (j-i);

	    corr_iso[time] += aiso[i]*aiso[j];
	    ++ntime[time];
	}
    }

    for(size_t i = 0; i < int_time_ind; ++i)
	corr_iso[i] /= ntime[i];

    std::cerr << "* loaded " << this_file_nsamples
              << " samples from '" << filename << "'"
              << std::endl;

    if (nsamples < this_file_nsamples){

	for (int k = 0; k < nsamples; ++k) {

	    accum[k] += corr_iso[k];
	    ++ntrajs[k];
	}

	for (int k = nsamples; k < int_time_ind; ++k) {
	    accum.push_back(corr_iso[k]);
	    ntrajs.push_back(1.0);
	}

	nsamples = int_time_ind;

    }else{
	for (int k = 0; k < int_time_ind; ++k) {
	    accum[k] += corr_iso[k];
	    ++ntrajs[k];
	}

    }

    ++ntraj;
}

} // namespace

int main(int argc, char **argv)
{
    if (argc < 3) {
        std::cerr << "usage: " << argv[0] << " <file_name(s)> <mb-spec.json>" << std::endl;
        return EXIT_FAILURE;
    }

    System sys;
    sys.SetUpFromJson("RAMAN", argv[argc - 1]);

    double int_max(0);
    {
    	int_max = sys.GetIntMax();
    	if (!std::is_same<decltype(int_max), double>::value) {
    	    std::cerr << "could not convert '" << int_max << "' to real number" << std::endl;
    	    return EXIT_FAILURE;
    	}
    }

    for (int i = 1; i < argc - 1; ++i)
        do_autocorr(argv[i], i - 1, int_max);

    if (ntraj == 0)
        exit(0);

    std::cout << std::setprecision(8) << std::scientific;

    for (int k = 0; k < nsamples; ++k) {
        std::cout << std::setw(18) << dt*k
                  << std::setw(18) << (accum[k]/ntrajs[k])
                  << '\n';
    }
}
