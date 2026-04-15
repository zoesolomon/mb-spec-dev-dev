#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <sstream>

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
    std::vector<double> v3x;
    std::vector<double> v3y;
    std::vector<double> v3z;

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
        double t, dx, dy, dz, dm;
        ifs >> t >> dx >> dy >> dz >> dm;

        if (v3x.size() == 0 && ntraj == 0)
            t0 = t;

        if (v3x.size() == 1 && ntraj == 0)
            dt = t - t0;

        if (ifs) {
            v3x.push_back(dx*debye);
            v3y.push_back(dy*debye);
            v3z.push_back(dz*debye);
        }
    }


    int this_file_nsamples = (int) v3x.size();
    
    int int_time_ind = (int) (int_max / dt);

    double corr_x[int_time_ind];
    double corr_y[int_time_ind];
    double corr_z[int_time_ind];
    double ntime[int_time_ind];

    std::fill(corr_x, corr_x + int_time_ind, 0.0);
    std::fill(corr_y, corr_y + int_time_ind, 0.0);
    std::fill(corr_z, corr_z + int_time_ind, 0.0);
    std::fill(ntime, ntime + int_time_ind, 0.0);

    for(int i = 0; i < this_file_nsamples; ++i){
	for(int j = i; j < std::min(this_file_nsamples, i + int_time_ind); ++j){
	    int time = (j-i);

	    corr_x[time] += v3x[i]*v3x[j];
	    corr_y[time] += v3y[i]*v3y[j];
	    corr_z[time] += v3z[i]*v3z[j];
	    ++ntime[time];
	}
    }

    for(int i = 0; i < int_time_ind; ++i){
	corr_x[i] /= ntime[i];
    	corr_y[i] /= ntime[i];
    	corr_z[i] /= ntime[i];
    }

    std::cerr << "* loaded " << int_time_ind
              << " samples from '" << filename << "'"
              << std::endl;

    if (nsamples < int_time_ind){

	for (int k = 0; k < nsamples; ++k) {
	    const double autocorr =
		(corr_x[k] + corr_y[k] + corr_z[k]);

	    accum[k] += autocorr;
	    ++ntrajs[k];
	}

	for (int k = nsamples; k < int_time_ind; ++k) {
	    const double autocorr =
		(corr_x[k] + corr_y[k] + corr_z[k]);

	    accum.push_back(autocorr);
	    ntrajs.push_back(1.0);
	}

	nsamples = int_time_ind;

    }else{
	for (int k = 0; k < int_time_ind; ++k) {
	    const double autocorr =
		(corr_x[k] + corr_y[k] + corr_z[k]);
	    
	    accum[k] += autocorr;
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
    sys.SetUpFromJson("INFRARED", argv[argc - 1]);

    double int_max(0);
    {
    	int_max = sys.GetIntMax();
    	if (!std::is_same<decltype(int_max), double>::value) {
    	    std::cerr << "could not convert '" << int_max << "' to real number" << std::endl;
    	    return EXIT_FAILURE;
    	}
    }


    for (int i = 1; i < argc - 1; ++i) {
        do_autocorr(argv[i], i - 1, int_max);
    }

    if (ntraj == 0)
        exit(0);

    std::cout << std::setprecision(8) << std::scientific;

    for (int k = 0; k < nsamples; ++k) {
        std::cout << std::setw(18) << dt*k
                  << std::setw(18) << (accum[k]/ntrajs[k])
                  << '\n';
    }
}
