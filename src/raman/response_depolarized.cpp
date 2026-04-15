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
    std::vector<double> bxx;
    std::vector<double> bxy;
    std::vector<double> bxz;
    std::vector<double> byy;
    std::vector<double> byz;
    std::vector<double> bzz;

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
        double t, my_axx, my_ayy, my_azz, my_axy, my_axz, my_ayz;
        ifs >> t >> my_axx >> my_axy >> my_axz >> my_ayy >> my_ayz >> my_azz;

	double a_iso = (my_axx + my_ayy + my_azz)/3.0;

        if (bxx.size() == 0 && ntraj == 0)
            t0 = t;

        if (bxx.size() == 1 && ntraj == 0)
            dt = t - t0;

        if (ifs) {
	    bxx.push_back((my_axx-a_iso));//*1e-30); // A^3 -> m^3
	    byy.push_back((my_ayy-a_iso));//*1e-30); // A^3 -> m^3
	    bzz.push_back((my_azz-a_iso));//*1e-30); // A^3 -> m^3
	    bxy.push_back((my_axy      ));//*1e-30); // A^3 -> m^3
	    bxz.push_back((my_axz      ));//*1e-30); // A^3 -> m^3
	    byz.push_back((my_ayz      ));//*1e-30); // A^3 -> m^3
        }
    }


    int this_file_nsamples = (int) bxx.size();
    int int_time_ind = (int) (int_max / dt);

    double corr_depolar[int_time_ind];
    double ntime[int_time_ind];

    std::fill(corr_depolar, corr_depolar + int_time_ind, 0.0);
    std::fill(ntime, ntime + int_time_ind, 0.0);

    for(int i = 0; i < this_file_nsamples; ++i){
	for(int j = i; j < std::min(this_file_nsamples, i + int_time_ind); ++j){
	    int time = (j-i);

	    // xx = xx*xx + xy*yx + xz*zx
	    double xx = bxx[i]*bxx[j] + bxy[i]*bxy[j] + bxz[i]*bxz[j];
		
	    // yy = yx*xy + yy*yy + yz*zy
	    double yy = bxy[i]*bxy[j] + byy[i]*byy[j] + byz[i]*byz[j];

	    // zz = zx*xz + zy*yz + zz*zz
	    double zz = bxz[i]*bxz[j] + byz[i]*byz[j] + bzz[i]*bzz[j];

	    corr_depolar[time] += xx + yy + zz;
	    ++ntime[time];
	}
    }

    for(size_t i = 0; i < int_time_ind; ++i)
	corr_depolar[i] /= ntime[i];

    std::cerr << "* loaded " << int_time_ind
              << " samples from '" << filename << "'"
              << std::endl;

    if (nsamples < int_time_ind){

	for (int k = 0; k < nsamples; ++k) {

	    accum[k] += corr_depolar[k];
	    ++ntrajs[k];
	}

	for (int k = nsamples; k < int_time_ind; ++k) {
	    accum.push_back(corr_depolar[k]);
	    ntrajs.push_back(1.0);
	}

	nsamples = int_time_ind;

    }else{
	for (int k = 0; k < int_time_ind; ++k) {
	    accum[k] += corr_depolar[k];
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
