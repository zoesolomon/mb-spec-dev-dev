#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

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

namespace {

namespace t = ttm;

typedef ttm::ttm4_dip ind_model;
static h2o::ltp2011 dip1b_model;

// ================= distance =================
inline double distance3(size_t imcon, const double* cell, const double* rcell,
                        const double* r1, const double* r2)
{
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

// ================= imaging =================
inline void image(const size_t imcon, const double* cell, const double* rcell,
                  const double* r1, const double* r2, double* out)
{
    double dx = r1[0] - r2[0];
    double dy = r1[1] - r2[1];
    double dz = r1[2] - r2[2];

    if (imcon == 1 || imcon == 2) {
        double fx = std::floor(dx/cell[0] + 0.5);
        double fy = std::floor(dy/cell[4] + 0.5);
        double fz = std::floor(dz/cell[8] + 0.5);

        out[0] = r2[0] + cell[0]*fx;
        out[1] = r2[1] + cell[4]*fy;
        out[2] = r2[2] + cell[8]*fz;
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

        out[0] = r1[0] - dx;
        out[1] = r1[1] - dy;
        out[2] = r1[2] - dz;
    }
}

inline void image_dimer(const size_t imcon, const double* cell,
                        const double* rcell, const double* r1, const double* r2,
                        double* dimer)
{
    std::fill(dimer, dimer + 18, 0.0);

    for (int i = 0; i < 3; ++i)
        dimer[i] = r1[i];

    for (int a = 1; a < 3; ++a)
        image(imcon, cell, rcell, r1, r1 + 3*a, dimer + 3*a);

    image(imcon, cell, rcell, r1, r2, dimer + 9);

    for (int a = 1; a < 3; ++a)
        image(imcon, cell, rcell, dimer + 9, r2 + 3*a, dimer + 9 + 3*a);
}

double f_switch2(const double& r, const double r2i, const double r2f)
{
    if (r > r2f) return 0.0;
    if (r > r2i) {
        double x = (r - r2i)/(r2f - r2i);
        return 1 + x*x*(2*x - 3);
    }
    return 1.0;
}

} // namespace

// ================= MAIN =================
int main(int argc, char** argv)
{
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <traj> <json>\n";
        return EXIT_FAILURE;
    }

    omp_set_nested(1);

    System sys;
    sys.SetUpFromJson("INFRARED", argv[2]);

    std::ifstream ifs(argv[1]);

    ind_model pot_ind;

    size_t imcon = sys.GetImcon();
    double r2i = sys.GetCutoff1();
    double r2f = sys.GetCutoff2();

    std::string mb = sys.GetMb();
    std::transform(mb.begin(), mb.end(), mb.begin(), ::tolower);

    size_t lineno = 0, frameno = 0;

    while (!ifs.eof()) {

        std::vector<opt::lammpstrj_frame> frames;

        while (frames.size() < 32 && !ifs.eof()) {
            opt::lammpstrj_frame f;
            opt::read_lammpstrj_frame(lineno, imcon, ifs, f, frameno++);
            if (f.natm == 0) break;
            frames.push_back(f);
        }

        int n_frames = frames.size();
        if (n_frames == 0) break;

        std::vector<std::array<double,3>> all_dip(n_frames);
        std::vector<double> all_mag(n_frames);
        std::vector<double> all_time(n_frames);

        #pragma omp parallel for schedule(dynamic)
        for (int f = 0; f < n_frames; ++f) {

            const auto& frame = frames[f];

            double total_dip[3] = {0,0,0};
            double dip2b_corr[3] = {0,0,0};
            double dip2b_ind[3]  = {0,0,0};
            double dip_ind[3]    = {0,0,0};

            size_t nmolec = frame.natm / 3;

            if (nmolec > 200) {

                #pragma omp parallel
                {
                    double local_corr[3] = {0,0,0};
                    double local_ind2b[3] = {0,0,0};
                    double local_ind[3] = {0,0,0};

                    double dimer[18];

                    #pragma omp for schedule(dynamic)
                    for (size_t i = 0; i < nmolec; ++i) {
                        for (size_t j = i+1; j < nmolec; ++j) {

                            double r = distance3(imcon, frame.cell, frame.rcell,
                                                 frame.xyz + 9*i, frame.xyz + 9*j);

                            if (r < r2f) {
                                double s = f_switch2(r, r2i, r2f);

                                image_dimer(imcon, frame.cell, frame.rcell,
                                            frame.xyz + 9*i, frame.xyz + 9*j, dimer);

                                double tmp[3];

                                h2o::dip_poly2b5::eval(dimer, dimer+9, tmp);
                                for (int k=0;k<3;k++) local_corr[k]+=s*tmp[k];

                                pot_ind(2, dimer, tmp);
                                for (int k=0;k<3;k++) {
                                    local_ind2b[k]+=(1-s)*tmp[k];
                                    local_ind[k]+=s*tmp[k];
                                }
                            }
                        }
                    }

                    #pragma omp critical
                    for (int k=0;k<3;k++) {
                        dip2b_corr[k]+=local_corr[k];
                        dip2b_ind[k]+=local_ind2b[k];
                        dip_ind[k]+=local_ind[k];
                    }
                }

            } else {
                double dimer[18];

                for (size_t i=0;i<nmolec;i++)
                for (size_t j=i+1;j<nmolec;j++) {

                    double r = distance3(imcon, frame.cell, frame.rcell,
                                         frame.xyz + 9*i, frame.xyz + 9*j);

                    if (r < r2f) {
                        double s = f_switch2(r,r2i,r2f);

                        image_dimer(imcon, frame.cell, frame.rcell,
                                    frame.xyz+9*i, frame.xyz+9*j, dimer);

                        double tmp[3];

                        h2o::dip_poly2b5::eval(dimer,dimer+9,tmp);
                        for (int k=0;k<3;k++) dip2b_corr[k]+=s*tmp[k];

                        pot_ind(2,dimer,tmp);
                        for (int k=0;k<3;k++) {
                            dip2b_ind[k]+=(1-s)*tmp[k];
                            dip_ind[k]+=s*tmp[k];
                        }
                    }
                }
            }

            for (int k=0;k<3;k++)
                total_dip[k]+=dip2b_corr[k]-dip2b_ind[k];

            double mag = std::sqrt(total_dip[0]*total_dip[0] +
                                   total_dip[1]*total_dip[1] +
                                   total_dip[2]*total_dip[2]);

            all_dip[f] = {total_dip[0], total_dip[1], total_dip[2]};
            all_mag[f] = mag;
            all_time[f] = frame.time;
        }

        for (int i=0;i<n_frames;i++)
            std::cout << all_time[i] << " "
                      << all_dip[i][0] << " "
                      << all_dip[i][1] << " "
                      << all_dip[i][2] << " "
                      << all_mag[i] << "\n";
    }

    return 0;
}