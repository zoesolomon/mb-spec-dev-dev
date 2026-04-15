#ifdef HAVE_CONFIG_H
#   include "config.h"
#endif // HAVE_CONFIG_H

#include <cmath>
#include <cassert>
#include <algorithm>

#include <iostream>
#include <iomanip>

#include "macros.h"

#include "ps.h"
#include "ttm4-periodic-es.h"
#include "ttm4-smear.h"

#include "mat_tools.h"

#include "constants.h"

////////////////////////////////////////////////////////////////////////////////

namespace {

//----------------------------------------------------------------------------//

const double CHARGECON = constants::CHARGECON;

//----------------------------------------------------------------------------//

// polarization parameters

const double polarO = 1.310;
const double polarH = 0.294;
const double polarM = 0.0;

const double polar_sqrt[4] = {
    std::sqrt(polarO),
    std::sqrt(polarH),
    std::sqrt(polarH),
    std::sqrt(polarM)
};

const double polar[4] = {
    polarO,
    polarH,
    polarH,
    polarM
};

const double polfacO = polarO;
const double polfacH = polarH;
const double polfacM = polarO;

const double AA[4] = {
    std::pow(polfacO, 1.0/6.0),
    std::pow(polfacH, 1.0/6.0),
    std::pow(polfacH, 1.0/6.0),
    std::pow(polfacM, 1.0/6.0)
};

const double aCC = 0.4;
const double aCD = 0.4;

const double aDD_inter    = 0.055;
const double aDD_intra_HH = 0.055;
const double aDD_intra_OH = 0.626;

//----------------------------------------------------------------------------//

// M-site positioning (TTM2.1-F)
const double gammaM = 0.426706882;

const double gamma1 = 1.0 - gammaM;
const double gamma2 = gammaM/2;

//----------------------------------------------------------------------------//

inline void compute_M_site_crd
    (const double* RESTRICT O,
     const double* RESTRICT H1,
     const double* RESTRICT H2,
     double* RESTRICT M)
{
    for (size_t i = 0; i < 3; ++i)
        M[i] = gamma1*O[i] + gamma2*(H1[i] + H2[i]);
}

//----------------------------------------------------------------------------//

inline void displacement(const double* cell,
      			 const double* r1, const double* r2, double* drij)
{
    double dx = r1[0] - r2[0];
    double dy = r1[1] - r2[1];
    double dz = r1[2] - r2[2];

    dx -= cell[0]*std::floor(dx/cell[0] + 0.5);
    dy -= cell[4]*std::floor(dy/cell[4] + 0.5);
    dz -= cell[8]*std::floor(dz/cell[8] + 0.5);

    drij[0] = dx;
    drij[1] = dy;
    drij[2] = dz;
}

//----------------------------------------------------------------------------//

inline double distance(const double* cell,
		       const double* r1, const double* r2)
{
    double drij[3];
    displacement(cell, r1, r2, drij);

    return std::sqrt(drij[0]*drij[0] + drij[1]*drij[1] + drij[2]*drij[2]);
}

//----------------------------------------------------------------------------//

inline void image(const double* cell, const double* r1, const double* r2,
		  double* crd_out)
{
    double dx = r1[0] - r2[0];
    double dy = r1[1] - r2[1];
    double dz = r1[2] - r2[2];

    double facx = std::floor(dx/cell[0] + 0.5);
    double facy = std::floor(dy/cell[4] + 0.5);
    double facz = std::floor(dz/cell[8] + 0.5);

    crd_out[0] = r2[0] + cell[0]*facx;
    crd_out[1] = r2[1] + cell[4]*facy;
    crd_out[2] = r2[2] + cell[8]*facz;
}

//----------------------------------------------------------------------------//

} // namespace

////////////////////////////////////////////////////////////////////////////////

namespace ttm {

//
// for CC-pol monomer geometry
//

const double ttm4_periodic_es::molecular_polarizability = 1.43016; // A^3
const double ttm4_periodic_es::molecular_dipole = 1.867895838174; // D

//----------------------------------------------------------------------------//

ttm4_periodic_es::ttm4_periodic_es()
: m_nw(0),
  m_mem(0)
{
}

//----------------------------------------------------------------------------//

ttm4_periodic_es::~ttm4_periodic_es()
{
    if (m_nw > 0)
        delete[] m_mem;
}

//----------------------------------------------------------------------------//

void ttm4_periodic_es::allocate(size_t nw)
{
    if (m_nw > nw)
        return;

    delete[] m_mem;

    // see operator() below
    m_mem = new double[95*nw + 144*nw*nw];

    m_nw = nw;
}

//----------------------------------------------------------------------------//

double ttm4_periodic_es::operator()
    (size_t nw, const double* RESTRICT crd, const double* RESTRICT cell)
{
    assert(crd || nw == 0);

    allocate(nw);

    // setup pointers

    const size_t natom = 4*nw; // O H H M
    const size_t natom3 = 3*natom;

    double* xyz    = m_mem;           // 12*nw
    double* charge = xyz + natom3;    //  4*nw
    double* dipole = charge + natom;  // 12*nw
    double* diag   = dipole + natom3; // 12*nw
    double* grdq   = diag + natom3;   // 12*nw
    double* grd4   = grdq + 27*nw;    // 27*nw
    double* Efq    = grd4 + natom3;   // 12*nw
    double* phi    = Efq + natom3;    //  4*nw
    double* ddt    = phi + natom;     // 144*nw*nw
    
    // zero out Efq/phi

    std::fill(Efq, Efq + natom3 + natom + natom3*natom3, 0.0);

    // compute M-sites and charges

    for (size_t n = 0; n < nw; ++n) {
        const size_t n4 = 4*n;
        const size_t n9 = 9*n;
        const size_t n12 = 12*n;

	// Apply PBC to xyz monomer
        std::copy(crd + n9, crd + n9 + 3, xyz + n12);
	for(size_t i = 1; i < 3; ++i)
	    image(cell, crd + n9, crd + n9 + 3*i, xyz + n12 + 3*i);
	
        compute_M_site_crd(xyz + n12, xyz + n12 + 3, xyz + n12 + 6,
                           xyz + n12 + 9);

	double q3[3];
	ps::dms_nasa(0.0, 0.0, 0.0, xyz + n12, q3, 0, false);

        // TTM2.1-F assignment
        const double tmp = gamma2/gamma1;

        charge[n4 + 0] = 0.0;                                     // O
        charge[n4 + 1] = CHARGECON*(q3[1] + tmp*(q3[1] + q3[2])); // H1
        charge[n4 + 2] = CHARGECON*(q3[2] + tmp*(q3[1] + q3[2])); // H2
        charge[n4 + 3] = CHARGECON*q3[0]/gamma1;                  // M

    }

    for (size_t i = 0; i < natom; ++i) {
        const size_t i3 = 3*i;

        for (size_t j = i + 1; j < natom; ++j) {
            const size_t j3 = 3*j;

            double Rij[3], Rsq(0);
	    displacement(cell, xyz + i3, xyz + j3, Rij);
            for (size_t k = 0; k < 3; ++k) {
		//PBC
		//Rij[k] = xyz[i3 + k] - xyz[j3 + k];
                Rsq += Rij[k]*Rij[k];
            }

            // charge-charge

            const bool ij_from_same_water = (i/4 == j/4);

            if (!ij_from_same_water) {
                double ts0, ts1;
                ttm4::smear01(std::sqrt(Rsq), AA[i%4]*AA[j%4], aCC, ts0, ts1);

                phi[i] += ts0*charge[j];
                phi[j] += ts0*charge[i];

                for (size_t k = 0; k < 3; ++k) {
                    Efq[i3 + k] += charge[j]*ts1*Rij[k];
                    Efq[j3 + k] -= charge[i]*ts1*Rij[k];
                }
            }

            // dipole-dipole tensor

            const double aDD = ij_from_same_water ?
                (i%4 == 0 ? aDD_intra_OH : aDD_intra_HH) : aDD_inter;

            double ts1, ts2;
            ttm4::smear2(std::sqrt(Rsq), AA[i%4]*AA[j%4], aDD, ts1, ts2);

            double dd3[3][3];

            dd3[0][0] = 3.0*ts2*Rij[0]*Rij[0] - ts1;
            dd3[1][1] = 3.0*ts2*Rij[1]*Rij[1] - ts1;
            dd3[2][2] = 3.0*ts2*Rij[2]*Rij[2] - ts1;
            dd3[0][1] = 3.0*ts2*Rij[0]*Rij[1];
            dd3[0][2] = 3.0*ts2*Rij[0]*Rij[2];
            dd3[1][2] = 3.0*ts2*Rij[1]*Rij[2];
            dd3[1][0] = dd3[0][1];
            dd3[2][0] = dd3[0][2];
            dd3[2][1] = dd3[1][2];

            const double aiaj = polar_sqrt[i%4]*polar_sqrt[j%4];
            for (size_t k = 0; k < 3; ++k)
                for (size_t l = 0; l < 3; ++l)
                    ddt[natom3*(i3 + k) + j3 + l] = -aiaj*dd3[k][l];
        }
    }

    // populate the diagonal of the dipole-dipole tensor

    for (size_t i = 0; i < natom3; ++i)
        ddt[natom3*i + i] = 1.0;

    // perform Cholesky decomposition ddt = L*L^T
    // storing L in the lower triangle of ddt and
    // its diagonal in diag

    for (size_t i = 0; i < natom3; ++i) {
        for (size_t j = i; j < natom3; ++j) {

            double sum = ddt[i*natom3 + j];
            for (size_t k = 0; k < i; ++k)
                sum -= ddt[i*natom3 + k]*ddt[j*natom3 + k];

            if (i == j) {
                assert(sum > 0.0);
                diag[i] = std::sqrt(sum);
            } else {
                ddt[j*natom3 + i] = sum/diag[i];
            }
        }
    }

    // solve L*y = sqrt(a)*Efq storing y in dipole[]

    for (size_t i = 0; i < natom3; ++i) {
        double sum = polar_sqrt[(i/3)%4]*Efq[i];
        for (size_t k = 0; k < i; ++k)
            sum -= ddt[i*natom3 + k]*dipole[k];

        dipole[i] = sum/diag[i];
    }

    // solve L^T*x = y

    for (size_t i = natom3; i > 0; --i) {
        const size_t i1 = i - 1;
        double sum = dipole[i1];
        for (size_t k = i; k < natom3; ++k)
            sum -= ddt[k*natom3 + i1]*dipole[k];

        dipole[i1] = sum/diag[i1];
    }

    for (size_t i = 0; i < natom3; ++i)
        dipole[i] *= polar_sqrt[(i/3)%4];

    // compute the energies

    double E_elec(0.0);
    for (size_t i = 0; i < natom; ++i)
        E_elec += phi[i]*charge[i];

    E_elec /= 2;

    double E_ind(0.0);
    for (size_t i = 0; i < natom3; ++i)
        E_ind -= dipole[i]*Efq[i];

    E_ind /= 2;

    return E_elec + E_ind;

}

//----------------------------------------------------------------------------//

double ttm4_periodic_es::lapack
    (size_t nw, const double* RESTRICT crd, const double* RESTRICT cell)
{
    assert(crd || nw == 0);

    allocate(nw);

    // setup pointers

    const size_t natom = 4*nw; // O H H M
    const size_t natom3 = 3*natom;

    double* xyz    = m_mem;           // 12*nw
    double* charge = xyz + natom3;    //  4*nw
    double* dipole = charge + natom;  // 12*nw
    double* diag   = dipole + natom3; // 12*nw
    double* grdq   = diag + natom3;   // 12*nw
    double* grd4   = grdq + 27*nw;    // 27*nw
    double* Efq    = grd4 + natom3;   // 12*nw
    double* phi    = Efq + natom3;    //  4*nw
    double* ddt    = phi + natom;     // 144*nw*nw

    // zero out Efq/phi

    std::fill(Efq, Efq + natom3 + natom + natom3*natom3, 0.0);

    // compute M-sites and charges

    for (size_t n = 0; n < nw; ++n) {
        const size_t n4 = 4*n;
        const size_t n9 = 9*n;
        const size_t n12 = 12*n;

	// Apply PBC to xyz monomer
        std::copy(crd + n9, crd + n9 + 3, xyz + n12);
	for(size_t i = 1; i < 3; ++i)
	    image(cell, crd + n9, crd + n9 + 3*i, xyz + n12 + 3*i);
	
        compute_M_site_crd(xyz + n12, xyz + n12 + 3, xyz + n12 + 6,
                           xyz + n12 + 9);

	double q3[3];
	ps::dms_nasa(0.0, 0.0, 0.0, xyz + n12, q3, 0, false);

        // TTM2.1-F assignment
        const double tmp = gamma2/gamma1;

        charge[n4 + 0] = 0.0;                                     // O
        charge[n4 + 1] = CHARGECON*(q3[1] + tmp*(q3[1] + q3[2])); // H1
        charge[n4 + 2] = CHARGECON*(q3[2] + tmp*(q3[1] + q3[2])); // H2
        charge[n4 + 3] = CHARGECON*q3[0]/gamma1;                  // M

    }

    std::fill(ddt, ddt + natom3*natom3, 0.0);
    for (size_t i = 0; i < natom; ++i) {
        const size_t i3 = 3*i;

        for (size_t j = i + 1; j < natom; ++j) {
            const size_t j3 = 3*j;

            double Rij[3], Rsq(0);
	    displacement(cell, xyz + i3, xyz + j3, Rij);
            for (size_t k = 0; k < 3; ++k) {
		//PBC
		//Rij[k] = xyz[i3 + k] - xyz[j3 + k];
                Rsq += Rij[k]*Rij[k];
            }

            // charge-charge

            const bool ij_from_same_water = (i/4 == j/4);

            if (!ij_from_same_water) {
                double ts0, ts1;
                ttm4::smear01(std::sqrt(Rsq), AA[i%4]*AA[j%4], aCC, ts0, ts1);

                phi[i] += ts0*charge[j];
                phi[j] += ts0*charge[i];

                for (size_t k = 0; k < 3; ++k) {
                    Efq[i3 + k] += charge[j]*ts1*Rij[k];
                    Efq[j3 + k] -= charge[i]*ts1*Rij[k];
                }
            }

            // dipole-dipole tensor

            const double aDD = ij_from_same_water ?
                (i%4 == 0 ? aDD_intra_OH : aDD_intra_HH) : aDD_inter;

            double ts1, ts2;
            ttm4::smear2(std::sqrt(Rsq), AA[i%4]*AA[j%4], aDD, ts1, ts2);

            double dd3[3][3];

            dd3[0][0] = 3.0*ts2*Rij[0]*Rij[0] - ts1;
            dd3[1][1] = 3.0*ts2*Rij[1]*Rij[1] - ts1;
            dd3[2][2] = 3.0*ts2*Rij[2]*Rij[2] - ts1;
            dd3[0][1] = 3.0*ts2*Rij[0]*Rij[1];
            dd3[0][2] = 3.0*ts2*Rij[0]*Rij[2];
            dd3[1][2] = 3.0*ts2*Rij[1]*Rij[2];
            dd3[1][0] = dd3[0][1];
            dd3[2][0] = dd3[0][2];
            dd3[2][1] = dd3[1][2];

            for (size_t k = 0; k < 3; ++k)
                for (size_t l = 0; l < 3; ++l)
                    ddt[natom3*(i3 + k) + j3 + l] = -dd3[k][l];
        }
    }

    // populate the diagonal of the dipole-dipole tensor

    for (size_t i = 0; i < natom3; ++i)
        ddt[natom3*i + i] = polar[i%4];

    // ddt will now be ddt^{-1}
//    tools::invert(natom3, ddt);

//    tools::mat_mult(natom3,ddt,natom3,Efq,1,dipole);

    // compute the energies

    double E_elec(0.0);
    for (size_t i = 0; i < natom; ++i)
        E_elec += phi[i]*charge[i];

    E_elec /= 2;

    double E_ind(0.0);
    for (size_t i = 0; i < natom3; ++i)
        E_ind -= dipole[i]*Efq[i];

    E_ind /= 2;

    return E_elec + E_ind;

}

//----------------------------------------------------------------------------//

void ttm4_periodic_es::dip(const size_t& nw, double* dip)
{
    // Must keep in line with operator() above

    const size_t natom_M = 4*nw;
    const size_t natom_M_3 = 3*natom_M;

    //xyz is already imaged at the monomer level
    double* xyz    = m_mem;           // 12*nw
    double* charge = xyz + natom_M_3;    //  4*nw
    double* dipole = charge + natom_M;  // 12*nw

    std::fill(dip, dip + 3, 0.0);

    for (size_t n = 0; n < nw; ++n) {

//	const size_t n3 = 3*n;
	const size_t n4 = 4*n;
	const size_t n12 = 3*n4;

	for (size_t i = 0; i < 4; ++i)
	    for (size_t k = 0; k < 3; ++k)
		dip[k] += charge[n4 + i]*xyz[n12 + 3*i + k]
		        + dipole[n12 + 3*i + k];
    }

    for (size_t k = 0; k < 3; ++k)
	dip[k] *= constants::DEBYE/constants::CHARGECON;


}

//----------------------------------------------------------------------------//

void ttm4_periodic_es::dipind(const size_t& nw, double* dipind)
{
    // Must keep in line with operator() above

    const size_t natom_M = 4*nw;
    const size_t natom_M_3 = 3*natom_M;

    double* dipole = m_mem + natom_M_3 + natom_M;  // 12*nw
    double* temp_dip = new double[nw * 3];

    std::fill(dipind, dipind + 3, 0.0);
    std::fill(temp_dip, temp_dip + (nw * 3), 0.0);

    for (size_t n = 0; n < nw; ++n) {

    //	const size_t n3 = 3*n;
        const size_t n4 = 4*n;
        const size_t n12 = 3*n4;

        for (size_t i = 0; i < 4; ++i) {
            for (size_t k = 0; k < 3; ++k) {
                dipind[k] += dipole[n12 + 3*i + k];
                temp_dip[n * 3 + k] += dipole[n12 + 3*i + k];
            }
            
        }
    }
    for (size_t k = 0; k < 3; ++k)
	dipind[k] *= constants::DEBYE/constants::CHARGECON;


    for (size_t num = 0; num < nw; num++) {
        //std::cerr << "O" << "\t";
        for (size_t coord = 0; coord < 3; coord++) {
            std::cerr << temp_dip[num * 3 + coord] * constants::DEBYE/constants::CHARGECON << "\t";
        }
        std::cerr << std::endl;
        //std::cerr << "H" << "\t" << "0.0\t0.0\t0.0" << std::endl;
        //std::cerr << "H" << "\t" << "0.0\t0.0\t0.0" << std::endl;
    }


}

//----------------------------------------------------------------------------//

void ttm4_periodic_es::ptensor(const size_t& nw, double* polind)
{
    // Must keep in line with operator() above

    const size_t natom_M = 4*nw;
    const size_t natom_M_3 = 3*natom_M;

///////////
    double* xyz    = m_mem;           // 12*nw
    double* charge = xyz + natom_M_3;    //  4*nw
    double* dipole = charge + natom_M;  // 12*nw
    double* diag   = dipole + natom_M_3; // 12*nw
    double* grdq   = diag + natom_M_3;   // 12*nw
    double* grd4   = grdq + 27*nw;    // 27*nw
    double* Efq    = grd4 + natom_M_3;   // 12*nw
    double* phi    = Efq + natom_M_3;    //  4*nw
    double* ddt    = phi + natom_M;     // 144*nw*nw
///////////

    for (size_t i = 0; i < natom_M_3; ++i) {
	ddt[i*natom_M_3 + i] = 1.0/diag[i];

	for (size_t j = i + 1; j < natom_M_3; ++j) {

	    double sum(0);
	    for (size_t k = i; k < j; ++k)
		sum -= ddt[j*natom_M_3 + k]*ddt[k*natom_M_3 + i];

	    ddt[j*natom_M_3 + i] = sum/diag[j];
	}
    }

    std::fill(polind, polind + 9, 0.0);

    for (size_t i = 0; i < natom_M; ++i) {
	for (size_t j = 0; j < natom_M; ++j) {
	    for (size_t a = 0; a < 3; ++a) {
		for (size_t b = 0; b < 3; ++b) {

		    double sum(0);

		    const size_t ia = 3*i + a;
		    const size_t jb = 3*j + b;

		    for (size_t k = std::max(ia, jb); k < natom_M_3; ++k)
			sum += ddt[k*natom_M_3 + ia]*ddt[k*natom_M_3 + jb];

		    polind[3*a + b] += sum*polar_sqrt[i%4]*polar_sqrt[j%4];
		}
	    }
	}
    }

}

#if 0
//----------------------------------------------------------------------------//

void ttm4_periodic_es::ptensor_lapack(const size_t& nw, double* polind)
{
    // Must keep in line with operator() above

    const size_t natom_M = 4*nw;
    const size_t natom_M_3 = 3*natom_M;

///////////
    double* xyz    = m_mem;           // 12*nw
    double* charge = xyz + natom_M_3;    //  4*nw
    double* dipole = charge + natom_M;  // 12*nw
    double* diag   = dipole + natom_M_3; // 12*nw
    double* grdq   = diag + natom_M_3;   // 12*nw
    double* grd4   = grdq + 27*nw;    // 27*nw
    double* Efq    = grd4 + natom_M_3;   // 12*nw
    double* phi    = Efq + natom_M_3;    //  4*nw
    double* ddt    = phi + natom_M;     // 144*nw*nw
///////////

    std::fill(polind, polind + 9, 0.0);

    for (size_t i = 0; i < natom_M; ++i) {
	for (size_t j = 0; j < natom_M; ++j) {
	    for (size_t a = 0; a < 3; ++a) {
		for (size_t b = 0; b < 3; ++b) {

		    double sum(0);

		    const size_t ia = 3*i + a;
		    const size_t jb = 3*j + b;

		    for (size_t k = std::max(ia, jb); k < natom_M_3; ++k)
			sum += ddt[k*natom_M_3 + ia]*ddt[k*natom_M_3 + jb];

		    polind[3*a + b] += sum*polar_sqrt[i%4]*polar_sqrt[j%4];

		}
	    }
	}
    }

}

//----------------------------------------------------------------------------//
#endif

} // namespace ttm

////////////////////////////////////////////////////////////////////////////////
