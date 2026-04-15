#ifdef HAVE_CONFIG_H
#   include "config.h"
#endif // HAVE_CONFIG_H

#include <cmath>
#include <cassert>
#include <algorithm>
#include <iostream>

#include "macros.h"

#include "ps.h"
#include "ttm4-dip.h"
#include "ttm4-smear.h"

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

} // namespace

////////////////////////////////////////////////////////////////////////////////

namespace ttm {

//
// for CC-pol monomer geometry
//

const double ttm4_dip::molecular_polarizability = 1.43016; // A^3
const double ttm4_dip::molecular_dipole = 1.867895838174; // D

//----------------------------------------------------------------------------//

ttm4_dip::ttm4_dip()
: m_nw(0),
  m_mem(0)
{
}

//----------------------------------------------------------------------------//

ttm4_dip::~ttm4_dip()
{
    if (m_nw > 0)
        delete[] m_mem;
}

//----------------------------------------------------------------------------//

void ttm4_dip::allocate(size_t nw)
{
    if (m_nw > nw)
        return;

    delete[] m_mem;

    // see operator() below
    m_mem = new double[56*nw + 144*nw*nw];

    m_nw = nw;
}

//----------------------------------------------------------------------------//

void ttm4_dip::operator()
    (size_t nw, const double* RESTRICT crd, double* dip) 
{
    assert(crd || nw == 0);

    allocate(nw);

    // setup pointers

    const size_t natom = 4*nw; // O H H M
    const size_t natom3 = 3*natom;

    double* xyz    = m_mem;           // 12*nw     //   3*natom
    double* charge = xyz + natom3;    //  4*nw     //     natom
    double* dipole = charge + natom;  // 12*nw     //   3*natom
    double* diag   = dipole + natom3; // 12*nw     //   3*natom
    double* Efq    = diag + natom3;   // 12*nw     //   3*natom
    double* phi    = Efq + natom3;    //  4*nw     //     natom
    double* ddt    = phi + natom;     //144*nw*nw  // 3*3*natom

    // zero out Efq/phi/ddt

    std::fill(Efq, Efq + natom3 + natom + natom3*natom3, 0.0);
    std::fill(dip, dip + 3, 0.0);

    // compute M-sites and charges

    for (size_t n = 0; n < nw; ++n) {
        const size_t n3 = 3*n;
        const size_t n4 = 4*n;
        const size_t n9 = 3*n3;
        const size_t n12 = 12*n;

        std::copy(crd + n9, crd + n9 + 9, xyz + n12);
        compute_M_site_crd(crd + n9, crd + n9 + 3, crd + n9 + 6,
                           xyz + n12 + 9);

        double q3[3];
	ps::dms_nasa(0.0, 0.0, 0.0, crd + n9, q3, 0, false);

        // TTM2.1-F assignment
        const double tmp = gamma2/gamma1;

	const double qO  = q3[0];
	const double qH1 = q3[1];
	const double qH2 = q3[2];

        charge[n4 + 0] = 0.0;                               // O
        charge[n4 + 1] = CHARGECON*(qH1 + tmp*(qH1 + qH2)); // H1
        charge[n4 + 2] = CHARGECON*(qH2 + tmp*(qH1 + qH2)); // H2
        charge[n4 + 3] = CHARGECON*qO/gamma1;               // M

    }

    for (size_t i = 0; i < natom; ++i) {
        const size_t i3 = 3*i;

        for (size_t j = i + 1; j < natom; ++j) {
            const size_t j3 = 3*j;

            double Rij[3], Rsq(0);
            for (size_t k = 0; k < 3; ++k) {
                Rij[k] = xyz[i3 + k] - xyz[j3 + k];
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

    std::fill(dip, dip + 3, 0.0);

    // Add N-body dipole to total dipole
    for (size_t i = 0; i < natom; ++i)
	for (size_t j = 0; j < 3; ++j)
	    dip[j] += dipole[3*i + j];

    for (size_t i = 0; i < 3; ++i)
	dip[i] *= constants::DEBYE/CHARGECON;

}

//----------------------------------------------------------------------------//

} // namespace ttm

////////////////////////////////////////////////////////////////////////////////
