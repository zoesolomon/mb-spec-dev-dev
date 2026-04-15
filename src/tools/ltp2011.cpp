#include <cmath>
#include <cassert>
#include <iostream>
#include <cstddef>
#include <algorithm>

#include "ltp2011.h"
#include "orient.h"
#include "mat_tools.h"
#include "constants.h"

////////////////////////////////////////////////////////////////////////////////

namespace {

//----------------------------------------------------------------------------//

double damp(const double& r)
{

    const double r_min = 3.5;
    const double r_max = 5.5;

    const double t = ( ( (r-r_min)/(r_max-r_min) ) - 0.5 );

    if( t <= -0.5)
        return 1.0;
    else if(t >= 0.5)
        return 0.0;
    else
	return 1.0 / ( std::exp(4.0*std::tan(M_PI*t)) + 1.0 );

}

//----------------------------------------------------------------------------//

double muOH(const double& r)
{

// Lorenzo Lodi:
// The following is a fit in even-tempered Gaussian functions of 
// 59 AQCC[7,5]/XP aug-cc-pV5Z dipoles for the ground X state of OH
// between 1 and 4 bohrs. Dipoles should be accurate to about 2.e-3. 
// Relativistic & core-corrections are NOT included.
// NB: There is a NASTY cancellation error with this fit which eats up 
// about 5 digits of accuracy. double precision reals are compulsory. 
// RMS of the fit is 5e-5 a.u. The largest residuals are about 0.7e-3 a.u.
// Note: Kahan summation algorithm (compensated summation) does not help here. 
// The problem is ill-conditioned
//       as sum_i |x_i| / |sum_i x_i| = ~ 3e+7
// It might be solved by using as basis orthogonalised Gaussians, 
// but I didn't deem it necessary at this time.
//
// Fit has reasonable behaviour for r-> +oo and also not too bad for r->0.

     const double a = 0.934; 
     const size_t n_coeffs(8);

     double powers_a[n_coeffs];

     powers_a[0] = a;
     for(size_t i = 0; i < n_coeffs; ++i)
	 powers_a[i] = powers_a[i-1]*a;
    
     double coeffs[n_coeffs];
     coeffs[0] =   -51632.29191732509;
     coeffs[1] =   336439.76799518000;
     coeffs[2] =  -948384.09499239320;
     coeffs[3] =  1499393.47990293540;
     coeffs[4] = -1436572.60883963100;
     coeffs[5] =   834824.61962008930;
     coeffs[6] =  -272815.32710332540;
     coeffs[7] =    38746.11983311285;
    
     const double r2 = r*r;

     double muOH(0.0);
     for(size_t i = 0; i < n_coeffs; ++i)
	 muOH += coeffs[i]*std::exp( -powers_a[i]*r2 );
    
     return muOH*r2;

}

//----------------------------------------------------------------------------//

} // namespace

////////////////////////////////////////////////////////////////////////////////

namespace h2o {

void ltp2011::cartesian(const size_t nw, const double* xyz,
			double* dip)
{    
    assert(nw == 1);

    // Determine intramolecular geometry
    double aHOH = tools::angle(xyz + 3, xyz, xyz + 6);
    double r1 = tools::distance(xyz + 3, xyz);
    double r2 = tools::distance(xyz + 6, xyz);

    // Place oxygen of molecule at origin
    double molec[9];
    for(size_t i = 0; i < 3; ++i){
	molec[i] = 0.0;
	molec[3 + i] = xyz[3 + i] - xyz[i];
	molec[6 + i] = xyz[6 + i] - xyz[i];
    }

    // Define molecular coordinate system. 
    // i.e., molecule in XY plane; O atm at origin, x-axis bisecting HOH angle
    //       H1 in 1st quadrant, H2 in 4th quadrant. Explicitly, 
    //
    // First, make things simpler (?) by rotating molecule into 2D

    double rot[9];
    double rot_3d_to_2d [9];
    {
	double x[3]; // H1-O
	tools::norm_vector(molec + 3, molec, x);

	double oh2[3]; // H2-O
	tools::norm_vector(molec + 6, molec, oh2);

	double z[3]; // points out of molecular plane
	tools::cross(x, oh2, z);

	double y[3]; // perpendicular to oh1 x out-of-plane
	tools::cross(x, z, y);

	// pack molecular frame
	for(size_t i = 0; i < 3; ++i){
	    rot_3d_to_2d [3*i + 0] = x[i];
	    rot_3d_to_2d [3*i + 1] = y[i];
	    rot_3d_to_2d [3*i + 2] = z[i];
	}

	// After invert, rotation matrix stored in rot_3d_to_2d
	tools::invert(3, rot_3d_to_2d );

	double molec_transpose[9]; // need xyz in rows
	tools::transpose(3, 3, molec, molec_transpose);

	// rotated molecule now stored in molec with xyz in rows and atm in col.
	tools::mat_mult(3, rot_3d_to_2d , 3, molec_transpose, 3, molec);

    }

    // Now, calulate rotation required to put bisector on x-axis
    {

	const double aHOH2 = aHOH/2.0*constants::deg2rad;
	const double val_h2_y = molec[5];
	
	// Atom H1 will be on the x-axis based on definition in prev. sxn.
	// If the bisector points down (val_h2_y is negative), rotate +(aHOH/2)
	// If the bisector points up   (val_h2_y is positive), rotate -(aHOH/2)
	double cos_a = std::cos(aHOH2); 
	double sin_a = std::sin(copysign(aHOH2, -1.0*val_h2_y)); 

	double rot_aHOH2[9] = { cos_a,  -sin_a,  0.0,
	                        sin_a,   cos_a,  0.0,
			          0.0,     0.0,  1.0};

	// Calculate total rotation matrix
	tools::mat_mult(3, rot_aHOH2, 3, rot_3d_to_2d, 3, rot);

    }

#ifdef DEBUG
    // reset molecule
    for(size_t i = 0; i < 3; ++i){
	molec[i] = 0.0;
	molec[3 + i] = xyz[3 + i] - xyz[i];
	molec[6 + i] = xyz[6 + i] - xyz[i];
    }
    double molec_transpose[9]; // need xyz in rows
    tools::transpose(3, 3, molec, molec_transpose);

    tools::mat_mult(3, rot, 3, molec_transpose, 3, molec);
#endif

    // rot will now store rot-1
    tools::invert(3, rot);

    // Calculate the molecular frame dipoles
    double muX(0), muY(0);
    h2o::ltp2011::operator()(r1, r2, std::cos(aHOH*constants::deg2rad), 
			     muX, muY);

    // Finally, rotate dipole back to lab frame and return
    double dip_mol[3] = {muX, muY, 0.0};
    tools::mat_mult(3, rot, 3, dip_mol, 1, dip);

}

void ltp2011::operator()(const double r1_in_Ang, const double r2_in_Ang, 
			 const double cos_theta, 
			 double& muX, double& muY)
{

//==============================================================================
// Greg Medders, adapted from: Lorenzo Lodi
// LTP2011 dipole-moment surface (DMS) for the ground-state water molecule: 
//         
// REF: L. Lodi, J. Tennyson and O.L Polyansky, 
//         "A global, high accuracy ab initio dipole moment
//          surface for the electronic ground state of the water molecule", 
//          J. Chem. Phys. 135, 034113 (2011)
//          dx.doi.org/10.1063/1.3604934
//
// INPUT:  r1/Angstrom, r2/Angstrom, cos( theta )
// OUTPUT: mu_Y(Debye), mu_X(Debye)
//
// Internally, this function uses a.u.
//
// mu_Y is the "small" component, perpendicular to the bond-angle bisector  
//     (=0 a.u. at equilibrium)
// mu_X is the "large" component, parallel      to the bond-angle bisector  
//     (=~0.730 a.u. at equilibrium)
// See the ref. for a more precise description of the orientation of the axes.
//    M O R E     D E T A I L S
// All surfaces are based on IC-MRCI+Q[8,10] energy-derivative dipoles in the 
// aug-cc-pCV6Z basis set, and incorporates a damping function to ensure a 
// reasonable behaviour upon dissociation r->+inf. LTP2011 includes 
// relativistic corrections. The Y component of the dipole has an 
// irregular behaviour for linear H-O-H geometries when one bond is stretched 
// beyond about 3 bohrs (due to surface intersection).
//
// NOTE: Numerically this fit is rather ill-conditioned, about 9 digits of 
//       accuracy are lost in the summation because of floating-point 
//       truncation error. Ugly as it is, it shouldn't be a problem in practice.
//==============================================================================

    const double r1 = r1_in_Ang/constants::Bohr_A;
    const double r2 = r2_in_Ang/constants::Bohr_A;

    const size_t n_r_X(8);
//    const size_t n_theta_X(8); 
    const size_t ncoeffs_X(200);

#define index_X(i,j) (indX[ncoeffs_X*(i-1)+(j-1)])    
    int indX[3*ncoeffs_X];
    std::fill(indX, indX + 3*ncoeffs_X, 0);
    index_X(1,  1) = 0; index_X(2,  1) = 0; index_X(3,  1) = 1; 
    index_X(1,  2) = 0; index_X(2,  2) = 0; index_X(3,  2) = 2; 
    index_X(1,  3) = 0; index_X(2,  3) = 0; index_X(3,  3) = 3; 
    index_X(1,  4) = 0; index_X(2,  4) = 0; index_X(3,  4) = 4; 
    index_X(1,  5) = 0; index_X(2,  5) = 0; index_X(3,  5) = 5; 
    index_X(1,  6) = 0; index_X(2,  6) = 0; index_X(3,  6) = 6; 
    index_X(1,  7) = 0; index_X(2,  7) = 0; index_X(3,  7) = 7; 
    index_X(1,  8) = 0; index_X(2,  8) = 0; index_X(3,  8) = 8; 
    index_X(1,  9) = 0; index_X(2,  9) = 2; index_X(3,  9) = 1; 
    index_X(1, 10) = 0; index_X(2, 10) = 2; index_X(3, 10) = 2; 
    index_X(1, 11) = 0; index_X(2, 11) = 2; index_X(3, 11) = 3; 
    index_X(1, 12) = 0; index_X(2, 12) = 2; index_X(3, 12) = 4; 
    index_X(1, 13) = 0; index_X(2, 13) = 2; index_X(3, 13) = 5; 
    index_X(1, 14) = 0; index_X(2, 14) = 2; index_X(3, 14) = 6; 
    index_X(1, 15) = 0; index_X(2, 15) = 2; index_X(3, 15) = 7; 
    index_X(1, 16) = 0; index_X(2, 16) = 2; index_X(3, 16) = 8; 
    index_X(1, 17) = 0; index_X(2, 17) = 4; index_X(3, 17) = 1; 
    index_X(1, 18) = 0; index_X(2, 18) = 4; index_X(3, 18) = 2; 
    index_X(1, 19) = 0; index_X(2, 19) = 4; index_X(3, 19) = 3; 
    index_X(1, 20) = 0; index_X(2, 20) = 4; index_X(3, 20) = 4; 
    index_X(1, 21) = 0; index_X(2, 21) = 4; index_X(3, 21) = 5; 
    index_X(1, 22) = 0; index_X(2, 22) = 4; index_X(3, 22) = 6; 
    index_X(1, 23) = 0; index_X(2, 23) = 4; index_X(3, 23) = 7; 
    index_X(1, 24) = 0; index_X(2, 24) = 4; index_X(3, 24) = 8; 
    index_X(1, 25) = 0; index_X(2, 25) = 6; index_X(3, 25) = 1; 
    index_X(1, 26) = 0; index_X(2, 26) = 6; index_X(3, 26) = 2; 
    index_X(1, 27) = 0; index_X(2, 27) = 6; index_X(3, 27) = 3; 
    index_X(1, 28) = 0; index_X(2, 28) = 6; index_X(3, 28) = 4; 
    index_X(1, 29) = 0; index_X(2, 29) = 6; index_X(3, 29) = 5; 
    index_X(1, 30) = 0; index_X(2, 30) = 6; index_X(3, 30) = 6; 
    index_X(1, 31) = 0; index_X(2, 31) = 6; index_X(3, 31) = 7; 
    index_X(1, 32) = 0; index_X(2, 32) = 6; index_X(3, 32) = 8; 
    index_X(1, 33) = 0; index_X(2, 33) = 8; index_X(3, 33) = 1; 
    index_X(1, 34) = 0; index_X(2, 34) = 8; index_X(3, 34) = 2; 
    index_X(1, 35) = 0; index_X(2, 35) = 8; index_X(3, 35) = 3; 
    index_X(1, 36) = 0; index_X(2, 36) = 8; index_X(3, 36) = 4; 
    index_X(1, 37) = 0; index_X(2, 37) = 8; index_X(3, 37) = 5; 
    index_X(1, 38) = 0; index_X(2, 38) = 8; index_X(3, 38) = 6; 
    index_X(1, 39) = 0; index_X(2, 39) = 8; index_X(3, 39) = 7; 
    index_X(1, 40) = 0; index_X(2, 40) = 8; index_X(3, 40) = 8; 
    index_X(1, 41) = 1; index_X(2, 41) = 0; index_X(3, 41) = 1; 
    index_X(1, 42) = 1; index_X(2, 42) = 0; index_X(3, 42) = 2; 
    index_X(1, 43) = 1; index_X(2, 43) = 0; index_X(3, 43) = 3; 
    index_X(1, 44) = 1; index_X(2, 44) = 0; index_X(3, 44) = 4; 
    index_X(1, 45) = 1; index_X(2, 45) = 0; index_X(3, 45) = 5; 
    index_X(1, 46) = 1; index_X(2, 46) = 0; index_X(3, 46) = 6; 
    index_X(1, 47) = 1; index_X(2, 47) = 0; index_X(3, 47) = 7; 
    index_X(1, 48) = 1; index_X(2, 48) = 0; index_X(3, 48) = 8; 
    index_X(1, 49) = 1; index_X(2, 49) = 2; index_X(3, 49) = 1; 
    index_X(1, 50) = 1; index_X(2, 50) = 2; index_X(3, 50) = 2; 
    index_X(1, 51) = 1; index_X(2, 51) = 2; index_X(3, 51) = 3; 
    index_X(1, 52) = 1; index_X(2, 52) = 2; index_X(3, 52) = 4; 
    index_X(1, 53) = 1; index_X(2, 53) = 2; index_X(3, 53) = 5; 
    index_X(1, 54) = 1; index_X(2, 54) = 2; index_X(3, 54) = 6; 
    index_X(1, 55) = 1; index_X(2, 55) = 2; index_X(3, 55) = 7; 
    index_X(1, 56) = 1; index_X(2, 56) = 2; index_X(3, 56) = 8; 
    index_X(1, 57) = 1; index_X(2, 57) = 4; index_X(3, 57) = 1; 
    index_X(1, 58) = 1; index_X(2, 58) = 4; index_X(3, 58) = 2; 
    index_X(1, 59) = 1; index_X(2, 59) = 4; index_X(3, 59) = 3; 
    index_X(1, 60) = 1; index_X(2, 60) = 4; index_X(3, 60) = 4; 
    index_X(1, 61) = 1; index_X(2, 61) = 4; index_X(3, 61) = 5; 
    index_X(1, 62) = 1; index_X(2, 62) = 4; index_X(3, 62) = 6; 
    index_X(1, 63) = 1; index_X(2, 63) = 4; index_X(3, 63) = 7; 
    index_X(1, 64) = 1; index_X(2, 64) = 4; index_X(3, 64) = 8; 
    index_X(1, 65) = 1; index_X(2, 65) = 6; index_X(3, 65) = 1; 
    index_X(1, 66) = 1; index_X(2, 66) = 6; index_X(3, 66) = 2; 
    index_X(1, 67) = 1; index_X(2, 67) = 6; index_X(3, 67) = 3; 
    index_X(1, 68) = 1; index_X(2, 68) = 6; index_X(3, 68) = 4; 
    index_X(1, 69) = 1; index_X(2, 69) = 6; index_X(3, 69) = 5; 
    index_X(1, 70) = 1; index_X(2, 70) = 6; index_X(3, 70) = 6; 
    index_X(1, 71) = 1; index_X(2, 71) = 6; index_X(3, 71) = 7; 
    index_X(1, 72) = 1; index_X(2, 72) = 6; index_X(3, 72) = 8; 
    index_X(1, 73) = 2; index_X(2, 73) = 0; index_X(3, 73) = 1; 
    index_X(1, 74) = 2; index_X(2, 74) = 0; index_X(3, 74) = 2; 
    index_X(1, 75) = 2; index_X(2, 75) = 0; index_X(3, 75) = 3; 
    index_X(1, 76) = 2; index_X(2, 76) = 0; index_X(3, 76) = 4; 
    index_X(1, 77) = 2; index_X(2, 77) = 0; index_X(3, 77) = 5; 
    index_X(1, 78) = 2; index_X(2, 78) = 0; index_X(3, 78) = 6; 
    index_X(1, 79) = 2; index_X(2, 79) = 0; index_X(3, 79) = 7; 
    index_X(1, 80) = 2; index_X(2, 80) = 0; index_X(3, 80) = 8; 
    index_X(1, 81) = 2; index_X(2, 81) = 2; index_X(3, 81) = 1; 
    index_X(1, 82) = 2; index_X(2, 82) = 2; index_X(3, 82) = 2; 
    index_X(1, 83) = 2; index_X(2, 83) = 2; index_X(3, 83) = 3; 
    index_X(1, 84) = 2; index_X(2, 84) = 2; index_X(3, 84) = 4; 
    index_X(1, 85) = 2; index_X(2, 85) = 2; index_X(3, 85) = 5; 
    index_X(1, 86) = 2; index_X(2, 86) = 2; index_X(3, 86) = 6; 
    index_X(1, 87) = 2; index_X(2, 87) = 2; index_X(3, 87) = 7; 
    index_X(1, 88) = 2; index_X(2, 88) = 2; index_X(3, 88) = 8; 
    index_X(1, 89) = 2; index_X(2, 89) = 4; index_X(3, 89) = 1; 
    index_X(1, 90) = 2; index_X(2, 90) = 4; index_X(3, 90) = 2; 
    index_X(1, 91) = 2; index_X(2, 91) = 4; index_X(3, 91) = 3; 
    index_X(1, 92) = 2; index_X(2, 92) = 4; index_X(3, 92) = 4; 
    index_X(1, 93) = 2; index_X(2, 93) = 4; index_X(3, 93) = 5; 
    index_X(1, 94) = 2; index_X(2, 94) = 4; index_X(3, 94) = 6; 
    index_X(1, 95) = 2; index_X(2, 95) = 4; index_X(3, 95) = 7; 
    index_X(1, 96) = 2; index_X(2, 96) = 4; index_X(3, 96) = 8; 
    index_X(1, 97) = 2; index_X(2, 97) = 6; index_X(3, 97) = 1; 
    index_X(1, 98) = 2; index_X(2, 98) = 6; index_X(3, 98) = 2; 
    index_X(1, 99) = 2; index_X(2, 99) = 6; index_X(3, 99) = 3; 
    index_X(1,100) = 2; index_X(2,100) = 6; index_X(3,100) = 4; 
    index_X(1,101) = 2; index_X(2,101) = 6; index_X(3,101) = 5; 
    index_X(1,102) = 2; index_X(2,102) = 6; index_X(3,102) = 6; 
    index_X(1,103) = 2; index_X(2,103) = 6; index_X(3,103) = 7; 
    index_X(1,104) = 2; index_X(2,104) = 6; index_X(3,104) = 8; 
    index_X(1,105) = 3; index_X(2,105) = 0; index_X(3,105) = 1; 
    index_X(1,106) = 3; index_X(2,106) = 0; index_X(3,106) = 2; 
    index_X(1,107) = 3; index_X(2,107) = 0; index_X(3,107) = 3; 
    index_X(1,108) = 3; index_X(2,108) = 0; index_X(3,108) = 4; 
    index_X(1,109) = 3; index_X(2,109) = 0; index_X(3,109) = 5; 
    index_X(1,110) = 3; index_X(2,110) = 0; index_X(3,110) = 6; 
    index_X(1,111) = 3; index_X(2,111) = 0; index_X(3,111) = 7; 
    index_X(1,112) = 3; index_X(2,112) = 0; index_X(3,112) = 8; 
    index_X(1,113) = 3; index_X(2,113) = 2; index_X(3,113) = 1; 
    index_X(1,114) = 3; index_X(2,114) = 2; index_X(3,114) = 2; 
    index_X(1,115) = 3; index_X(2,115) = 2; index_X(3,115) = 3; 
    index_X(1,116) = 3; index_X(2,116) = 2; index_X(3,116) = 4; 
    index_X(1,117) = 3; index_X(2,117) = 2; index_X(3,117) = 5; 
    index_X(1,118) = 3; index_X(2,118) = 2; index_X(3,118) = 6; 
    index_X(1,119) = 3; index_X(2,119) = 2; index_X(3,119) = 7; 
    index_X(1,120) = 3; index_X(2,120) = 2; index_X(3,120) = 8; 
    index_X(1,121) = 3; index_X(2,121) = 4; index_X(3,121) = 1; 
    index_X(1,122) = 3; index_X(2,122) = 4; index_X(3,122) = 2; 
    index_X(1,123) = 3; index_X(2,123) = 4; index_X(3,123) = 3; 
    index_X(1,124) = 3; index_X(2,124) = 4; index_X(3,124) = 4; 
    index_X(1,125) = 3; index_X(2,125) = 4; index_X(3,125) = 5; 
    index_X(1,126) = 3; index_X(2,126) = 4; index_X(3,126) = 6; 
    index_X(1,127) = 3; index_X(2,127) = 4; index_X(3,127) = 7; 
    index_X(1,128) = 3; index_X(2,128) = 4; index_X(3,128) = 8; 
    index_X(1,129) = 4; index_X(2,129) = 0; index_X(3,129) = 1; 
    index_X(1,130) = 4; index_X(2,130) = 0; index_X(3,130) = 2; 
    index_X(1,131) = 4; index_X(2,131) = 0; index_X(3,131) = 3; 
    index_X(1,132) = 4; index_X(2,132) = 0; index_X(3,132) = 4; 
    index_X(1,133) = 4; index_X(2,133) = 0; index_X(3,133) = 5; 
    index_X(1,134) = 4; index_X(2,134) = 0; index_X(3,134) = 6; 
    index_X(1,135) = 4; index_X(2,135) = 0; index_X(3,135) = 7; 
    index_X(1,136) = 4; index_X(2,136) = 0; index_X(3,136) = 8; 
    index_X(1,137) = 4; index_X(2,137) = 2; index_X(3,137) = 1; 
    index_X(1,138) = 4; index_X(2,138) = 2; index_X(3,138) = 2; 
    index_X(1,139) = 4; index_X(2,139) = 2; index_X(3,139) = 3; 
    index_X(1,140) = 4; index_X(2,140) = 2; index_X(3,140) = 4; 
    index_X(1,141) = 4; index_X(2,141) = 2; index_X(3,141) = 5; 
    index_X(1,142) = 4; index_X(2,142) = 2; index_X(3,142) = 6; 
    index_X(1,143) = 4; index_X(2,143) = 2; index_X(3,143) = 7; 
    index_X(1,144) = 4; index_X(2,144) = 2; index_X(3,144) = 8; 
    index_X(1,145) = 4; index_X(2,145) = 4; index_X(3,145) = 1; 
    index_X(1,146) = 4; index_X(2,146) = 4; index_X(3,146) = 2; 
    index_X(1,147) = 4; index_X(2,147) = 4; index_X(3,147) = 3; 
    index_X(1,148) = 4; index_X(2,148) = 4; index_X(3,148) = 4; 
    index_X(1,149) = 4; index_X(2,149) = 4; index_X(3,149) = 5; 
    index_X(1,150) = 4; index_X(2,150) = 4; index_X(3,150) = 6; 
    index_X(1,151) = 4; index_X(2,151) = 4; index_X(3,151) = 7; 
    index_X(1,152) = 4; index_X(2,152) = 4; index_X(3,152) = 8; 
    index_X(1,153) = 5; index_X(2,153) = 0; index_X(3,153) = 1; 
    index_X(1,154) = 5; index_X(2,154) = 0; index_X(3,154) = 2; 
    index_X(1,155) = 5; index_X(2,155) = 0; index_X(3,155) = 3; 
    index_X(1,156) = 5; index_X(2,156) = 0; index_X(3,156) = 4; 
    index_X(1,157) = 5; index_X(2,157) = 0; index_X(3,157) = 5; 
    index_X(1,158) = 5; index_X(2,158) = 0; index_X(3,158) = 6; 
    index_X(1,159) = 5; index_X(2,159) = 0; index_X(3,159) = 7; 
    index_X(1,160) = 5; index_X(2,160) = 0; index_X(3,160) = 8; 
    index_X(1,161) = 5; index_X(2,161) = 2; index_X(3,161) = 1; 
    index_X(1,162) = 5; index_X(2,162) = 2; index_X(3,162) = 2; 
    index_X(1,163) = 5; index_X(2,163) = 2; index_X(3,163) = 3; 
    index_X(1,164) = 5; index_X(2,164) = 2; index_X(3,164) = 4; 
    index_X(1,165) = 5; index_X(2,165) = 2; index_X(3,165) = 5; 
    index_X(1,166) = 5; index_X(2,166) = 2; index_X(3,166) = 6; 
    index_X(1,167) = 5; index_X(2,167) = 2; index_X(3,167) = 7; 
    index_X(1,168) = 5; index_X(2,168) = 2; index_X(3,168) = 8; 
    index_X(1,169) = 6; index_X(2,169) = 0; index_X(3,169) = 1; 
    index_X(1,170) = 6; index_X(2,170) = 0; index_X(3,170) = 2; 
    index_X(1,171) = 6; index_X(2,171) = 0; index_X(3,171) = 3; 
    index_X(1,172) = 6; index_X(2,172) = 0; index_X(3,172) = 4; 
    index_X(1,173) = 6; index_X(2,173) = 0; index_X(3,173) = 5; 
    index_X(1,174) = 6; index_X(2,174) = 0; index_X(3,174) = 6; 
    index_X(1,175) = 6; index_X(2,175) = 0; index_X(3,175) = 7; 
    index_X(1,176) = 6; index_X(2,176) = 0; index_X(3,176) = 8; 
    index_X(1,177) = 6; index_X(2,177) = 2; index_X(3,177) = 1; 
    index_X(1,178) = 6; index_X(2,178) = 2; index_X(3,178) = 2; 
    index_X(1,179) = 6; index_X(2,179) = 2; index_X(3,179) = 3; 
    index_X(1,180) = 6; index_X(2,180) = 2; index_X(3,180) = 4; 
    index_X(1,181) = 6; index_X(2,181) = 2; index_X(3,181) = 5; 
    index_X(1,182) = 6; index_X(2,182) = 2; index_X(3,182) = 6; 
    index_X(1,183) = 6; index_X(2,183) = 2; index_X(3,183) = 7; 
    index_X(1,184) = 6; index_X(2,184) = 2; index_X(3,184) = 8; 
    index_X(1,185) = 7; index_X(2,185) = 0; index_X(3,185) = 1; 
    index_X(1,186) = 7; index_X(2,186) = 0; index_X(3,186) = 2; 
    index_X(1,187) = 7; index_X(2,187) = 0; index_X(3,187) = 3; 
    index_X(1,188) = 7; index_X(2,188) = 0; index_X(3,188) = 4; 
    index_X(1,189) = 7; index_X(2,189) = 0; index_X(3,189) = 5; 
    index_X(1,190) = 7; index_X(2,190) = 0; index_X(3,190) = 6; 
    index_X(1,191) = 7; index_X(2,191) = 0; index_X(3,191) = 7; 
    index_X(1,192) = 7; index_X(2,192) = 0; index_X(3,192) = 8; 
    index_X(1,193) = 8; index_X(2,193) = 0; index_X(3,193) = 1; 
    index_X(1,194) = 8; index_X(2,194) = 0; index_X(3,194) = 2; 
    index_X(1,195) = 8; index_X(2,195) = 0; index_X(3,195) = 3; 
    index_X(1,196) = 8; index_X(2,196) = 0; index_X(3,196) = 4; 
    index_X(1,197) = 8; index_X(2,197) = 0; index_X(3,197) = 5; 
    index_X(1,198) = 8; index_X(2,198) = 0; index_X(3,198) = 6; 
    index_X(1,199) = 8; index_X(2,199) = 0; index_X(3,199) = 7; 
    index_X(1,200) = 8; index_X(2,200) = 0; index_X(3,200) = 8; 

    const size_t n_r_Y(9); 
//    const size_t n_theta_Y(8); 
    const size_t ncoeffs_Y (225);

    int indY[3*ncoeffs_Y];
    std::fill(indY, indY + 3*ncoeffs_Y, 0);
#define index_Y(i,j) (indY[ncoeffs_Y*(i-1)+(j-1)])    
    index_Y(1,  1) = 0; index_Y(2,  1) = 1; index_Y(3,  1) = 0;
    index_Y(1,  2) = 0; index_Y(2,  2) = 1; index_Y(3,  2) = 1;
    index_Y(1,  3) = 0; index_Y(2,  3) = 1; index_Y(3,  3) = 2;
    index_Y(1,  4) = 0; index_Y(2,  4) = 1; index_Y(3,  4) = 3;
    index_Y(1,  5) = 0; index_Y(2,  5) = 1; index_Y(3,  5) = 4;
    index_Y(1,  6) = 0; index_Y(2,  6) = 1; index_Y(3,  6) = 5;
    index_Y(1,  7) = 0; index_Y(2,  7) = 1; index_Y(3,  7) = 6;
    index_Y(1,  8) = 0; index_Y(2,  8) = 1; index_Y(3,  8) = 7;
    index_Y(1,  9) = 0; index_Y(2,  9) = 1; index_Y(3,  9) = 8;
    index_Y(1, 10) = 0; index_Y(2, 10) = 3; index_Y(3, 10) = 0;
    index_Y(1, 11) = 0; index_Y(2, 11) = 3; index_Y(3, 11) = 1;
    index_Y(1, 12) = 0; index_Y(2, 12) = 3; index_Y(3, 12) = 2;
    index_Y(1, 13) = 0; index_Y(2, 13) = 3; index_Y(3, 13) = 3;
    index_Y(1, 14) = 0; index_Y(2, 14) = 3; index_Y(3, 14) = 4;
    index_Y(1, 15) = 0; index_Y(2, 15) = 3; index_Y(3, 15) = 5;
    index_Y(1, 16) = 0; index_Y(2, 16) = 3; index_Y(3, 16) = 6;
    index_Y(1, 17) = 0; index_Y(2, 17) = 3; index_Y(3, 17) = 7;
    index_Y(1, 18) = 0; index_Y(2, 18) = 3; index_Y(3, 18) = 8;
    index_Y(1, 19) = 0; index_Y(2, 19) = 5; index_Y(3, 19) = 0;
    index_Y(1, 20) = 0; index_Y(2, 20) = 5; index_Y(3, 20) = 1;
    index_Y(1, 21) = 0; index_Y(2, 21) = 5; index_Y(3, 21) = 2;
    index_Y(1, 22) = 0; index_Y(2, 22) = 5; index_Y(3, 22) = 3;
    index_Y(1, 23) = 0; index_Y(2, 23) = 5; index_Y(3, 23) = 4;
    index_Y(1, 24) = 0; index_Y(2, 24) = 5; index_Y(3, 24) = 5;
    index_Y(1, 25) = 0; index_Y(2, 25) = 5; index_Y(3, 25) = 6;
    index_Y(1, 26) = 0; index_Y(2, 26) = 5; index_Y(3, 26) = 7;
    index_Y(1, 27) = 0; index_Y(2, 27) = 5; index_Y(3, 27) = 8;
    index_Y(1, 28) = 0; index_Y(2, 28) = 7; index_Y(3, 28) = 0;
    index_Y(1, 29) = 0; index_Y(2, 29) = 7; index_Y(3, 29) = 1;
    index_Y(1, 30) = 0; index_Y(2, 30) = 7; index_Y(3, 30) = 2;
    index_Y(1, 31) = 0; index_Y(2, 31) = 7; index_Y(3, 31) = 3;
    index_Y(1, 32) = 0; index_Y(2, 32) = 7; index_Y(3, 32) = 4;
    index_Y(1, 33) = 0; index_Y(2, 33) = 7; index_Y(3, 33) = 5;
    index_Y(1, 34) = 0; index_Y(2, 34) = 7; index_Y(3, 34) = 6;
    index_Y(1, 35) = 0; index_Y(2, 35) = 7; index_Y(3, 35) = 7;
    index_Y(1, 36) = 0; index_Y(2, 36) = 7; index_Y(3, 36) = 8;
    index_Y(1, 37) = 0; index_Y(2, 37) = 9; index_Y(3, 37) = 0;
    index_Y(1, 38) = 0; index_Y(2, 38) = 9; index_Y(3, 38) = 1;
    index_Y(1, 39) = 0; index_Y(2, 39) = 9; index_Y(3, 39) = 2;
    index_Y(1, 40) = 0; index_Y(2, 40) = 9; index_Y(3, 40) = 3;
    index_Y(1, 41) = 0; index_Y(2, 41) = 9; index_Y(3, 41) = 4;
    index_Y(1, 42) = 0; index_Y(2, 42) = 9; index_Y(3, 42) = 5;
    index_Y(1, 43) = 0; index_Y(2, 43) = 9; index_Y(3, 43) = 6;
    index_Y(1, 44) = 0; index_Y(2, 44) = 9; index_Y(3, 44) = 7;
    index_Y(1, 45) = 0; index_Y(2, 45) = 9; index_Y(3, 45) = 8;
    index_Y(1, 46) = 1; index_Y(2, 46) = 1; index_Y(3, 46) = 0;
    index_Y(1, 47) = 1; index_Y(2, 47) = 1; index_Y(3, 47) = 1;
    index_Y(1, 48) = 1; index_Y(2, 48) = 1; index_Y(3, 48) = 2;
    index_Y(1, 49) = 1; index_Y(2, 49) = 1; index_Y(3, 49) = 3;
    index_Y(1, 50) = 1; index_Y(2, 50) = 1; index_Y(3, 50) = 4;
    index_Y(1, 51) = 1; index_Y(2, 51) = 1; index_Y(3, 51) = 5;
    index_Y(1, 52) = 1; index_Y(2, 52) = 1; index_Y(3, 52) = 6;
    index_Y(1, 53) = 1; index_Y(2, 53) = 1; index_Y(3, 53) = 7;
    index_Y(1, 54) = 1; index_Y(2, 54) = 1; index_Y(3, 54) = 8;
    index_Y(1, 55) = 1; index_Y(2, 55) = 3; index_Y(3, 55) = 0;
    index_Y(1, 56) = 1; index_Y(2, 56) = 3; index_Y(3, 56) = 1;
    index_Y(1, 57) = 1; index_Y(2, 57) = 3; index_Y(3, 57) = 2;
    index_Y(1, 58) = 1; index_Y(2, 58) = 3; index_Y(3, 58) = 3;
    index_Y(1, 59) = 1; index_Y(2, 59) = 3; index_Y(3, 59) = 4;
    index_Y(1, 60) = 1; index_Y(2, 60) = 3; index_Y(3, 60) = 5;
    index_Y(1, 61) = 1; index_Y(2, 61) = 3; index_Y(3, 61) = 6;
    index_Y(1, 62) = 1; index_Y(2, 62) = 3; index_Y(3, 62) = 7;
    index_Y(1, 63) = 1; index_Y(2, 63) = 3; index_Y(3, 63) = 8;
    index_Y(1, 64) = 1; index_Y(2, 64) = 5; index_Y(3, 64) = 0;
    index_Y(1, 65) = 1; index_Y(2, 65) = 5; index_Y(3, 65) = 1;
    index_Y(1, 66) = 1; index_Y(2, 66) = 5; index_Y(3, 66) = 2;
    index_Y(1, 67) = 1; index_Y(2, 67) = 5; index_Y(3, 67) = 3;
    index_Y(1, 68) = 1; index_Y(2, 68) = 5; index_Y(3, 68) = 4;
    index_Y(1, 69) = 1; index_Y(2, 69) = 5; index_Y(3, 69) = 5;
    index_Y(1, 70) = 1; index_Y(2, 70) = 5; index_Y(3, 70) = 6;
    index_Y(1, 71) = 1; index_Y(2, 71) = 5; index_Y(3, 71) = 7;
    index_Y(1, 72) = 1; index_Y(2, 72) = 5; index_Y(3, 72) = 8;
    index_Y(1, 73) = 1; index_Y(2, 73) = 7; index_Y(3, 73) = 0;
    index_Y(1, 74) = 1; index_Y(2, 74) = 7; index_Y(3, 74) = 1;
    index_Y(1, 75) = 1; index_Y(2, 75) = 7; index_Y(3, 75) = 2;
    index_Y(1, 76) = 1; index_Y(2, 76) = 7; index_Y(3, 76) = 3;
    index_Y(1, 77) = 1; index_Y(2, 77) = 7; index_Y(3, 77) = 4;
    index_Y(1, 78) = 1; index_Y(2, 78) = 7; index_Y(3, 78) = 5;
    index_Y(1, 79) = 1; index_Y(2, 79) = 7; index_Y(3, 79) = 6;
    index_Y(1, 80) = 1; index_Y(2, 80) = 7; index_Y(3, 80) = 7;
    index_Y(1, 81) = 1; index_Y(2, 81) = 7; index_Y(3, 81) = 8;
    index_Y(1, 82) = 2; index_Y(2, 82) = 1; index_Y(3, 82) = 0;
    index_Y(1, 83) = 2; index_Y(2, 83) = 1; index_Y(3, 83) = 1;
    index_Y(1, 84) = 2; index_Y(2, 84) = 1; index_Y(3, 84) = 2;
    index_Y(1, 85) = 2; index_Y(2, 85) = 1; index_Y(3, 85) = 3;
    index_Y(1, 86) = 2; index_Y(2, 86) = 1; index_Y(3, 86) = 4;
    index_Y(1, 87) = 2; index_Y(2, 87) = 1; index_Y(3, 87) = 5;
    index_Y(1, 88) = 2; index_Y(2, 88) = 1; index_Y(3, 88) = 6;
    index_Y(1, 89) = 2; index_Y(2, 89) = 1; index_Y(3, 89) = 7;
    index_Y(1, 90) = 2; index_Y(2, 90) = 1; index_Y(3, 90) = 8;
    index_Y(1, 91) = 2; index_Y(2, 91) = 3; index_Y(3, 91) = 0;
    index_Y(1, 92) = 2; index_Y(2, 92) = 3; index_Y(3, 92) = 1;
    index_Y(1, 93) = 2; index_Y(2, 93) = 3; index_Y(3, 93) = 2;
    index_Y(1, 94) = 2; index_Y(2, 94) = 3; index_Y(3, 94) = 3;
    index_Y(1, 95) = 2; index_Y(2, 95) = 3; index_Y(3, 95) = 4;
    index_Y(1, 96) = 2; index_Y(2, 96) = 3; index_Y(3, 96) = 5;
    index_Y(1, 97) = 2; index_Y(2, 97) = 3; index_Y(3, 97) = 6;
    index_Y(1, 98) = 2; index_Y(2, 98) = 3; index_Y(3, 98) = 7;
    index_Y(1, 99) = 2; index_Y(2, 99) = 3; index_Y(3, 99) = 8;
    index_Y(1,100) = 2; index_Y(2,100) = 5; index_Y(3,100) = 0;
    index_Y(1,101) = 2; index_Y(2,101) = 5; index_Y(3,101) = 1;
    index_Y(1,102) = 2; index_Y(2,102) = 5; index_Y(3,102) = 2;
    index_Y(1,103) = 2; index_Y(2,103) = 5; index_Y(3,103) = 3;
    index_Y(1,104) = 2; index_Y(2,104) = 5; index_Y(3,104) = 4;
    index_Y(1,105) = 2; index_Y(2,105) = 5; index_Y(3,105) = 5;
    index_Y(1,106) = 2; index_Y(2,106) = 5; index_Y(3,106) = 6;
    index_Y(1,107) = 2; index_Y(2,107) = 5; index_Y(3,107) = 7;
    index_Y(1,108) = 2; index_Y(2,108) = 5; index_Y(3,108) = 8;
    index_Y(1,109) = 2; index_Y(2,109) = 7; index_Y(3,109) = 0;
    index_Y(1,110) = 2; index_Y(2,110) = 7; index_Y(3,110) = 1;
    index_Y(1,111) = 2; index_Y(2,111) = 7; index_Y(3,111) = 2;
    index_Y(1,112) = 2; index_Y(2,112) = 7; index_Y(3,112) = 3;
    index_Y(1,113) = 2; index_Y(2,113) = 7; index_Y(3,113) = 4;
    index_Y(1,114) = 2; index_Y(2,114) = 7; index_Y(3,114) = 5;
    index_Y(1,115) = 2; index_Y(2,115) = 7; index_Y(3,115) = 6;
    index_Y(1,116) = 2; index_Y(2,116) = 7; index_Y(3,116) = 7;
    index_Y(1,117) = 2; index_Y(2,117) = 7; index_Y(3,117) = 8;
    index_Y(1,118) = 3; index_Y(2,118) = 1; index_Y(3,118) = 0;
    index_Y(1,119) = 3; index_Y(2,119) = 1; index_Y(3,119) = 1;
    index_Y(1,120) = 3; index_Y(2,120) = 1; index_Y(3,120) = 2;
    index_Y(1,121) = 3; index_Y(2,121) = 1; index_Y(3,121) = 3;
    index_Y(1,122) = 3; index_Y(2,122) = 1; index_Y(3,122) = 4;
    index_Y(1,123) = 3; index_Y(2,123) = 1; index_Y(3,123) = 5;
    index_Y(1,124) = 3; index_Y(2,124) = 1; index_Y(3,124) = 6;
    index_Y(1,125) = 3; index_Y(2,125) = 1; index_Y(3,125) = 7;
    index_Y(1,126) = 3; index_Y(2,126) = 1; index_Y(3,126) = 8;
    index_Y(1,127) = 3; index_Y(2,127) = 3; index_Y(3,127) = 0;
    index_Y(1,128) = 3; index_Y(2,128) = 3; index_Y(3,128) = 1;
    index_Y(1,129) = 3; index_Y(2,129) = 3; index_Y(3,129) = 2;
    index_Y(1,130) = 3; index_Y(2,130) = 3; index_Y(3,130) = 3;
    index_Y(1,131) = 3; index_Y(2,131) = 3; index_Y(3,131) = 4;
    index_Y(1,132) = 3; index_Y(2,132) = 3; index_Y(3,132) = 5;
    index_Y(1,133) = 3; index_Y(2,133) = 3; index_Y(3,133) = 6;
    index_Y(1,134) = 3; index_Y(2,134) = 3; index_Y(3,134) = 7;
    index_Y(1,135) = 3; index_Y(2,135) = 3; index_Y(3,135) = 8;
    index_Y(1,136) = 3; index_Y(2,136) = 5; index_Y(3,136) = 0;
    index_Y(1,137) = 3; index_Y(2,137) = 5; index_Y(3,137) = 1;
    index_Y(1,138) = 3; index_Y(2,138) = 5; index_Y(3,138) = 2;
    index_Y(1,139) = 3; index_Y(2,139) = 5; index_Y(3,139) = 3;
    index_Y(1,140) = 3; index_Y(2,140) = 5; index_Y(3,140) = 4;
    index_Y(1,141) = 3; index_Y(2,141) = 5; index_Y(3,141) = 5;
    index_Y(1,142) = 3; index_Y(2,142) = 5; index_Y(3,142) = 6;
    index_Y(1,143) = 3; index_Y(2,143) = 5; index_Y(3,143) = 7;
    index_Y(1,144) = 3; index_Y(2,144) = 5; index_Y(3,144) = 8;
    index_Y(1,145) = 4; index_Y(2,145) = 1; index_Y(3,145) = 0;
    index_Y(1,146) = 4; index_Y(2,146) = 1; index_Y(3,146) = 1;
    index_Y(1,147) = 4; index_Y(2,147) = 1; index_Y(3,147) = 2;
    index_Y(1,148) = 4; index_Y(2,148) = 1; index_Y(3,148) = 3;
    index_Y(1,149) = 4; index_Y(2,149) = 1; index_Y(3,149) = 4;
    index_Y(1,150) = 4; index_Y(2,150) = 1; index_Y(3,150) = 5;
    index_Y(1,151) = 4; index_Y(2,151) = 1; index_Y(3,151) = 6;
    index_Y(1,152) = 4; index_Y(2,152) = 1; index_Y(3,152) = 7;
    index_Y(1,153) = 4; index_Y(2,153) = 1; index_Y(3,153) = 8;
    index_Y(1,154) = 4; index_Y(2,154) = 3; index_Y(3,154) = 0;
    index_Y(1,155) = 4; index_Y(2,155) = 3; index_Y(3,155) = 1;
    index_Y(1,156) = 4; index_Y(2,156) = 3; index_Y(3,156) = 2;
    index_Y(1,157) = 4; index_Y(2,157) = 3; index_Y(3,157) = 3;
    index_Y(1,158) = 4; index_Y(2,158) = 3; index_Y(3,158) = 4;
    index_Y(1,159) = 4; index_Y(2,159) = 3; index_Y(3,159) = 5;
    index_Y(1,160) = 4; index_Y(2,160) = 3; index_Y(3,160) = 6;
    index_Y(1,161) = 4; index_Y(2,161) = 3; index_Y(3,161) = 7;
    index_Y(1,162) = 4; index_Y(2,162) = 3; index_Y(3,162) = 8;
    index_Y(1,163) = 4; index_Y(2,163) = 5; index_Y(3,163) = 0;
    index_Y(1,164) = 4; index_Y(2,164) = 5; index_Y(3,164) = 1;
    index_Y(1,165) = 4; index_Y(2,165) = 5; index_Y(3,165) = 2;
    index_Y(1,166) = 4; index_Y(2,166) = 5; index_Y(3,166) = 3;
    index_Y(1,167) = 4; index_Y(2,167) = 5; index_Y(3,167) = 4;
    index_Y(1,168) = 4; index_Y(2,168) = 5; index_Y(3,168) = 5;
    index_Y(1,169) = 4; index_Y(2,169) = 5; index_Y(3,169) = 6;
    index_Y(1,170) = 4; index_Y(2,170) = 5; index_Y(3,170) = 7;
    index_Y(1,171) = 4; index_Y(2,171) = 5; index_Y(3,171) = 8;
    index_Y(1,172) = 5; index_Y(2,172) = 1; index_Y(3,172) = 0;
    index_Y(1,173) = 5; index_Y(2,173) = 1; index_Y(3,173) = 1;
    index_Y(1,174) = 5; index_Y(2,174) = 1; index_Y(3,174) = 2;
    index_Y(1,175) = 5; index_Y(2,175) = 1; index_Y(3,175) = 3;
    index_Y(1,176) = 5; index_Y(2,176) = 1; index_Y(3,176) = 4;
    index_Y(1,177) = 5; index_Y(2,177) = 1; index_Y(3,177) = 5;
    index_Y(1,178) = 5; index_Y(2,178) = 1; index_Y(3,178) = 6;
    index_Y(1,179) = 5; index_Y(2,179) = 1; index_Y(3,179) = 7;
    index_Y(1,180) = 5; index_Y(2,180) = 1; index_Y(3,180) = 8;
    index_Y(1,181) = 5; index_Y(2,181) = 3; index_Y(3,181) = 0;
    index_Y(1,182) = 5; index_Y(2,182) = 3; index_Y(3,182) = 1;
    index_Y(1,183) = 5; index_Y(2,183) = 3; index_Y(3,183) = 2;
    index_Y(1,184) = 5; index_Y(2,184) = 3; index_Y(3,184) = 3;
    index_Y(1,185) = 5; index_Y(2,185) = 3; index_Y(3,185) = 4;
    index_Y(1,186) = 5; index_Y(2,186) = 3; index_Y(3,186) = 5;
    index_Y(1,187) = 5; index_Y(2,187) = 3; index_Y(3,187) = 6;
    index_Y(1,188) = 5; index_Y(2,188) = 3; index_Y(3,188) = 7;
    index_Y(1,189) = 5; index_Y(2,189) = 3; index_Y(3,189) = 8;
    index_Y(1,190) = 6; index_Y(2,190) = 1; index_Y(3,190) = 0;
    index_Y(1,191) = 6; index_Y(2,191) = 1; index_Y(3,191) = 1;
    index_Y(1,192) = 6; index_Y(2,192) = 1; index_Y(3,192) = 2;
    index_Y(1,193) = 6; index_Y(2,193) = 1; index_Y(3,193) = 3;
    index_Y(1,194) = 6; index_Y(2,194) = 1; index_Y(3,194) = 4;
    index_Y(1,195) = 6; index_Y(2,195) = 1; index_Y(3,195) = 5;
    index_Y(1,196) = 6; index_Y(2,196) = 1; index_Y(3,196) = 6;
    index_Y(1,197) = 6; index_Y(2,197) = 1; index_Y(3,197) = 7;
    index_Y(1,198) = 6; index_Y(2,198) = 1; index_Y(3,198) = 8;
    index_Y(1,199) = 6; index_Y(2,199) = 3; index_Y(3,199) = 0;
    index_Y(1,200) = 6; index_Y(2,200) = 3; index_Y(3,200) = 1;
    index_Y(1,201) = 6; index_Y(2,201) = 3; index_Y(3,201) = 2;
    index_Y(1,202) = 6; index_Y(2,202) = 3; index_Y(3,202) = 3;
    index_Y(1,203) = 6; index_Y(2,203) = 3; index_Y(3,203) = 4;
    index_Y(1,204) = 6; index_Y(2,204) = 3; index_Y(3,204) = 5;
    index_Y(1,205) = 6; index_Y(2,205) = 3; index_Y(3,205) = 6;
    index_Y(1,206) = 6; index_Y(2,206) = 3; index_Y(3,206) = 7;
    index_Y(1,207) = 6; index_Y(2,207) = 3; index_Y(3,207) = 8;
    index_Y(1,208) = 7; index_Y(2,208) = 1; index_Y(3,208) = 0;
    index_Y(1,209) = 7; index_Y(2,209) = 1; index_Y(3,209) = 1;
    index_Y(1,210) = 7; index_Y(2,210) = 1; index_Y(3,210) = 2;
    index_Y(1,211) = 7; index_Y(2,211) = 1; index_Y(3,211) = 3;
    index_Y(1,212) = 7; index_Y(2,212) = 1; index_Y(3,212) = 4;
    index_Y(1,213) = 7; index_Y(2,213) = 1; index_Y(3,213) = 5;
    index_Y(1,214) = 7; index_Y(2,214) = 1; index_Y(3,214) = 6;
    index_Y(1,215) = 7; index_Y(2,215) = 1; index_Y(3,215) = 7;
    index_Y(1,216) = 7; index_Y(2,216) = 1; index_Y(3,216) = 8;
    index_Y(1,217) = 8; index_Y(2,217) = 1; index_Y(3,217) = 0;
    index_Y(1,218) = 8; index_Y(2,218) = 1; index_Y(3,218) = 1;
    index_Y(1,219) = 8; index_Y(2,219) = 1; index_Y(3,219) = 2;
    index_Y(1,220) = 8; index_Y(2,220) = 1; index_Y(3,220) = 3;
    index_Y(1,221) = 8; index_Y(2,221) = 1; index_Y(3,221) = 4;
    index_Y(1,222) = 8; index_Y(2,222) = 1; index_Y(3,222) = 5;
    index_Y(1,223) = 8; index_Y(2,223) = 1; index_Y(3,223) = 6;
    index_Y(1,224) = 8; index_Y(2,224) = 1; index_Y(3,224) = 7;
    index_Y(1,225) = 8; index_Y(2,225) = 1; index_Y(3,225) = 8;

    double dX[ncoeffs_X];
    double dY[ncoeffs_Y];
#define d_X(i) (dX[i-1])
#define d_Y(i) (dY[i-1])

    d_X(  1) = -4.3606625464226454E+01; d_Y(  1) = -6.0025296063971437E+03;
    d_X(  2) =  3.3027196259894445E+03; d_Y(  2) =  3.9291648647678376E+04;
    d_X(  3) = -1.2727332561105239E+04; d_Y(  3) = -1.0569378589871156E+05;
    d_X(  4) =  2.0478327360204032E+04; d_Y(  4) =  1.5020926617361890E+05;
    d_X(  5) = -1.7320391032929223E+04; d_Y(  5) = -1.2289367399421075E+05;
    d_X(  6) =  8.0857296986251222E+03; d_Y(  6) =  5.9210063162055972E+04;
    d_X(  7) = -1.9717556345513731E+03; d_Y(  7) = -1.6309896166413906E+04;
    d_X(  8) =  1.9607508538273396E+02; d_Y(  8) =  2.3058900193641894E+03;
    d_X(  9) =  3.0335075152822174E+02; d_Y(  9) = -1.2187808235920966E+02;
    d_X( 10) = -1.1709356490854843E+03; d_Y( 10) = -1.0438762294771004E+02;
    d_X( 11) =  1.9956778959599906E+03; d_Y( 11) = -2.6954003638798604E+03;
    d_X( 12) = -1.9712161529225559E+03; d_Y( 12) =  1.6112550047979807E+04;
    d_X( 13) =  1.2340220565548516E+03; d_Y( 13) = -3.5922088374701911E+04;
    d_X( 14) = -4.8970182591694174E+02; d_Y( 14) =  4.1234354044479434E+04;
    d_X( 15) =  1.1206177641027170E+02; d_Y( 15) = -2.6668830950832620E+04;
    d_X( 16) = -1.1105175444354245E+01; d_Y( 16) =  9.7672387842139869E+03;
    d_X( 17) =  1.6925709138547063E+01; d_Y( 17) = -1.8779090769474860E+03;
    d_X( 18) = -4.3912662965647542E+01; d_Y( 18) =  1.4585975086590042E+02;
    d_X( 19) =  2.4288671317884109E+01; d_Y( 19) = -2.1822372815887320E+01;
    d_X( 20) =  3.4365336338240013E+01; d_Y( 20) = -2.2110734126872558E+01;
    d_X( 21) = -5.5705027501924633E+01; d_Y( 21) =  3.7571119498786572E+02;
    d_X( 22) =  3.1458201995923446E+01; d_Y( 22) = -8.9280859359521492E+02;
    d_X( 23) = -8.0773946006338520E+00; d_Y( 23) =  9.9854772253553529E+02;
    d_X( 24) =  7.8402541108334844E-01; d_Y( 24) = -6.1509045913653972E+02;
    d_X( 25) =  1.3694677677176514E-01; d_Y( 25) =  2.1314747796376469E+02;
    d_X( 26) = -1.2184309628182746E+00; d_Y( 26) = -3.8767443468488636E+01;
    d_X( 27) =  3.9898597473516020E+00; d_Y( 27) =  2.8598573219642276E+00;
    d_X( 28) = -6.4466157734350418E+00; d_Y( 28) = -6.9625815262816104E-01;
    d_X( 29) =  5.6885223634108115E+00; d_Y( 29) =  6.6073020690482736E+00;
    d_X( 30) = -2.7918382966599893E+00; d_Y( 30) = -2.7882436566596880E+01;
    d_X( 31) =  7.1398484320161515E-01; d_Y( 31) =  6.1709270655274850E+01;
    d_X( 32) = -7.3916051669584704E-02; d_Y( 32) = -7.7645710806545594E+01;
    d_X( 33) =  2.6790328253412099E-04; d_Y( 33) =  5.7561024705354612E+01;
    d_X( 34) = -4.0414235845958046E-03; d_Y( 34) = -2.4882164185464717E+01;
    d_X( 35) =  1.5002012491777350E-02; d_Y( 35) =  5.7990635062305955E+00;
    d_X( 36) = -2.5395644888476454E-02; d_Y( 36) = -5.6292930951167364E-01;
    d_X( 37) =  2.3034273439407116E-02; d_Y( 37) =  1.0287016772296820E-03;
    d_X( 38) = -1.1575981466648955E-02; d_Y( 38) = -1.0860989693355805E-02;
    d_X( 39) =  3.0378956698768889E-03; d_Y( 39) =  4.7514016035023587E-02;
    d_X( 40) = -3.2474704221385764E-04; d_Y( 40) = -1.0624429525728374E-01;
    d_X( 41) =  1.7490126454435995E+02; d_Y( 41) =  1.3119592284215287E-01;
    d_X( 42) = -6.6062814319734498E+03; d_Y( 42) = -9.2928998177512767E-02;
    d_X( 43) =  2.4632070084248153E+04; d_Y( 43) =  3.7497752367926296E-02;
    d_X( 44) = -3.9170928839035245E+04; d_Y( 44) = -7.9895325261531980E-03;
    d_X( 45) =  3.2936121345909283E+04; d_Y( 45) =  6.9510952312157315E-04;
    d_X( 46) = -1.5323433860934820E+04; d_Y( 46) =  1.1107348110055093E+04;
    d_X( 47) =  3.7287418976416811E+03; d_Y( 47) = -7.3812023665533910E+04;
    d_X( 48) = -3.7025900128080684E+02; d_Y( 48) =  2.0230408994767399E+05;
    d_X( 49) = -3.8204879799564151E+02; d_Y( 49) = -2.9410739202800079E+05;
    d_X( 50) =  1.4816809813710424E+03; d_Y( 50) =  2.4753687850269134E+05;
    d_X( 51) = -2.5410134350081353E+03; d_Y( 51) = -1.2375368727427488E+05;
    d_X( 52) =  2.5251484275723342E+03; d_Y( 52) =  3.5882111890754779E+04;
    d_X( 53) = -1.5866169642140449E+03; d_Y( 53) = -5.4833204488015035E+03;
    d_X( 54) =  6.2936423488600121E+02; d_Y( 54) =  3.3225098766572773E+02;
    d_X( 55) = -1.4335551639070036E+02; d_Y( 55) =  1.1210052264027763E+02;
    d_X( 56) =  1.4095056124002440E+01; d_Y( 56) =  3.4399633045124938E+03;
    d_X( 57) = -1.3412212213906969E+01; d_Y( 57) = -2.0218419913111138E+04;
    d_X( 58) =  3.4861250143381767E+01; d_Y( 58) =  4.4955127263088827E+04;
    d_X( 59) = -1.9057659324804263E+01; d_Y( 59) = -5.1673000975081872E+04;
    d_X( 60) = -2.8186887285377452E+01; d_Y( 60) =  3.3546952728580436E+04;
    d_X( 61) =  4.5464201636816142E+01; d_Y( 61) = -1.2360422032248898E+04;
    d_X( 62) = -2.5773491567353631E+01; d_Y( 62) =  2.3971254975670017E+03;
    d_X( 63) =  6.6598208994000743E+00; d_Y( 63) = -1.8849571401080175E+02;
    d_X( 64) = -6.5198890475585358E-01; d_Y( 64) =  1.5801327026620129E+01;
    d_X( 65) = -3.8146332383576009E-02; d_Y( 65) =  3.0778086309463106E+01;
    d_X( 66) =  3.4491055845546725E-01; d_Y( 66) = -3.5038283671477348E+02;
    d_X( 67) = -1.1504040370477924E+00; d_Y( 67) =  8.1883162195943078E+02;
    d_X( 68) =  1.8871106806586795E+00; d_Y( 68) = -9.2646558530730545E+02;
    d_X( 69) = -1.6847004003047914E+00; d_Y( 69) =  5.8468471543774649E+02;
    d_X( 70) =  8.3441442563798773E-01; d_Y( 70) = -2.0984500327254500E+02;
    d_X( 71) = -2.1497937103777076E-01; d_Y( 71) =  4.0027073412467871E+01;
    d_X( 72) =  2.2393755010853056E-02; d_Y( 72) = -3.1483292751508998E+00;
    d_X( 73) = -2.1892849125475004E+02; d_Y( 73) =  2.8459849010778271E-01;
    d_X( 74) =  5.6962800355550971E+03; d_Y( 74) = -2.6711004988846980E+00;
    d_X( 75) = -2.0588561699539132E+04; d_Y( 75) =  1.1089900367103837E+01;
    d_X( 76) =  3.2364759842469120E+04; d_Y( 76) = -2.4258537179810446E+01;
    d_X( 77) = -2.7053799841186432E+04; d_Y( 77) =  3.0265838116033592E+01;
    d_X( 78) =  1.2543623220686441E+04; d_Y( 78) = -2.2283342391293445E+01;
    d_X( 79) = -3.0456876672294966E+03; d_Y( 79) =  9.5735821221169317E+00;
    d_X( 80) =  3.0198138552824821E+02; d_Y( 80) = -2.2183009214204503E+00;
    d_X( 81) =  1.9542895348778075E+02; d_Y( 81) =  2.1411987383908127E-01;
    d_X( 82) = -7.6245485172340705E+02; d_Y( 82) = -8.8954806908923638E+03;
    d_X( 83) =  1.3183009633221081E+03; d_Y( 83) =  6.0029031672360856E+04;
    d_X( 84) = -1.3213562538592669E+03; d_Y( 84) = -1.6752066997622303E+05;
    d_X( 85) =  8.3564896815556131E+02; d_Y( 85) =  2.4866939134875033E+05;
    d_X( 86) = -3.3223581216292450E+02; d_Y( 86) = -2.1454483436766942E+05;
    d_X( 87) =  7.5487293678030255E+01; d_Y( 87) =  1.1058141888443311E+05;
    d_X( 88) = -7.3736087346769636E+00; d_Y( 88) = -3.3344258355232014E+04;
    d_X( 89) =  3.9449696341125673E+00; d_Y( 89) =  5.3747705450634530E+03;
    d_X( 90) = -1.0362686658220809E+01; d_Y( 90) = -3.5245815839042189E+02;
    d_X( 91) =  5.8361834900792928E+00; d_Y( 91) = -4.7782641193291056E+01;
    d_X( 92) =  8.1898774691727567E+00; d_Y( 92) = -1.7770035056206107E+03;
    d_X( 93) = -1.3488763912672312E+01; d_Y( 93) =  1.0278888960119511E+04;
    d_X( 94) =  7.7195530796107050E+00; d_Y( 94) = -2.2787458652314963E+04;
    d_X( 95) = -2.0111792595416773E+00; d_Y( 95) =  2.6208096604709513E+04;
    d_X( 96) =  1.9856014939796296E-01; d_Y( 96) = -1.7057952242937361E+04;
    d_X( 97) =  2.6992892313586481E-03; d_Y( 97) =  6.3111728824057036E+03;
    d_X( 98) = -2.3169406981168095E-02; d_Y( 98) = -1.2312067912722996E+03;
    d_X( 99) =  7.6920653129892269E-02; d_Y( 99) =  9.7611219354381319E+01;
    d_X(100) = -1.2699177163437980E-01; d_Y(100) = -3.9142006385295645E+00;
    d_X(101) =  1.1419975744036037E-01; d_Y(101) = -1.5359912640847142E+01;
    d_X(102) = -5.6911577327582563E-02; d_Y(102) =  1.2776397874860595E+02;
    d_X(103) =  1.4731428730556217E-02; d_Y(103) = -2.9297032442760428E+02;
    d_X(104) = -1.5391529138923943E-03; d_Y(104) =  3.3596087026596683E+02;
    d_X(105) =  1.3622336107641650E+02; d_Y(105) = -2.1782781495290965E+02;
    d_X(106) = -2.7620215237450120E+03; d_Y(106) =  8.1111138066315789E+01;
    d_X(107) =  9.6961566179726942E+03; d_Y(107) = -1.6200769375291998E+01;
    d_X(108) = -1.5071634521499458E+04; d_Y(108) =  1.3474396226224599E+00;
    d_X(109) =  1.2525535485496084E+04; d_Y(109) = -2.8478794804414065E-02;
    d_X(110) = -5.7877475885727608E+03; d_Y(110) =  2.6593452089218772E-01;
    d_X(111) =  1.4022418853302088E+03; d_Y(111) = -1.0932653638288912E+00;
    d_X(112) = -1.3882010495019813E+02; d_Y(112) =  2.3736962936600321E+00;
    d_X(113) = -5.1784446024105819E+01; d_Y(113) = -2.9437369227216550E+00;
    d_X(114) =  2.0347186203791898E+02; d_Y(114) =  2.1552315808708045E+00;
    d_X(115) = -3.5517945792744649E+02; d_Y(115) = -9.2075081022630911E-01;
    d_X(116) =  3.5968370839488307E+02; d_Y(116) =  2.1211273777407769E-01;
    d_X(117) = -2.2941598322942446E+02; d_Y(117) = -2.0351663807559817E-02;
    d_X(118) =  9.1615227223264810E+01; d_Y(118) =  4.0244729236675139E+03;
    d_X(119) = -2.0803746205627249E+01; d_Y(119) = -2.7579194260415024E+04;
    d_X(120) =  2.0215951358823077E+00; d_Y(120) =  7.8286179604332443E+04;
    d_X(121) = -5.1223059493923984E-01; d_Y(121) = -1.1841790641238638E+05;
    d_X(122) =  1.3793274425197524E+00; d_Y(122) =  1.0438173383473718E+05;
    d_X(123) = -8.6111511230274118E-01; d_Y(123) = -5.5170356933581454E+04;
    d_X(124) = -9.3618549735765555E-01; d_Y(124) =  1.7149894029817311E+04;
    d_X(125) =  1.6723633643267704E+00; d_Y(125) = -2.8724756908097334E+03;
    d_X(126) = -9.7717597698039071E-01; d_Y(126) =  1.9822578233739478E+02;
    d_X(127) =  2.5754294537479439E-01; d_Y(127) =  1.0548770133245853E+01;
    d_X(128) = -2.5639662240507732E-02; d_Y(128) =  4.7160384631398483E+02;
    d_X(129) = -4.8664524863777615E+01; d_Y(129) = -2.6927523238212161E+03;
    d_X(130) =  8.2253513436822686E+02; d_Y(130) =  5.9529358189060440E+03;
    d_X(131) = -2.8104242401845022E+03; d_Y(131) = -6.8453123481617804E+03;
    d_X(132) =  4.3215344250088574E+03; d_Y(132) =  4.4601981106036874E+03;
    d_X(133) = -3.5712787276469517E+03; d_Y(133) = -1.6534222457115247E+03;
    d_X(134) =  1.6446990040199869E+03; d_Y(134) =  3.2343774421433045E+02;
    d_X(135) = -3.9761250447597183E+02; d_Y(135) = -2.5735116318981454E+01;
    d_X(136) =  3.9302009433029525E+01; d_Y(136) =  3.6760240764022001E-01;
    d_X(137) =  7.4642192204291860E+00; d_Y(137) =  3.1648925950707962E+00;
    d_X(138) = -2.9549045720626509E+01; d_Y(138) = -2.1074204175226441E+01;
    d_X(139) =  5.2096102541524488E+01; d_Y(139) =  4.7397684800605475E+01;
    d_X(140) = -5.3340409897036352E+01; d_Y(140) = -5.4950297427654732E+01;
    d_X(141) =  3.4353054530474992E+01; d_Y(141) =  3.6442449369556925E+01;
    d_X(142) = -1.3800578378807131E+01; d_Y(142) = -1.3973705057789175E+01;
    d_X(143) =  3.1368956403484844E+00; d_Y(143) =  2.8879156570819760E+00;
    d_X(144) = -3.0360358130678833E-01; d_Y(144) = -2.4945287573154928E-01;
    d_X(145) =  2.4705941986514546E-02; d_Y(145) = -1.1241905260007698E+03;
    d_X(146) = -6.9117806702358564E-02; d_Y(146) =  7.8203937269640610E+03;
    d_X(147) =  5.0433248065949954E-02; d_Y(147) = -2.2550301375257903E+04;
    d_X(148) =  3.2425384490856857E-02; d_Y(148) =  3.4682736691500708E+04;
    d_X(149) = -7.1011085503201343E-02; d_Y(149) = -3.1134366985278044E+04;
    d_X(150) =  4.3066331371839794E-02; d_Y(150) =  1.6797128474809051E+04;
    d_X(151) = -1.1522000939407917E-02; d_Y(151) = -5.3466702394040858E+03;
    d_X(152) =  1.1545216566446470E-03; d_Y(152) =  9.2113318815936827E+02;
    d_X(153) =  1.0455117395260515E+01; d_Y(153) = -6.5819873444136647E+01;
    d_X(154) = -1.5379526273858755E+02; d_Y(154) = -1.3629199079396130E+00;
    d_X(155) =  5.1262413047375753E+02; d_Y(155) = -6.6848897861893420E+01;
    d_X(156) = -7.8025562897946952E+02; d_Y(156) =  3.7931898119598372E+02;
    d_X(157) =  6.4134095430272794E+02; d_Y(157) = -8.3698743594581538E+02;
    d_X(158) = -2.9441601383245484E+02; d_Y(158) =  9.6148168185669510E+02;
    d_X(159) =  7.1026892747306960E+01; d_Y(159) = -6.2592613445075858E+02;
    d_X(160) = -7.0097630764862515E+00; d_Y(160) =  2.3178244914996230E+02;
    d_X(161) = -5.5181072854595925E-01; d_Y(161) = -4.5266862800091985E+01;
    d_X(162) =  2.1988689794610252E+00; d_Y(162) =  3.5923807281924383E+00;
    d_X(163) = -3.9118084258312251E+00; d_Y(163) = -8.8188292516697686E-03;
    d_X(164) =  4.0477572843179388E+00; d_Y(164) = -2.2668415176542567E-01;
    d_X(165) = -2.6330324078872316E+00; d_Y(165) =  1.2977989413177795E+00;
    d_X(166) =  1.0650938147472431E+00; d_Y(166) = -2.8690230287728440E+00;
    d_X(167) = -2.4261849843153982E-01; d_Y(167) =  3.3521400971742992E+00;
    d_X(168) =  2.3404996829526681E-02; d_Y(168) = -2.2606551635301457E+00;
    d_X(169) = -1.3330355566063190E+00; d_Y(169) =  8.8515537020313673E-01;
    d_X(170) =  1.7597899896022227E+01; d_Y(170) = -1.8718353979136282E-01;
    d_X(171) = -5.7368093491015827E+01; d_Y(171) =  1.6558293980750705E-02;
    d_X(172) =  8.6504254282326798E+01; d_Y(172) =  1.9837545634243168E+02;
    d_X(173) = -7.0750176940476308E+01; d_Y(173) = -1.3997062144072215E+03;
    d_X(174) =  3.2382040446805647E+01; d_Y(174) =  4.0931169349134307E+03;
    d_X(175) = -7.7965304444907169E+00; d_Y(175) = -6.3858066957452929E+03;
    d_X(176) =  7.6827845837534881E-01; d_Y(176) =  5.8197297979717314E+03;
    d_X(177) =  1.6229968998813238E-02; d_Y(177) = -3.1917314434027830E+03;
    d_X(178) = -6.4908083799454452E-02; d_Y(178) =  1.0346443125482040E+03;
    d_X(179) =  1.1618507169545289E-01; d_Y(179) = -1.8198254073487243E+02;
    d_X(180) = -1.2124387765841604E-01; d_Y(180) =  1.3322791574812072E+01;
    d_X(181) =  7.9581487575241194E-02; d_Y(181) =  1.1165936308050561E-01;
    d_X(182) = -3.2415275442599034E-02; d_Y(182) =  4.6716401640872505E+00;
    d_X(183) =  7.4020781937269575E-03; d_Y(183) = -2.6730580795888727E+01;
    d_X(184) = -7.1146892846662979E-04; d_Y(184) =  5.9012280279981155E+01;
    d_X(185) =  9.2649525643705655E-02; d_Y(185) = -6.7650478121726849E+01;
    d_X(186) = -1.1242874381538730E+00; d_Y(186) =  4.3866386197446985E+01;
    d_X(187) =  3.5950660241542902E+00; d_Y(187) = -1.6145447026892441E+01;
    d_X(188) = -5.3760392178762846E+00; d_Y(188) =  3.1251545236183915E+00;
    d_X(189) =  4.3774577026751516E+00; d_Y(189) = -2.4475078288219265E-01;
    d_X(190) = -1.9981707526528956E+00; d_Y(190) = -2.1572413034670078E+01;
    d_X(191) =  4.8021225340536289E-01; d_Y(191) =  1.5418085563717500E+02;
    d_X(192) = -4.7250502079100087E-02; d_Y(192) = -4.5630573530583376E+02;
    d_X(193) = -2.6924375660369038E-03; d_Y(193) =  7.2028687205061055E+02;
    d_X(194) =  3.0639958403565398E-02; d_Y(194) = -6.6435370200832926E+02;
    d_X(195) = -9.6417377890139655E-02; d_Y(195) =  3.6897713056122154E+02;
    d_X(196) =  1.4317113229849029E-01; d_Y(196) = -1.2123709316119232E+02;
    d_X(197) = -1.1613665480567770E-01; d_Y(197) =  2.1641503394243188E+01;
    d_X(198) =  5.2889764456063852E-02; d_Y(198) = -1.6107065097376747E+00;
    d_X(199) = -1.2689809454080206E-02; d_Y(199) = -4.8973249873083269E-03;
    d_X(200) =  1.2468169536632245E-03; d_Y(200) = -1.1852809969847300E-01;
                                        d_Y(201) =  7.0847529039930990E-01;
                                        d_Y(202) = -1.5740582223622930E+00;
                                        d_Y(203) =  1.7975475049534229E+00;
                                        d_Y(204) = -1.1537548422919350E+00;
                                        d_Y(205) =  4.1766798107112368E-01;
                                        d_Y(206) = -7.8846735007957136E-02;
                                        d_Y(207) =  5.9433544921461134E-03;
                                        d_Y(208) =  1.3200627540448489E+00;
                                        d_Y(209) = -9.5383998224375830E+00;
                                        d_Y(210) =  2.8500648410562967E+01;
                                        d_Y(211) = -4.5392119414633939E+01;
                                        d_Y(212) =  4.2237271808823294E+01;
                                        d_Y(213) = -2.3668495437070739E+01;
                                        d_Y(214) =  7.8486672857396815E+00;
                                        d_Y(215) = -1.4145069394761105E+00;
                                        d_Y(216) =  1.0634681558898884E-01;
                                        d_Y(217) = -3.4744476473319408E-02;
                                        d_Y(218) =  2.5316013291309014E-01;
                                        d_Y(219) = -7.6152525136575100E-01;
                                        d_Y(220) =  1.2199452241563369E+00;
                                        d_Y(221) = -1.1412953127550940E+00;
                                        d_Y(222) =  6.4282295895550057E-01;
                                        d_Y(223) = -2.1419872523397476E-01;
                                        d_Y(224) =  3.8777767747524089E-02;
                                        d_Y(225) = -2.9272691210419748E-03;

    //    ------------------------------------------------------------    //

    assert( std::fabs(cos_theta) < 1.0 );
    assert( r1 > 0.0 );
    assert( r2 > 0.0 );

//    const size_t max_param = std::max(ncoeffs_X, ncoeffs_Y);
    const size_t max_exp = std::max(n_r_X, n_r_Y) + 1;

    // Fitting Variables

    const double theta = std::acos(cos_theta);

    const double x1 = r1 + r2;
    const double x2 = r2 - r1;
    const double x3 = M_PI - theta;

    // Pre-compute powers of the variables

    double powers_x1[max_exp];
    double powers_x2[max_exp];
    double powers_x3[max_exp];

    powers_x1[0] = 1.0;
    powers_x2[0] = 1.0;
    powers_x3[0] = 1.0;

    for(size_t i = 1; i < max_exp; ++i){
	powers_x1[i] =  x1*powers_x1[i-1];
	powers_x2[i] =  x2*powers_x2[i-1];
    }

    for(size_t i = 1; i < max_exp; ++i)
	powers_x3[i] =  x3*powers_x3[i-1];

    muX = 0.0;
    for(size_t i = 0; i < ncoeffs_X; ++i)
	muX += dX[i] * powers_x1[index_X(1,i+1)] 
	             * powers_x2[index_X(2,i+1)] 
		     * powers_x3[index_X(3,i+1)];

    muY = 0.0;
    for(size_t i = 0; i < ncoeffs_Y; ++i)
	muY += dY[i] * powers_x1[index_Y(1,i+1)] 
	             * powers_x2[index_Y(2,i+1)] 
		     * powers_x3[index_Y(3,i+1)];

    // add damping/asymptotic terms
    const double damp1 = damp(r1);
    const double damp2 = damp(r2);
    const double value1 = muOH(r1)*( 1.0 - damp2 ) ;
    const double value2 = muOH(r2)*( 1.0 - damp1 ) ;
    
    const double cos_theta_half = std::sqrt( (1.0 + cos_theta)/2.0  );
    const double sin_theta_half = std::sqrt( (1.0 - cos_theta)/2.0  );
    
    muX = muX*damp1*damp2 + cos_theta_half*( value1 + value2  );
    muY = muY*damp1*damp2 + sin_theta_half*( value1 - value2  );

    muX *= constants::au2debye;
    muY *= constants::au2debye;
    
}

//----------------------------------------------------------------------------//

} // namespace h2o

////////////////////////////////////////////////////////////////////////////////
