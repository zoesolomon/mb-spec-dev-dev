#include <cmath>
#include <cstdlib>
#include <iostream>
#include "mat_tools.h"
#include <mkl_lapacke.h>
//#include "mkl_lapacke.h"

extern "C" {
    // LU decomoposition of a general matrix
//    void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
//
//    // generate inverse of a matrix given its LU decomposition
//    void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, 
//		 int* lwork, int* INFO);
    	
//    int dsyev(MKL_INT matrix_order, char jobz, char uplo, MKL_INT n,
//		   double* a, MKL_INT lda, double* w);
}

namespace tools {

double sgn(const double& x)
{
    return ( x > 0 ) - ( x < 0 );
}

void invert(int N, double* A)
{
    int *IPIV = new int[N+1];
    int LWORK = N*N;
    double *WORK = new double[LWORK];
    //int INFO;

    lapack_int n = N, info;
    info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n,n,A,n,IPIV);
    info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, A, n,IPIV);

    delete IPIV;
    delete WORK;
}

int eigenvalues(size_t N, double* H, double* eig)
{
    // H is N*N array
    
//    double* w = new double[N];
    double* a = new double[N*N];
    std::fill(a, a + N*N, 0.0);

    for(size_t i = 0; i < N; ++i){
	const size_t iN = i*N;
	for(size_t j = i; j < N; ++j){
	    a[iN+j] = H[iN+j];
	}
    }

    lapack_int n = N, lda = N, info;//, LAPACK_ROW_MAJOR = 101;
    info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'N', 'U', n, a, lda, eig);

    if(info != 0){
	std::cerr << "Failed to converge eigenvalues" << std::endl;
    }
    
    delete a;
    return info;
}

int eigenvalvect(size_t N, double* H, double* eigval, double* eigvect)
{
    // H is N*N array
    
    // Pack eigvect with Hamiltonian
    for(size_t i = 0; i < N; ++i){
	const size_t iN = i*N;
	for(size_t j = i; j < N; ++j){
	    eigvect[iN+j] = H[iN+j];
	}
    }

    lapack_int n = N, lda = N, info;//, LAPACK_ROW_MAJOR = 101;
    info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', n, eigvect, lda, eigval);

    // Upon return, eigvect now holds eigenvectors

    if(info != 0){
	std::cerr << "Failed to converge eigenvalues" << std::endl;
    }
    
    return info;
}

void transpose(const size_t& m, const size_t& n,
       	       const double* A, double* At)
{
    for(size_t i = 0; i < n; ++i){
	for(size_t j = 0; j < m; ++j){
	    At[i*m + j] = A[j*n + i];
	}
    }
}

double vector_2norm(const size_t& m, const double* x)
{
    double val(0);
    for(size_t i = 0; i < m; ++i){
	val += x[i]*x[i];
    }
    return std::sqrt(val);
}

double vector_2norm(const size_t& m, const size_t& n, 
			   const size_t& row, const double* A)
{
    double val(0);
    for(size_t i = 0; i < m; ++i){
	val += A[i*n + row]*A[i*n + row];
    }
    return std::sqrt(val);
}

void mat_mult(const size_t& m, const double* A,
	      const size_t& n, const double* B, 
	      const size_t& o, double* C)
{
    std::fill(C, C + m*o, 0.0);

    // A: mxn
    // B: nxo
    // C: mxo

    for(size_t i = 0; i < m; ++i)
	for(size_t j = 0; j < o; ++j)
	    for(size_t k = 0; k < n; ++k)
		C[i*o + j] += A[i*n + k] * B[k*o + j];
}

int check_convergence(const size_t& m, const double* A, const double& tol)
{
    for(size_t j = 0; j < m-1; ++j){//columns
	double tmp(0);
	for(size_t i = j+1; i < m; ++i)//rows
	    tmp += A[i*m + j];
	if(tmp > tol) // check inf-norm of columns for lower diagonal
	    return 1;
    }
    return 0;
}

void print_m_by_n(const size_t& m, const size_t& n, const double* A)
{

    for(size_t i = 0; i < m; ++i){
	for(size_t j = 0; j < n; ++j){
	    std::cout << A[n*i + j] << "  ";
	}
	std::cout << std::endl;
    }
}

void make_eye(const size_t& m, const size_t& n, double* A)
{
    for(size_t i = 0; i < m; ++i){
	for(size_t j = 0; j < n; ++j){
	    A[n*i + j] = 0.0;
	}
	A[n*i + i] = 1.0;
    }
}

} // namespace tools
