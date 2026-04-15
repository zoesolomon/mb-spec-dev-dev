#ifndef MAT_TOOLS_H
#define MAT_TOOLS_H

namespace tools {

    double sgn(const double& );

void transpose(const size_t& , const size_t& ,
	       const double* , double* );

void invert(int, double* );

int eigenvalues(size_t, double*, double*);
int eigenvalues(size_t, double*, double*, double*);

double vector_2norm(const size_t& , const double* );

double vector_2norm(const size_t& , const size_t& , 
		    const size_t& , const double* );

void mat_mult(const size_t& , const double* , 
	      const size_t& , const double* , 
	      const size_t& , double* );

int check_convergence(const size_t& , const double* , const double& );

void print_m_by_n(const size_t&, const size_t&, const double*);

void make_eye(const size_t&, const size_t&, double*);

} // namespace tools

#endif // MAT_TOOLS_H
