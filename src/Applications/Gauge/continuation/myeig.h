#include "lapacke.h"
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_cblas.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_math.h>
#include<gsl/gsl_sf.h>
#include<gsl/gsl_eigen.h>

int lapack_gsl_sgeig_wrapper(gsl_matrix* A, gsl_matrix* B, int no_of_evals, gsl_vector* evals, gsl_matrix* evecs, int m){

  //check for square matrices

  // Array dimension
  int dim = A->size1; 
  int iwork[5*dim];
  double work[dim];
  return LAPACKE_dsygvx( LAPACK_ROW_MAJOR, 1, 'V', 'V', 'U', dim, \
		 A->data, dim, B->data, dim, -100, 100.,\
		 1, no_of_evals, 1e-10, &m ,evals->data ,evecs->data, \
			 dim, iwork);
	
}
