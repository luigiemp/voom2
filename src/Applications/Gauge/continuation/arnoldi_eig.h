#ifndef _ARNOLDI_H
#define _ARNOLDI_H 1
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_cblas.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_math.h>
#include<gsl/gsl_sf.h>
#include<gsl/gsl_eigen.h>
#include<stdio.h>
template<class ARFLOAT>
void arpack_matrix_wrapper( ARFLOAT* &A, ARFLOAT* &B, gsl_matrix* a, gsl_matrix* b);
int arnoldi_eig(gsl_matrix* a, gsl_matrix* mass, int no_of_eigs, gsl_vector* eigs, gsl_matrix* evecs, double sigma);
#endif
