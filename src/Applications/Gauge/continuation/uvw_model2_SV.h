#ifndef _UVW_SV_H
#define _UVW_SV_H 1
#include<gsl/gsl_vector.h>
#include"uvw_model2_FV.h"
#include"myfunctions.h"

int local_hessian(double* indvar,double* arg,int element, double** hessian, double* nodes);
int local_hessian_2constr(double* indvar,double* arg,int element, double** hessian, double* nodes);
int local_hessian_nonaxi(double* indvar,double* arg,int element, double hessian[11][11], double* nodes);
int Hessian(double* indvar,double* arg,gsl_matrix* out);
int Hessian_2constr(double* indvar,double* arg,gsl_matrix* out);
int Hessian_nonaxi(double* indvar,double* arg,gsl_matrix* out);

#endif


