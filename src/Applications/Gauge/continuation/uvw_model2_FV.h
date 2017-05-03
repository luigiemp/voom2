#ifndef _UVW_FV_H
#define _UVW_FV_H 1
#include<gsl/gsl_vector.h>
#include"myfunctions.h"

extern int NOE;
extern double* nodes;

double W(double Phi,double beta,double m1,double m2);
double Wprime(double Phi,double beta,double m1,double m2);
double W2prime(double Phi,double beta,double m1,double m2);
double stiffness(double Phi,double m1,double m2,double k1,double k2);
double stiffness_prime(double Phi,double m1,double m2,double k1,double k2);

int firstvariation(double* indvar,double* arg,gsl_vector* out);
int local_firstvariation(double* indvar,double* arg,int element, double* out, double* nodes);
double Hermite(double xi, int index, int derv, double* nodes, int element);
double Lagrange(double xi, int index, int derv, double* nodes, int element);
int assign_data(double* data_in, int e, char var, double* data_out, double* nodes);
#endif


