#ifndef _MYFUNCS_H
#define _MYFUNCS_H 1
double W(double Phi,double beta,double m1,double m2);
double Wprime(double Phi,double beta,double m1,double m2);
double W2prime(double Phi,double beta,double m1,double m2);
double stiffness(double Phi,double m1,double m2,double k1,double k2);
double stiffness_prime(double Phi,double m1,double m2,double k1,double k2);
#include<iostream>
#include<gsl/gsl_math.h>
#include<gsl/gsl_sf.h>
#define PI (3.141592653589793)
#endif
