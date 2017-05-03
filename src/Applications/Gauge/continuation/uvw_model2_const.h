#ifndef _UVW_CONST_H
#define _UVW_CONST_H 1
#include<gsl/gsl_vector.h>
#include"myfunctions.h"
extern int PLM_COUNT;
extern int u_len;
extern int v_len;
extern int w_len;
extern int phi_len;
extern double* plm;
extern double* plm_1;
extern double* plm_2;
extern double* plm_3;
extern double* plm_4;
extern double* SQFAC_i_1;
extern double* SQFAC_i_2;
extern double* SQFAC_i_3;
extern double* SQFAC_i_4;

double Plm(int l,int m,double x);
double Plm_t(int l,int m,double x);
double Plm_tt(int l,int m, double x);
double Ylm(int l,int m,double t,double p);
double Ylm_t(int l,int m,double t,double p);
double Ylm_tt(int l,int m,double t,double p);
/*
double rad(double* c,double t,double p,int dim);
double rad_t(double* c,double t,double p,int dim);
double rad_tt(double* c,double t,double p,int dim);
double r_i(int i,double t,double p);
double rp_i(int i,double t,double p);
double rt_i(int i,double t,double p);
double rtt_i(int i,double t,double p);
double rtp_i(int i,double t,double p);
double rpp_i(int i,double t,double p);
*/
double W(double Phi,double beta,double m1,double m2);
double Wprime(double Phi,double beta,double m1,double m2);
double W2prime(double Phi,double beta,double m1,double m2);
double stiffness(double Phi,double m1,double m2,double k1,double k2);
double stiffness_prime(double Phi,double m1,double m2,double k1,double k2);

int constraint(double* indvar,double* arg,gsl_matrix* out);
#endif


