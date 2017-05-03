#include"myfunctions.h"


double W(double Phi,double beta,double m1,double m2){
	return beta*pow((Phi-m1)*(Phi-m2),2);
}
double Wprime(double Phi,double beta,double m1,double m2){
	return 4*beta*(Phi-m1)*(Phi-m2)*(Phi-(m1+m2)/2.);
}
double W2prime(double Phi,double beta,double m1,double m2){
	return 4*beta*((Phi-m2)*(Phi-(m1+m2)/2.)+(Phi-m1)*(Phi-(m1+m2)/2.)+(Phi-m1)*(Phi-m2));
}
double stiffness(double Phi,double m1,double m2,double k1,double k2){
	return (( k2 - k1 )/( exp( (2* Phi-(m1+m2))/(2*0.1)) + 1 )+k1);
}
double stiffness_prime(double Phi,double m1,double m2,double k1,double k2){
	return (-(k2-k1)/pow(exp(1/2*(2*Phi-m1-m2)/.1)+1,2)/.1*exp(1/2*(2*Phi-m1-m2)/.1));
}
