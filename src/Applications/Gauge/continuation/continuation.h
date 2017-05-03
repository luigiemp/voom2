#ifndef _CONTINUATION_H
#define _CONTINUATION_H 1
#include<math.h>
#include<iostream>
#include<gsl/gsl_linalg.h>
using namespace std;
class Abstract{
	public:
		virtual int F(double* args,gsl_vector* out)=0;
		virtual int Jac(double* args,gsl_matrix* out)=0;
		virtual void array2vector(double*,gsl_vector* )=0;
		virtual void vector2array(gsl_vector*,double*)=0;

};
class Quadrature{
	public:
		void lebedevQuad(int (*fun)(double*,double*,gsl_vector*),double* arg,gsl_vector* out);
		void gaussLegendreQuad(int (*fun)(double*,double*,gsl_vector*),double* arg,gsl_vector* out);
		void adaptiveQuad(int (*fun)(double*,double*,gsl_vector*),double* arg, double a,\
				double b, double tau,gsl_vector* out);
		void gaussLegendreQuad_hess(int (*fun)(double*,double*,gsl_matrix*),double* arg,gsl_matrix* out);
		void adaptiveQuad_hess(int (*fun)(double*,double*,gsl_matrix*),double* arg, double a,\
				double b, double tau,gsl_matrix* out);
		void gaussLegendreQuad_hess_fou(int (*fun)(double*, double*, gsl_matrix*),double* args, gsl_matrix* out);
		
	public:
		Quadrature();
		int order;
		void computeQuad(double* args,double* out);
		
};
class LinearSolver{
	protected:
		double luSolve(gsl_matrix* A,gsl_vector* b,gsl_vector* x);
		double pseudoSolve(gsl_matrix* A, gsl_vector* b,gsl_vector* x);
	public:
		const char* solver;
		double linSolve(gsl_matrix* A,gsl_vector* b,gsl_vector* x);
};
class NonlinearSolver:public LinearSolver, public Abstract{
	protected:
		void quasiNewtonMethod(double* args);
		void newtonMethod(double* args);
		void array2vector(double*,gsl_vector* );
		void vector2array(gsl_vector*,double*);

	public:
		//Problem* problem;
		const char* predictorType;
		const char* correctorType;
		int noOfIterations;
		int maxIterations;
		double tolerance;
		double eps;
		gsl_vector* tangent;
		double* solution;
		gsl_matrix* jacobian;
		int noOfUnknownVar;
		int noOfEqs;
		int noOfInArgs;
		int (*ptr_F)(double*,double*,gsl_vector*);
		int noOfIndVar;
		double* indVar;
		
		NonlinearSolver();
		int setNoOfUnknownVar(int n);
		int setNoOfUnknownVar(int var, int eq);
		int setNoOfInArgs(int n);
		//void setNoOfIndVar(int n);
		//int setTangentSize(int n);
		//int (*ptr_F)(double* args,gsl_vector* out);
		int (*ptr_Jac)(double* args,gsl_matrix* out);
		int F(double* args,gsl_vector* out);
		int Jac(double* args,gsl_matrix* out);
		void predict(double* args);
		void correct(double* args);
		~NonlinearSolver();
		//int myFun(double* indvar,double* args,gsl_vector* out);
};
class Continuer:public NonlinearSolver{
	public:
		const char* continuerType;
		Quadrature Quad;
		int continuationPara;
		bool performQuad;
		double stepSize;
		double minStepSize;
		double maxStepSize;
		bool saveBPData;
		const char* BPFileName;
		const char* DataFileName;
		double LP;
		double BP;
		int posOfPara;
		int noOfPara;
		

		Continuer();
		//Problem* problem;
		int setNoOfUnknownVar(int n);
		int setNoOfUnknownVar(int var, int eq);
		int setNoOfInArgs(int n);
		void array2vector(double*,gsl_vector* );
		void vector2array(gsl_vector*,double*);
		void setInitialSolution(double* args);
		int F(double* args,gsl_vector* out);
		int Jac(double* args,gsl_matrix* out);
		char checkSingularPt(double* sol_tosv,gsl_vector* tgt_tosv);
		int switchBranch(double* sol,gsl_vector* tgt);
		int Continue(void);
		int load_BP(const char* filename,int BPnum);
		friend int post_process_sol(Continuer& c); //provides an interface for the user to implement a post processing routine
		~Continuer();
};
int savearray(const char* filename,double* array,int dim);
int loadarray(const char* filename,double* array,int dim);
int save_BP(const char* filename,double* array,gsl_vector*,int dim,int BPnum);
void pprintMat(gsl_matrix* M);
void printarray(double* A,int dim,char* format);
void printMatrix(gsl_matrix*);
/*
class EnergyProblem:public Problem{
	public:
		int r_len;
		int phi_len;
		int myFun(double*,double*,gsl_vector*);
};

double Plm(int l,int m,double x);
double Plm_t(int l,int m,double x);
double Plm_tt(int l,int m, double x);
double Ylm(int l,int m,double t,double p);
double Ylm_t(int l,int m,double t,double p);
double Ylm_tt(int l,int m,double t,double p);
double rad(double* c,double t,double p,int dim);
double rad_t(double* c,double t,double p,int dim);
double rad_tt(double* c,double t,double p,int dim);
double r_i(int i,double t,double p);
double rp_i(int i,double t,double p);
double rt_i(int i,double t,double p);
double rtt_i(int i,double t,double p);
double rtp_i(int i,double t,double p);
double rpp_i(int i,double t,double p);
double W(double Phi,double beta,double m1,double m2);
double Wprime(double Phi,double beta,double m1,double m2);
double W2prime(double Phi,double beta,double m1,double m2);
double stiffness(double Phi,double m1,double m2,double k1,double k2);
double stiffness_prime(double Phi,double m1,double m2,double k1,double k2);
*/
#endif
