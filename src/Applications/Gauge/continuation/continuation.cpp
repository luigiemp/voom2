#include<iostream>
#include<gsl/gsl_linalg.h>
#include"continuation.h"
#include<gsl/gsl_cblas.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_math.h>
#include<gsl/gsl_sf.h>
#include<gsl/gsl_eigen.h>
#include<string.h>
#include<sstream>
#include"kbhit.h"
#include"gausslegendre.h"
//#include"arnoldi_eig.h"

#define PI (3.141592653589793)
#define QUAD 32
extern FILE* fid_toerase, *fid_toerase2;
using namespace std;
Quadrature::Quadrature(){
  order=101;
}
void Quadrature::lebedevQuad(int (*fun)(double*,double*,gsl_vector*), double* args,gsl_vector* out){
  string lebfile("Lebedev/OurLebedev");
  stringstream ss;
  ss<<order;
  lebfile=lebfile+ss.str()+".dat";
	
  FILE* pfile;
  double t,p,weight;
  double* tp=new double[2];
  gsl_vector* temp=gsl_vector_calloc(out->size);
  gsl_vector_set_zero(out);
  pfile=fopen(lebfile.c_str(),"r");
  //in OurLebedev files the first element is phi .in. [-180,180]
  // second element is theta .in. [0,180] and the third element is the weight
  while(fscanf(pfile,"%lf %lf %lf",&p,&t,&weight)!=EOF){
    if(t==0.0) t=1e-4;
    if(t==180.0) t=180.0-1e-4;
    tp[0]=t*PI/180;
    tp[1]=p*PI/180;
    fun(tp,args,temp);
    //out<--out+weight*temp
    //-------------should not include st term in the function--------------
    gsl_blas_daxpy(weight,temp,out);
  }
  gsl_blas_dscal(4*M_PI,out);
  gsl_vector_free(temp);
  delete[] tp;
  fclose(pfile);
};
void Quadrature::gaussLegendreQuad(int (*fun)(double*, double*, gsl_vector*),double* args, gsl_vector* out){
  gsl_vector* temp = gsl_vector_calloc(out->size); //temp =[0,0,...,0]_size
  gsl_vector_set_zero(out);
  /*
  int n = 1024/2; //THIS NEEDS TO BE MODIFIED
  double tp[2]={0.,0.};
  for (int i=0; i< n; i++){
    tp[0] = (1+ x1024[i])*PI/2; // in degrees
    fun(tp,args,temp);
    //out <-- out+weight*temp
    gsl_blas_daxpy(w1024[i]*M_PI/2,temp,out);
		
    //Now for PI-theta. The weights from PI-theta and theta are the same
    tp[0]=(1-x1024[i])*PI/2;
    fun(tp,args,temp);
    //out <-- out+weight*temp
    gsl_blas_daxpy(w1024[i]*M_PI/2,temp,out);
		
    }*/
    
    int n = 4/2; //THIS NEEDS TO BE MODIFIED
    double tp[2]={0.,0.};
    for (int i=0; i< n; i++){
    tp[0] = (1+ x4[i])*1.0/2; // in degrees
    fun(tp,args,temp);
    //out <-- out+weight*temp
    gsl_blas_daxpy(w4[i]*1.0/2,temp,out);
		
    //Now for PI-theta. The weights from PI-theta and theta are the same
    tp[0]=(1-x4[i])*1.0/2;
    fun(tp,args,temp);
    //out <-- out+weight*temp
    gsl_blas_daxpy(w4[i]*1.0/2,temp,out);
		
    }
  
    //gsl_blas_dscal(2*M_PI,out); //this step accounts for integration in phi. Very simple for axisymmetric version
  gsl_vector_free(temp);
};

//------------------- atha: quadratur for FEMS type hessian : domain [0,1] ---------
void Quadrature::gaussLegendreQuad_hess(int (*fun)(double*, double*, gsl_matrix*),double* args, gsl_matrix* out){
  gsl_matrix* temp = gsl_matrix_calloc(out->size1,out->size2); //temp =[0,0,...,0]_size
  gsl_matrix_scale(out,0); //resetting out =0
	
  int n = 4/2; //THIS NEEDS TO BE MODIFIED
  double tp[2]={0.,0.};
  for (int i=0; i< n; i++){
    tp[0] = (1+ x4[i])*1./2; // in degrees
    fun(tp,args,temp);
    //out <-- out+weight*temp
    gsl_matrix_scale(temp,w4[i]*1./2); //temp <--- weight*temp
    gsl_matrix_add(out,temp); //out <-- out+temp 
		
    //Now for PI-theta. The weights from PI-theta and theta are the same
    tp[0]=(1-x4[i])*1./2;
    fun(tp,args,temp);
    //out <-- out+weight*temp
    gsl_matrix_scale(temp,w4[i]*1./2); //temp <--- weight*temp
    gsl_matrix_add(out,temp); //out <-- out+temp 
		
  }
  //gsl_matrix_scale(out,2*PI); //this step accounts for integration in phi. Very simple for axisymmetric version
  gsl_matrix_free(temp);
};

//------------------- iti: quadratur for FEMS type hessian : domain [0,1] ---------
//--------------atha : Quad for fourier hessian : dom [0,pi] --------------
void Quadrature::gaussLegendreQuad_hess_fou(int (*fun)(double*, double*, gsl_matrix*),double* args, gsl_matrix* out){
	gsl_matrix* temp = gsl_matrix_calloc(out->size1,out->size2); //temp =[0,0,...,0]_size
	gsl_matrix_scale(out,0); //resetting out =0
	
	int n = 1024/2; //THIS NEEDS TO BE MODIFIED
	double tp[2]={0.,0.};
	
	for (int i=0; i< n; i++){
		tp[0] = (1+ x1024[i])*PI/2 ; // in degrees
		fun(tp,args,temp);
		//out <-- out+weight*temp
		gsl_matrix_scale(temp,w1024[i]*PI/2); //temp <--- weight*temp
		gsl_matrix_add(out,temp); //out <-- out+temp 
		
		//Now for PI-theta. The weights from PI-theta and theta are the same
		tp[0]=(1-x1024[i])*PI/2;
		fun(tp,args,temp);
		//out <-- out+weight*temp
		gsl_matrix_scale(temp,w1024[i]*PI/2); //temp <--- weight*temp
		gsl_matrix_add(out,temp); //out <-- out+temp 
		
	}
	
	//gsl_matrix_scale(out,2*PI); //this step accounts for integration in phi. Very simple for axisymmetric version
	gsl_matrix_free(temp);
};
//--------------iti : Quad for fourier hessian : dom [0,pi] --------------
void Quadrature::adaptiveQuad_hess(int (*fun)(double*, double*, gsl_matrix*),double* args, double a, double b, double tau, gsl_matrix* out){
  //printf("----Adaptive Quadrature-----\n");
  gsl_matrix* temp = gsl_matrix_calloc(out->size1,out->size2); //temp =[0,0,...,0]_size
  gsl_matrix_scale(out,0);
  gsl_matrix* out1 = gsl_matrix_calloc(out->size1,out->size2); //out1 =[0,0,...,0]_size
  // Lower order quad 
  int n = 12/2; 
  double tp[2]={0.,0.};
  gsl_matrix_scale(temp,0); //temp =0 
  for (int i=0; i< n; i++){
    tp[0] = ((a+b)/2.0 + (b-a)/2.0*x12[i]); // in degrees
    fun(tp,args,temp);
    //out1 <-- out1+weight*temp
		
    gsl_matrix_scale(temp,w12[i]*(b-a)/2); //temp <--- weight*temp
    gsl_matrix_add(out1,temp); //out1 <-- out1+temp 
		
    //Now for PI-theta. The weights from PI-theta and theta are the same
    tp[0]=((a+b)/2-x12[i]*(b-a)/2);
    fun(tp,args,temp);
    //out1 <-- out1+weight*temp
		
    gsl_matrix_scale(temp,w12[i]*(b-a)/2); //temp <--- weight*temp
    gsl_matrix_add(out1,temp); //out1 <-- out1+temp 
		
  }
  gsl_matrix_scale(out1,2*M_PI); //this step accounts for integration in phi. Very simple for axisymmetric version
  // Higher order quad 
  n = 16/2; 
  gsl_matrix_scale(temp,0);
  for (int i=0; i< n; i++){
    tp[0] = ((a+b)/2+ x16[i]*(b-a)/2); // in degrees
    fun(tp,args,temp);
    //out <-- out+weight*temp
		
    gsl_matrix_scale(temp,w16[i]*(b-a)/2); //temp <--- weight*temp
    gsl_matrix_add(out,temp); //out <-- out+temp 
		
    //Now for PI-theta. The weights from PI-theta and theta are the same
    tp[0]=((a+b)/2-x16[i]*(b-a)/2);
    fun(tp,args,temp);
    //out <-- out+weight*temp
    gsl_matrix_scale(temp,w16[i]*(b-a)/2); //temp <--- weight*temp
    gsl_matrix_add(out,temp); //out <-- out+temp 
		
  }
  gsl_matrix_scale(out,2*M_PI); //this step accounts for integration in phi. Very simple for axisymmetric version
  gsl_matrix_memcpy(temp,out1);
  gsl_matrix_sub(temp,out);
  //gsl_matrix_add(temp,out);
  double max_matrix = gsl_matrix_max(temp);
  double min_matrix = gsl_matrix_min(temp);
  double abs_max = max(abs(max_matrix),abs(min_matrix));
  //printf("max error = %e\n", tau);
  if(abs_max > tau && (b-a)> 1e-4 ){
    adaptiveQuad_hess(fun, args, a, (a+b)/2, tau/2, out1);
    adaptiveQuad_hess(fun, args, (a+b)/2, b, tau/2, out);
    //out<--out1+out
    gsl_matrix_add(out,out1);
  }
  gsl_matrix_free(temp);
  gsl_matrix_free(out1);
  //gsl_vector_free(temp);
  return;

};

void Quadrature::adaptiveQuad(int (*fun)(double*, double*, gsl_vector*),double* args, double a, double b, double tau, gsl_vector* out){
  gsl_vector* temp = gsl_vector_calloc(out->size); //temp =[0,0,...,0]_size
  gsl_vector_set_zero(out);
  // out1 - vector to store lower order quad
  // out - vector stores higher order quad
  gsl_vector* out1 = gsl_vector_calloc(out->size); //out1 =[0,0,...,0]_size
  // Lower order quad 
  int n = 12/2; 
  double tp[2]={0.,0.};
  gsl_vector_set_zero(temp);
  for (int i=0; i< n; i++){
    tp[0] = ((a+b)/2.0 + (b-a)/2.0*x12[i]); // in degrees
    fun(tp,args,temp);
    //out <-- out+weight*temp
    gsl_blas_daxpy(w12[i]*(b-a)/2,temp,out1);
		
    //Now for PI-theta. The weights from PI-theta and theta are the same
    tp[0]=((a+b)/2-x12[i]*(b-a)/2);
    fun(tp,args,temp);
    //out <-- out+weight*temp
    gsl_blas_daxpy(w12[i]*(b-a)/2,temp,out1);
		
  }
  gsl_blas_dscal(2*M_PI,out1); //this step accounts for integration in phi. Very simple for axisymmetric version
  // Higher order quad 
  n = 16/2; 
  gsl_vector_set_zero(temp);
  for (int i=0; i< n; i++){
    tp[0] = ((a+b)/2+ x16[i]*(b-a)/2); // in degrees
    fun(tp,args,temp);
    //out <-- out+weight*temp
    gsl_blas_daxpy(w16[i]*(b-a)/2,temp,out);
		
    //Now for PI-theta. The weights from PI-theta and theta are the same
    tp[0]=((a+b)/2-x16[i]*(b-a)/2);
    fun(tp,args,temp);
    //out <-- out+weight*temp
    gsl_blas_daxpy(w16[i]*(b-a)/2,temp,out);
		
  }
  gsl_blas_dscal(2*M_PI,out); //this step accounts for integration in phi. Very simple for axisymmetric version
  gsl_vector_memcpy(temp,out1);
  gsl_blas_daxpy(-1,out,temp);
  if(gsl_blas_dnrm2(temp) > tau && (b-a)> 1e-4 ){
    adaptiveQuad(fun, args, a, (a+b)/2, tau/2, out1);
    adaptiveQuad(fun, args, (a+b)/2, b, tau/2, out);
    //out<--out1+out
    gsl_blas_daxpy(1.0,out1,out);
  }
  gsl_vector_free(temp);
  gsl_vector_free(out1);
  //gsl_vector_free(temp);
  return;

};
//* SOLVER*/
double LinearSolver::luSolve(gsl_matrix* A, gsl_vector* b,gsl_vector* x){
  //I don't want to change the matrix A when I leave this fun.
  gsl_matrix* temp=gsl_matrix_alloc(A->size1,A->size2);
  gsl_matrix_memcpy(temp,A);
	
	
  gsl_permutation* p=gsl_permutation_alloc(b->size);
  int s;
  gsl_linalg_LU_decomp(temp,p,&s); //After this line A is an "LU-ed" matrix
  gsl_linalg_LU_solve(temp,p,b,x);
		
  /* ****************SV SOLVE*****************************/
  /*
    gsl_matrix* V=gsl_matrix_alloc(A->size2,A->size2);
    gsl_vector* work=gsl_vector_alloc(A->size2);
    gsl_vector* S=gsl_vector_alloc(A->size2);
    gsl_linalg_SV_decomp(A,V,S,work);
    gsl_linalg_SV_solve(A,V,S,b,x);
  */
  //double determinant=gsl_linalg_LU_det(A,s);
  //gsl_matrix_memcpy(A,temp);
  gsl_matrix_free(temp);
  return 0.;
}
double LinearSolver::pseudoSolve(gsl_matrix* A, gsl_vector* b,gsl_vector* x){
  /* If A is a MxN matrix (M>=N), then we can use QR factorization to solve
     the matrix equation Ax=b in the least square sense. It turns out 
     ( http://tutorial.math.lamar.edu/Classes/LinAlg/LeastSquares.aspx ) that
     if A has rank N (# of cols) then the least square solution solves A'Ax=A'b
  */
  gsl_matrix* temp = gsl_matrix_calloc(A->size1,A->size2);
  gsl_matrix_memcpy(temp,A);
  gsl_vector* tau = gsl_vector_calloc(min(A->size1,A->size2));
  gsl_vector* residual = gsl_vector_calloc(b->size);

  gsl_linalg_QR_decomp(temp,tau);
  gsl_linalg_QR_lssolve(temp,tau,b,x,residual);
  cout << "Residual Norm = " << gsl_blas_dnrm2(residual)<<endl;
  gsl_vector_free(residual);
  gsl_matrix_free(temp);
  gsl_vector_free(tau);


}
double LinearSolver::linSolve(gsl_matrix* A,gsl_vector* b,gsl_vector* x){
  if(strcmp(solver,"LU")==0) {
    //cout << "Using LU decomposition for Inverse" <<endl;
    return luSolve(A,b,x);
  }
  else if(strcmp(solver,"Pseudo")==0) {
    //cout << "Using Pseudo Inverse" << endl;
    return pseudoSolve(A,b,x);
  }
  else {
    cout<<"Unkown Solver: "<<solver<<endl;
    return NULL;
  }
}
NonlinearSolver::NonlinearSolver(){
  predictorType="generic";
  correctorType="quasi-newton";
  solver="LU";
  maxIterations=50;
  noOfIterations=0;
  tolerance=1e-5;
  eps=1e-6;
  indVar=NULL;
	
}
int NonlinearSolver::setNoOfUnknownVar(int n){
  noOfUnknownVar=n;
  noOfEqs=noOfUnknownVar; //default: # equations = # unknowns
  jacobian=gsl_matrix_calloc(n,n);
  gsl_matrix_set_identity(jacobian);
  tangent=gsl_vector_calloc(n);
  return 0;
}
int NonlinearSolver::setNoOfUnknownVar(int var, int eq){
  noOfUnknownVar=var;
  noOfEqs=eq; 
  jacobian=gsl_matrix_calloc(eq,var);
  //gsl_matrix_set_identity(jacobian);
  tangent=gsl_vector_calloc(var);
  return 0;
}
int NonlinearSolver::setNoOfInArgs(int n){
  noOfInArgs=n;
  solution=new double[n];
  return 0;
}
NonlinearSolver::~NonlinearSolver(){
  delete[] solution;
  gsl_vector_free(tangent);
  gsl_matrix_free(jacobian);
}
/*
  void NonlinearSolver::setNoOfIndVar(int n){
  noOfIndVar=2;
  indVar[0]=0.;
  indVar[1]=0.;
  }
*/
int NonlinearSolver::F(double* args, gsl_vector* out){
  //cout<<"In Nonlinear::F "<<endl;
  //copying args into a new variable in
  double* in=new double[noOfInArgs];
  memcpy(in,args,sizeof(double)*noOfInArgs);
  ptr_F(indVar,in,out);
  delete[] in;
  return 0;
}
int NonlinearSolver::Jac(double* args, gsl_matrix* out){
  //copying args into a new variable in
  //cout<<"In Nonlinear:: Jac "<<endl;
  double* in=new double[noOfInArgs];
  memcpy(in,args,sizeof(double)*noOfInArgs);
	
  int row=(int)out->size1;
  int col=(int)out->size2;
  gsl_vector* temp_v=gsl_vector_calloc(row);
  gsl_vector* f0=gsl_vector_calloc(row);
  double temp_x;
  //problem->myFun(problem->indVar,args,f0);
  //cout<<"In Nonlinear::Jac\n";
  F(args,f0);
  for(int i=0;i<col;i++){
    temp_x=in[i];
    in[i]=in[i]+eps;
    F(in,temp_v);
    //temp_v<---temp_v-f0
    gsl_blas_daxpy(-1,f0,temp_v);
    //temp_v<--(temp_v)/eps;
    gsl_blas_dscal(1/eps,temp_v);
    //set the ith column as (Fpert-F0)/eps
    gsl_matrix_set_col(out,i,temp_v);
    //reset the perturbation
    in[i]=temp_x;
		
  }
  gsl_vector_free(temp_v);
  gsl_vector_free(f0);
  //gsl_matrix_memcpy(out,out_temp);
  //gsl_matrix_free(out_temp);
  delete[] in;
  return 0;
}
void NonlinearSolver::newtonMethod(double* args){
  //cout<<"Entering Newton Method\n";
  double* x_np1=new double[noOfInArgs];
  gsl_vector* v_np1=gsl_vector_calloc(jacobian->size2);
  gsl_vector* F_np1=gsl_vector_calloc(jacobian->size1);
  gsl_vector* Dv=gsl_vector_calloc(jacobian->size2);
  //gsl_matrix* jac=gsl_matrix_calloc(jacobian->size1,jacobian->size2);
  memcpy(x_np1,args,sizeof(double)*noOfInArgs);
  array2vector(args,v_np1);
  F(args,F_np1);
  double nrm=gsl_blas_dnrm2(F_np1);
  if(gsl_isnan(nrm)) {cout<<"Nrm returned is nan!"<<endl;exit(-1);}
  int ITERATIONS=0;
  while(nrm>tolerance ||ITERATIONS > 20){
    Jac(x_np1,jacobian);
    linSolve(jacobian,F_np1,Dv);
    //v_np1=v_np1-Dv
    //the -ve sign is because "Dv"= - inv(jac)*F_np1
    gsl_blas_daxpy(-1.0,Dv,v_np1);
    vector2array(v_np1,x_np1);
    F(x_np1,F_np1);
    nrm=gsl_blas_dnrm2(F_np1);
    if(gsl_isnan(nrm)) {cout<<"Norm in the interation was a nan!"<<endl; exit(-1);}
    cout<<"Iterating Norm= "<<nrm<<endl;
  }
  memcpy(solution,x_np1,sizeof(double)*noOfInArgs);
	
  gsl_vector_free(v_np1);
  gsl_vector_free(F_np1);
  gsl_vector_free(Dv);
  delete[] x_np1;
  //gsl_matrix_free(jac);
}

/*
 * Implementing the quasi newton scheme using the SR1 method. 
 * See wikipedia entry on quasi newton scheme for the method.
 * Actually, see the wiki entry on Broyden.
  */

//---------- Atha: Broyden without inverting the Hessian -------
void NonlinearSolver::quasiNewtonMethod(double* Args){
  cout<<"beginning lin corrector stage"<<endl;	
  //First, I don't want to change Args:
  double* args=new double[noOfInArgs];
  memcpy(args,Args,sizeof(double)*noOfInArgs);
  noOfIterations=0;
  gsl_vector* x_k=gsl_vector_alloc(jacobian->size1);
  gsl_vector* F_k=gsl_vector_alloc(jacobian->size1);
  gsl_vector* y_k=gsl_vector_alloc(jacobian->size1);
  gsl_vector* Dx_k=gsl_vector_alloc(jacobian->size1);
  gsl_matrix* B_k=gsl_matrix_alloc(jacobian->size1,jacobian->size2);
  	
  double denominator;
  gsl_vector* temp_x=gsl_vector_alloc(jacobian->size1);
  gsl_vector* temp_F=gsl_vector_alloc(jacobian->size1);

  Jac(args,B_k); 
  cout<<"end jacobian"<<endl;	
  gsl_matrix_memcpy(jacobian,B_k);
    
  array2vector(args,x_k); F(args,F_k);
  
  int ITERATIONS=0;
  double nrm=gsl_blas_dnrm2(F_k);
  if(gsl_isnan(nrm)) {printf("Failed: even before the loop!"); exit(EXIT_FAILURE);}
  while(nrm>tolerance){
    if(ITERATIONS > 30|| nrm>1e6) {printf("Failed to converge. Exiting!"); exit(EXIT_FAILURE);}
    ITERATIONS++;
    cout<<"in the loop"<<endl;
    //Dx_k=-inv(B_k)*F_k, here alpha_k=1
    luSolve(B_k,F_k,Dx_k);
    gsl_blas_dscal(-1.0,Dx_k);
    //store x_k in temp_k
    gsl_vector_memcpy(temp_x,x_k);
    //x_k<--x_k+Dx_k
    gsl_vector_add(x_k,Dx_k); //x_k+1
    //  y_k=F_k+1-F_k 
    //temp_k<---F_k
    gsl_vector_memcpy(temp_F,F_k);
    //Compute F_k+1
    vector2array(x_k,args); F(args,F_k); //F_k+1, args="x_k+1"
    //F_k<--F_k-F_k+1
    gsl_blas_daxpy(-1.0,F_k,temp_F);
    //F_k<-- -F_k
    gsl_blas_dscal(-1.0,temp_F);
    //y_k<-- F_k
    gsl_vector_memcpy(y_k,temp_F);
    //Updating B i.e, computing B_k+1
    //temp_x<--: -B_k*Dx_k 
    gsl_blas_dgemv(CblasNoTrans,-1.0,B_k,Dx_k,0.,temp_x);
    //y_k<--: y_k - B_kDx_k
    gsl_blas_daxpy(1.0,temp_x,y_k);
    //denominator<--Dx_k^T * Dx_k
    gsl_blas_ddot(Dx_k,Dx_k,&denominator);
    //a rank one update to B_k
    gsl_blas_dger(1/denominator,y_k,Dx_k,B_k);
    nrm=gsl_blas_dnrm2(F_k);
    if(gsl_isnan(nrm)) {printf("Failed to converge due to nan!\n"); exit(EXIT_FAILURE);}
    cout<<"Iterating Norm= "<<nrm<<endl;
  }
  memcpy(solution,args,sizeof(double)*noOfInArgs);
  gsl_vector_free(x_k);
  gsl_vector_free(F_k);
  gsl_vector_free(y_k);
  gsl_vector_free(Dx_k);
  gsl_vector_free(temp_x);
  gsl_vector_free(temp_F);
  gsl_matrix_free(B_k);
  delete[] args;
  cout<<"end corrector"<<endl;	
}
//---------- Iti: Broyden without inverting the hessian -----
/*
//---------Atha: Broyden with inversion: Slower----------
void NonlinearSolver::quasiNewtonMethod(double* Args){
  cout<<"beginning lin corrector stage"<<endl;	
  //First, I don't want to change Args:
  double* args=new double[noOfInArgs];
  memcpy(args,Args,sizeof(double)*noOfInArgs);
  noOfIterations=0;
  gsl_vector* x_k=gsl_vector_alloc(jacobian->size1);
  gsl_vector* F_k=gsl_vector_alloc(jacobian->size1);
  gsl_vector* y_k=gsl_vector_alloc(jacobian->size1);
  gsl_vector* Dx_k=gsl_vector_alloc(jacobian->size1);
  gsl_matrix* H_k=gsl_matrix_alloc(jacobian->size1,jacobian->size2);
  gsl_matrix* H_temp=gsl_matrix_alloc(jacobian->size1,jacobian->size2);
  gsl_permutation* p=gsl_permutation_alloc(jacobian->size1);
  int s;
	
  double denominator;
  double numerator;
  gsl_vector* temp_x=gsl_vector_alloc(jacobian->size1);
  gsl_vector* temp_F=gsl_vector_alloc(jacobian->size1);

  Jac(args,H_temp); 
  cout<<"end jacobian"<<endl;	
  gsl_matrix_memcpy(jacobian,H_temp);
  //H_k from now is actually the inverse of the jacobian
  gsl_linalg_LU_decomp(H_temp,p,&s);
  gsl_linalg_LU_invert(H_temp,p,H_k);
  //gsl_matrix_set_identity(H_k);
	
  array2vector(args,x_k); F(args,F_k);
  //gsl_vector_fprintf(stdout,F_k,"%f");cout<<endl;
  int ITERATIONS=0;
  double nrm=0;
  while(gsl_blas_dnrm2(F_k)>tolerance){
    if(ITERATIONS > 500|| nrm>1e10) {printf("Failed to converge. Exiting!"); exit(EXIT_FAILURE);}
    ITERATIONS++;
    cout<<"in the loop"<<endl;
    //Dx_k=-alpha_k*H_k*F_k, here alpha_k=1
    gsl_blas_dgemv(CblasNoTrans,-1.0,H_k,F_k,0.,Dx_k);
    //store x_k in temp_k
    gsl_vector_memcpy(temp_x,x_k);
    //x_k<--x_k+Dx_k
    gsl_vector_add(x_k,Dx_k); //x_k+1
    //  y_k=F_k+1-F_k
    //temp_k<---F_k
    gsl_vector_memcpy(temp_F,F_k);
    //Compute F_k+1
    vector2array(x_k,args); F(args,F_k); //F_k+1, args="x_k+1"
    //F_k<--F_k-F_k+1
    gsl_blas_daxpy(-1.0,F_k,temp_F);
    //F_k<-- -F_k
    gsl_blas_dscal(-1.0,temp_F);
    //y_k<-- F_k
    gsl_vector_memcpy(y_k,temp_F);
    //Updating H i.e, computing H_k+1
    //temp_x<--H_k*y_k 
    gsl_blas_dgemv(CblasNoTrans,1.0,H_k,y_k,0.,temp_x);
    //denominator<--Dx_k^T * H_k *y_k
    gsl_blas_ddot(Dx_k,temp_x,&denominator);
    //First saving Dx_k in temp_x
    gsl_vector_memcpy(temp_x,Dx_k);
    //Dx_k<--Dx_k-H_k*y_k
    gsl_blas_dgemv(CblasNoTrans,-1.0,H_k,y_k,1.0,Dx_k);
    //numerator<--(Dx_k-H_k*y_k)*Dx_k^T
    gsl_blas_ddot(Dx_k,temp_x,&numerator);
    //H_k<---(1+numerator/denominator)*H_k
    gsl_matrix_scale(H_k,1+numerator/denominator);
    nrm=gsl_blas_dnrm2(F_k);
    cout<<"Iterating Norm= "<<nrm<<endl;
  }
  memcpy(solution,args,sizeof(double)*noOfInArgs);
  gsl_vector_free(x_k);
  gsl_vector_free(F_k);
  gsl_vector_free(y_k);
  gsl_vector_free(Dx_k);
  gsl_vector_free(temp_x);
  gsl_vector_free(temp_F);
  gsl_matrix_free(H_k);
  gsl_matrix_free(H_temp);
  gsl_permutation_free(p);
  delete[] args;
  cout<<"end corrector"<<endl;	
}
//-----------ithi : Broyden with inversion---------------
*/
void NonlinearSolver::correct(double* args){
  if(strcmp(correctorType,"quasi-newton")==0)
    quasiNewtonMethod(args);
  else if(strcmp(correctorType,"newton")==0)
    newtonMethod(args);
  else
    cout<<"Unknown Corrector!"<<endl;
}
void NonlinearSolver::predict(double* args){
  if(strcmp(predictorType,"generic")==0){
    //------------------- atha: usual predictor ----------------
    
    cout<<"entering predictor stage"<<endl;
    //input: args, outputs: jacobian
    Jac(args, jacobian); 
    cout<<"computed jacobian"<<endl;
    //jacobian->size1=no. of rows
    gsl_vector* vecPredict=gsl_vector_calloc(jacobian->size1);
		
    //vecPredict=[0,0,0,...,1]
		
    gsl_vector_set(vecPredict,(int)(jacobian->size1-1), 1.0);
    memcpy(solution,args,sizeof(double)*noOfInArgs);
    cout<<"beginning lin solve"<<endl;
    linSolve(jacobian,vecPredict,tangent);
    cout<<"done lin solve"<<endl;
    // If the tangent points in the opposite direction of increasign arc
    //   length, then reverse the direction
    
    //
    //  if(gsl_vector_get(tangent,tangent->size-1)<0)
    //gsl_blas_dscal(-1,tangent);
    //
    //Normalize the tangent. Since the pseudo-arc equation in F
    //is tangent_(n-1).tangent_n = 1, tangent_n will not
    // have a norm =1.  
    
    gsl_vector_fprintf(stdout,tangent,"%f");
    
    gsl_blas_dscal(1/gsl_blas_dnrm2(tangent),tangent);
    gsl_vector_free(vecPredict);
    
  }
  else{
    cout<<"Predictor Type : "<<predictorType<<" not available"<<endl;
  }
}
void NonlinearSolver::array2vector(double* args,gsl_vector* out){
  //cout<<"In Nonlinear::array2vector\n";
  if(out->size!=noOfUnknownVar) 
    cout<<"No. Of variables not equal to size of vector!"<<endl;
  for(int i=0;i<noOfUnknownVar;i++)
    gsl_vector_set(out,i,args[i]);
  //gsl_vector_set(out,noOfUnknownVar,args[posOfPara]);
}
void NonlinearSolver::vector2array(gsl_vector* vecIn, double* out){
  if(vecIn->size!=noOfUnknownVar)
    cout<<"Size of vector not equal to no. Of variables!"<<endl;
  for(int i=0;i<noOfUnknownVar;i++)
    out[i]=gsl_vector_get(vecIn,i);
  //out[posOfPara]=gsl_vector_get(vecIn,noOfUnknownVar);
}


/* CONTINUER*/
Continuer::Continuer(){
  performQuad=true;
  stepSize=0.01;
  minStepSize=1e-4;
  maxStepSize=10.;
  saveBPData=true;
  BPFileName="BP_Filename";
  DataFileName="Data_Filename.txt";
  continuerType="pseudo-arc";
  LP=0;
  BP=0;
	
}
void Continuer::array2vector(double* args,gsl_vector* out){
  //cout<<"In Continuer::array2vector\n";
  if(out->size!=noOfUnknownVar+1) 
    cout<<"No. Of variables not equal to size of vector!"<<endl;
  for(int i=0;i<noOfUnknownVar;i++)
    gsl_vector_set(out,i,args[i]);
  gsl_vector_set(out,noOfUnknownVar,args[posOfPara]);
}
void Continuer::vector2array(gsl_vector* vecIn, double* out){
  if(vecIn->size!=noOfUnknownVar+1)
    cout<<"Size of vector not equal to no. Of variables!"<<endl;
  for(int i=0;i<noOfUnknownVar;i++)
    out[i]=gsl_vector_get(vecIn,i);
  out[posOfPara]=gsl_vector_get(vecIn,noOfUnknownVar);
}
int Continuer::F(double* args,gsl_vector* out){
  //cout<< "In Continuer::F\n";
  //if(continuerType=="pseudo-arc" && performQuad==true){
  gsl_vector* v1=gsl_vector_calloc(noOfUnknownVar+1);
  gsl_vector* v2=gsl_vector_calloc(noOfUnknownVar+1);
  array2vector(args,v1);
  array2vector(solution,v2);
  //v1<--v1-v2
		
  gsl_blas_daxpy(-1.,v2,v1);
  //Pseudo-arc equation
  double v1DotTgt;
  gsl_blas_ddot(v1,tangent,&v1DotTgt);
  double N=v1DotTgt-stepSize;
  //temporarily set last equation(value) to zero so that
  //while doin the quadrature in the next step we are not
  //accumulating anything
  out->data[noOfEqs]=0.;
  // NOW I NEED TO PERFROM QUADRATURE
  if(performQuad==true){
    //Quad.lebedevQuad(ptr_F,args,out);
    Quad.gaussLegendreQuad(ptr_F,args,out);
    //Quad.adaptiveQuad(ptr_F, args, 0.0, 1.0, 1e-6, out);

  }
  else{
    ptr_F(NULL,args,out);
  }
  //set the last equation(value) to be the pseudo-arc equation
  out->data[noOfEqs]=N;
  gsl_vector_free(v1);
  gsl_vector_free(v2);
	
}

//-------- Atha: Center Difference------
int Continuer::Jac(double* args,gsl_matrix* out){
  // ---------- Central difference Jacobian -----------
  // This is more accurate than the forward difference
  // implemented in the (commented) code below
  //cout<<"In Continuer::Jac\n"<<noOfUnknownVar<<endl;
  int row=(int)out->size1;
  int col=(int)out->size2;
  gsl_vector* temp_v=gsl_vector_calloc(row);
  gsl_vector* temp_u=gsl_vector_calloc(row);
  
  double temp_x;
  //Part of the Jacobian by varying unknowns
  for(int i=0;i<noOfUnknownVar;i++){
    temp_x=args[i];
    args[i]=args[i]+eps;
    F(args,temp_v);
    //2, since we have added eps to the variable two lines above
    args[i]=args[i]-2*eps; 
    F(args,temp_u);
    //temp_v<---temp_v-temp_u
    gsl_blas_daxpy(-1,temp_u,temp_v);
    //temp_v<--(temp_v)/(2 eps);
    gsl_blas_dscal(1/(2.*eps),temp_v);
    //set the ith column as (Fpert-F0)/eps
    gsl_matrix_set_col(out,i,temp_v);
    //reset the perturbation
    args[i]=temp_x;
    //cout<<"HERE "<<i<<endl;
  }
  //Part of the jacobian by varying parameter
  temp_x=args[posOfPara];
  args[posOfPara]=args[posOfPara]+eps;
  F(args,temp_v);
  //2, since we have added eps to the variable two lines above
  args[posOfPara]=args[posOfPara]-2.*eps;
  F(args,temp_u);
  //temp_v<---temp_v-f0
  gsl_blas_daxpy(-1,temp_u,temp_v);
  //temp_v<--(temp_v)/eps;
  gsl_blas_dscal(1/(2.*eps),temp_v);
  //set the next, i.e, "noOfUnknownVar"th position
  gsl_matrix_set_col(out,noOfUnknownVar,temp_v);
  //reset the perturbation
  args[posOfPara]=temp_x;
	
  //The last row of th Jacobian contains the tangent vector
  gsl_matrix_set_row(out,out->size1-1,tangent);
	
  gsl_vector_free(temp_v);
  gsl_vector_free(temp_u);
  return 0;
	
  //return NonlinearSolver::Jac(args,out);
}
//------Iti : Center Difference
/*
//-----------Atha : Forward Difference--------
int Continuer::Jac(double* args,gsl_matrix* out){
  //cout<<"In Continuer::Jac\n"<<noOfUnknownVar<<endl;
  int row=(int)out->size1;
  int col=(int)out->size2;
  gsl_vector* temp_v=gsl_vector_calloc(row);
  gsl_vector* f0=gsl_vector_calloc(row);
  double temp_x;
  F(args,f0);
  //Part of the Jacobian by varying unknowns
  for(int i=0;i<noOfUnknownVar;i++){
    temp_x=args[i];
    args[i]=args[i]+eps;
    F(args,temp_v);
    //temp_v<---temp_v-f0
    gsl_blas_daxpy(-1,f0,temp_v);
    //temp_v<--(temp_v)/eps;
    gsl_blas_dscal(1/eps,temp_v);
    //set the ith column as (Fpert-F0)/eps
    gsl_matrix_set_col(out,i,temp_v);
    //reset the perturbation
    args[i]=temp_x;
    //cout<<"HERE "<<i<<endl;
  }
  //Part of the jacobian by varying parameter
  temp_x=args[posOfPara];
  args[posOfPara]=args[posOfPara]+eps;
  F(args,temp_v);
  //temp_v<---temp_v-f0
  gsl_blas_daxpy(-1,f0,temp_v);
  //temp_v<--(temp_v)/eps;
  gsl_blas_dscal(1/eps,temp_v);
  //set the next, i.e, "noOfUnknownVar"th position
  gsl_matrix_set_col(out,noOfUnknownVar,temp_v);
  //reset the perturbation
  args[posOfPara]=temp_x;
	
  //The last row of th Jacobian contains the tangent vector
  gsl_matrix_set_row(out,out->size1-1,tangent);
	
  gsl_vector_free(temp_v);
  gsl_vector_free(f0);
  return 0;
	
  //return NonlinearSolver::Jac(args,out);
}
//--------Iti : forward difference -----
*/
int Continuer::setNoOfUnknownVar(int n){
  //+1 to account for the continuation parameter
  double ret=NonlinearSolver::setNoOfUnknownVar(n+1);
  //resets noOfUnknown variables fron n+1 to n
  noOfUnknownVar=n;
  return ret;
}
int Continuer::setNoOfUnknownVar(int var, int eq){
  //+1 to account for the continuation parameter
  double ret=NonlinearSolver::setNoOfUnknownVar(var+1,eq+1);
  //resets noOfUnknown variables fron n+1 to n
  noOfUnknownVar=var;
  noOfEqs=eq;
  return ret;
}
int Continuer::setNoOfInArgs(int n){
  return NonlinearSolver::setNoOfInArgs(n);
}
char Continuer::checkSingularPt(double* solution_toSave,gsl_vector* tgt_BP){
  gsl_matrix* jac_b;
  if (jacobian->size1 != jacobian->size2){ 
    cout << "Non symmetric matrix!"<<endl;    
    // Construct a sym matrix with size = # columns of jacobian
    jac_b=gsl_matrix_calloc(jacobian->size2,jacobian->size2);
    gsl_matrix* jac_bb = gsl_matrix_calloc(jac_b->size1, jac_b->size2);
    

    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,jacobian,jacobian,0.0,jac_b);
    //exit(0);
    gsl_matrix* V=gsl_matrix_alloc(jac_b->size1,jac_b->size2);
    gsl_vector* S=gsl_vector_alloc(jac_b->size2);
    gsl_vector* work=gsl_vector_alloc(jac_b->size2);
    //gsl_vector* tgt_BP = gsl_vector_alloc(V->size1);

    gsl_matrix_memcpy(jac_bb,jac_b);
    gsl_linalg_SV_decomp(jac_bb,V,S,work);	
    gsl_matrix_transpose(V);
    gsl_matrix_get_row(tgt_BP,V,V->size1-1);
    //gsl_vector_fprintf(stdout,S,"%le");	
    //cout << endl;
    //gsl_vector_fprintf(stdout,tgt, "%le");
    //pprintMat(c.jacobian);
    if (S->data[S->size-1] <1e-4) {
      
      gsl_vector* sol_temp=gsl_vector_alloc(tgt_BP->size);
      array2vector(solution,sol_temp);
      //sol<--sol+stepSize*tgt_BP
      //gsl_blas_daxpy(stepSize,tgt_BP,sol_temp);
      vector2array(sol_temp,solution_toSave);
		
      //save_BP(BPFileName,solution_toSave,tgt_BP,noOfInArgs);
			
      gsl_vector_free(sol_temp);
      gsl_vector_free(work);
      gsl_vector_free(S);
      gsl_matrix_free(V);
      gsl_matrix_free(jac_b);
      gsl_matrix_free(jac_bb);
      return 'B';
    }
  }
  else{
    jac_b=gsl_matrix_calloc(jacobian->size1,jacobian->size2);
    gsl_matrix_memcpy(jac_b,jacobian);
  }
    gsl_matrix_view temp=gsl_matrix_submatrix(jac_b,0,0,jac_b->size1-1,jac_b->size2-1);
    gsl_matrix* jac_l=gsl_matrix_calloc(jac_b->size1-1,jac_b->size2-1);
	
    gsl_matrix_memcpy(jac_l,&temp.matrix);
	
    gsl_permutation* p_b=gsl_permutation_alloc(jac_b->size1);
    int s;
    cout<<"computing LU decomp\n";
    gsl_linalg_LU_decomp(jac_b,p_b,&s);
    cout<<"computed LU decomp\n";
    double BP_next=(double)gsl_linalg_LU_det(jac_b,s);
	
    gsl_permutation* p_l=gsl_permutation_alloc(jac_l->size1);
    gsl_linalg_LU_decomp(jac_l,p_l,&s);
    double LP_next=(double)gsl_linalg_LU_det(jac_l,s);
    printf("\033[0;32m %le %le\n \033[m",LP_next,BP_next);
    //cout << "here-------"<<1e600<<endl;
    //cout << LP_next << " " << BP_next <<endl;
    if(LP==0 && BP==0){
      //Which is the case for the very first step
      LP=LP_next;
      BP=BP_next;
		
      gsl_matrix_free(jac_l);
      gsl_matrix_free(jac_b);
      gsl_permutation_free(p_b);
      gsl_permutation_free(p_l);
		
      return 'f';
    }
    else if(LP*LP_next<0 && BP*BP_next<0){
      //Which is the case when we move past a BP
      cout<<"\033[1;33m Branch Pt\033[m Detected"<<endl;
      LP=LP_next;
      BP=BP_next;
		
      //Matrix and Vector definitions needed for SVD
      gsl_matrix* V=gsl_matrix_alloc(jacobian->size2,jacobian->size2);
      gsl_vector* S=gsl_vector_alloc(jacobian->size2);
      gsl_vector* work=gsl_vector_alloc(jacobian->size2);
      gsl_matrix_memcpy(jac_b,jacobian);
      //jac_b=USV^t; U is stored in jac_b
      gsl_linalg_SV_decomp(jac_b,V,S,work);
      /*Places the last row of the matrix V^T in tgt_BP (gives 
	the bifurcation direction*/
      gsl_vector_fprintf(stdout,S,"%le");	
      gsl_matrix_transpose(V);
      gsl_matrix_get_row(tgt_BP,V,V->size1-1); //<--------------------------------
		
      //Now, we find the solution
      gsl_vector* sol_temp=gsl_vector_alloc(tgt_BP->size);
      array2vector(solution,sol_temp);
      //sol<--sol+stepSize*tgt_BP
      //gsl_blas_daxpy(stepSize,tgt_BP,sol_temp);
      vector2array(sol_temp,solution_toSave);
		
      //save_BP(BPFileName,solution_toSave,tgt_BP,noOfInArgs);
			
      gsl_vector_free(sol_temp);
      gsl_vector_free(work);
      gsl_vector_free(S);
      gsl_matrix_free(V);
      gsl_matrix_free(jac_l);
      gsl_matrix_free(jac_b);
      gsl_permutation_free(p_b);
      gsl_permutation_free(p_l);
		
      return 'B';
    }
    else if(LP*LP_next<0 && BP*BP_next>0){
      //Which happens at a LP
      cout<<"\033[1;34m Limit Pt\033[m Detected"<<endl;
      LP=LP_next;
      BP=BP_next;
      gsl_matrix_free(jac_l);
      gsl_matrix_free(jac_b);
      gsl_permutation_free(p_b);
      gsl_permutation_free(p_l);
      return 'L';
    }
    else{
      LP=LP_next;
      BP=BP_next;
      gsl_matrix_free(jac_b);
      gsl_permutation_free(p_b);
      gsl_permutation_free(p_l);
      return '-';
    }
  

}
void Continuer::setInitialSolution(double* args){
  //arraycpy(args,solution,noOfInArgs);
	
  memcpy(solution,args,sizeof(double)*noOfInArgs);
  //gsl_vector_set_zero(tangent);
  cout<<"JACOBIAN"<<endl;
  Jac(args,jacobian);
  cout<<"JACOBIAN DONE"<<endl;
  /*The jacobian will, at this stage, always have 0's in the last row because
    a trivial Augmented Pseudo arc equation. We will remove the last row to 
    find the tangential direction to move
  */
  //gsl_matrix_transpose(jacobian);
  gsl_matrix_view jac_clip_view=gsl_matrix_submatrix(jacobian,0,0,jacobian->size1,jacobian->size2);
  gsl_matrix* jac_clip=gsl_matrix_alloc(jacobian->size1,jacobian->size2);
  gsl_matrix_memcpy(jac_clip,&jac_clip_view.matrix);
  /*Now, performing an SVD to find the tangential direction*/
  /*------------ Matrices, vectors needed for svd------------*/
	
  gsl_matrix* U=gsl_matrix_alloc(jac_clip->size2,jac_clip->size2);
  gsl_vector* S=gsl_vector_alloc(jac_clip->size2);
  gsl_vector* work=gsl_vector_alloc(jac_clip->size2);
  
  //--------- SVD to get the tangent vector --------
  gsl_linalg_SV_decomp(jac_clip,U,S,work);
  gsl_matrix_get_row(tangent,jac_clip,jac_clip->size2-1);

  //gsl_vector_fprintf(stdout,tangent, "%le");
  
  // //--------- arnoldi method to get the zero eigen direction--------
  // gsl_matrix* mass = gsl_matrix_calloc(jac_clip->size2,jac_clip->size2);
  // gsl_vector* eig=gsl_vector_calloc(jac_clip->size1);
  // gsl_matrix* evec = gsl_matrix_calloc(jac_clip->size2,jac_clip->size2);
  // gsl_matrix_set_identity(mass);
  // arnoldi_eig(jac_clip,mass,10,eig,evec,1e-20);
  //  for(int i=0; i<10;i++) printf("%le\n",eig->data[i]);
  // gsl_matrix_get_col(tangent,evec,0);
  // //somehow the arnoldi eig gives a vector pointing in the "wrong direction"
  // //so, multiplying it by the last element of the vector, to point in the inc 's' direction
  // gsl_blas_dscal(tangent->data[tangent->size-1],tangent); 
  
  //------------iti: arnodli
  
  cout<<"SVD done"<<endl;
  gsl_matrix_free(U);
  gsl_vector_free(S);
  gsl_vector_free(work);
  gsl_matrix_free(jac_clip);
  //gsl_matrix_transpose(jacobian);
	
}
int Continuer::Continue(){
  gsl_vector* tgt_BP=gsl_vector_calloc(tangent->size);
  double* solution_toSave=new double[noOfInArgs];
  memcpy(solution_toSave,solution,sizeof(double)*noOfInArgs);
  int BPnum=0;
  int key;
  do{
    if(kbhit())
      {
	key = fgetc(stdin);
	if (key == 'i'){
	  stepSize += .01;
	  cout<< "\n \033[1;35m Step Size increased to "<< stepSize<<"\033[m\n\n";
	}
	else if(key=='d'){
	  stepSize -= .01;
	  cout<< "\n \033[1;36m Step Size decreased to "<< stepSize<<"\033[m\n\n";
	}
	else if(key=='r'){
	  stepSize=-stepSize;
	  cout<< "\n \033[1;36m Step Size reversed \033[m\n\n";
	}
	else if(key =='n'){
	  stepSize = 2*stepSize;
	  cout<< "\n \033[1;36m Step Size doubled to"<<stepSize<<" \033[m\n\n";
	}
	else if(key =='h'){
	  stepSize = stepSize/2.0;
	  cout<< "\n \033[1;33m Step Size halved to"<<stepSize<<" \033[m\n\n";
	}
	fflush(stdin);
      } 
    else{
			
      predict(solution);
      //setInitialSolution(solution);//use the svd way to predict
      correct(solution);
      printf("------------\n");
      
      char sing_pt=checkSingularPt(solution_toSave,tgt_BP);
      
      // -------- RESET jacobian to it old value--------//
      
      cout<<solution[posOfPara]<<" "<<sing_pt<<endl;
      		
      if(sing_pt=='L'){ 
	
      }
      else if(sing_pt=='B'){
	BPnum=BPnum+1;
	
	save_BP(BPFileName,solution_toSave,tgt_BP,noOfInArgs,BPnum);
      }
    cout<<"saving...\n";
    savearray(DataFileName,solution,noOfInArgs);
    cout<<"saved.\n";
    }
     //post_process_sol(*this);
  }while(key!='x');
};
int Continuer::load_BP(const char* filename,int BPnum){
  string fname(filename);
  stringstream out;
  out<<BPnum;
  fname=fname+"_"+out.str()+".txt";
  string fnameT(filename);
  fnameT=fnameT+"_"+out.str()+"_Tgt"+".txt";
  loadarray(fname.c_str(),solution,noOfInArgs);
	
  cout<<endl;
  FILE* pfile;
  pfile=fopen(fnameT.c_str(),"r");
  gsl_vector_fscanf(pfile,tangent);
	
  fclose(pfile);
  gsl_vector* sol_temp=gsl_vector_alloc(tangent->size);
  //solution=solution+stepSize*tangent
  array2vector(solution,sol_temp);
  gsl_blas_daxpy(stepSize,tangent,sol_temp);
  vector2array(sol_temp,solution);
  gsl_vector_free(sol_temp);
  return 0;
}
Continuer::~Continuer(){
}
int save_BP(const char* filename,double* array,gsl_vector* tgt,int dim,int BPnum){
  string fname(filename);
  stringstream out;
  out<<BPnum;
  fname=fname+"_"+out.str()+".txt";
  savearray(fname.c_str(),array,dim);
  FILE* pfile;

  string fnameT(filename);
  fnameT=fnameT+"_"+out.str()+"_Tgt"+".txt";
  pfile=fopen(fnameT.c_str(),"a");
  gsl_vector_fprintf(pfile,tgt,"%le");
  fclose(pfile);
}
int savearray(const char* filename,double* array,int dim){
  FILE* pfile;
  pfile=fopen(filename,"a");
  for(int i=0;i<dim;i++)
    fprintf(pfile,"%1.20f ", array[i]);
  fprintf(pfile,"\n");
  fclose(pfile);
  return 0;
}
int loadarray(const char* filename,double* out,int dim){
  //double* array=new double[dim];
  FILE* pfile;
  pfile=fopen(filename,"r");
  if (pfile!=NULL){
    for(int i=0;i<dim;i++){
      double temp;
      int retval=fscanf(pfile,"%le ", &temp);
      out[i]=temp;
    }
    fclose(pfile);
    //memcpy(out,array,dim);
    //printarray(out,dim,"%f|~~");
	
  }
  else{
    cout<<"File: "<<filename<<" does not exist"<<endl;
  }
  //delete[] array;
  return 0;
}
void pprintMat(gsl_matrix* M){
  cout<<"["<<endl;
  double val;
  for(int i=0;i<M->size1;i++){
    cout<<"[";
    for(int j=0;j<M->size2;j++){
      val=gsl_matrix_get(M,i,j);
      if(fabs(val)<1e-2){ 
	printf("%1.15e, ",val);
      }
      else{
	//printf("\033[1;33m%1.2e \033[m ",val);
	printf("%1.15e, ",val);
      }
    }
    cout<<"],";
  }
  cout<<"]";
}
void printarray(double* A,int dim,char* format){
  for(int i=0;i<dim;i++)
    printf(format,A[i]);
  cout<<endl;
}

