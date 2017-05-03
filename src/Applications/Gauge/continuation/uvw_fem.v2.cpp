#include<iostream>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_eigen.h>
#include<gsl/gsl_cblas.h>
#include<gsl/gsl_blas.h>
#include<stdio.h>
#include <unistd.h>
#include<string.h>
#include"uvw_model2_FV.h"
#include"uvw_model2_SV.h"
#include"uvw_model2_const.h"
#include"continuation.h"
//#include"myeig.h"
#include"arnoldi_eig.h"
#define PI (3.141592653589793)
#define NO_OF_EIGS 10
void print_help();
double* nodes;
int NOE=100;

int main(int argc,char** argv){
  /* Parsing Command line options */
  char* in_filename=NULL;   //Atmost 50 characters
  char* out_filename="outfile.txt";
  int opt;
  int DataPt=1;
  char* cont_para = "null";  //a string to store continuation para;
  int para;                //corresponding int value
  double strt_step= 0.01;
  bool flag_trivial=false;
  bool flag_profile=false;
  bool flag_switch=false;
  bool flag_cont_para=false;
  bool flag_hessian=false;
  int vec_pos=0;
  
  while ( (opt=getopt(argc,argv,"tsckfohDSHpVM"))!=-1){
    switch (opt){
    case 'M':
      NOE = atoi(argv[optind]);
      break;
    case 't': //Trivial solution
      flag_trivial=true;
      break;
    case 's': //Switch branch
      flag_switch=true;
      break;
    case 'f': //set filename
      in_filename=argv[optind];
      break;
    case 'c':  //Stifness cont
      flag_cont_para=true;
      cont_para=argv[optind];
      break;
    case 'D': //sets the data point the value of the argument
      DataPt=atoi(argv[optind]);
      break;
    case 'S': //Step length
      //printf("%s\n",optarg);
      
      strt_step=atof(argv[optind]);
      printf("Step Size: %f\n",strt_step);
      break;
    case 'o':
      out_filename=argv[optind];
      break;
    case 'p':
      flag_profile=true;
      break;
    case 'H':
      flag_hessian=true;
      break;
    case 'V':
      vec_pos=atoi(argv[optind]);
      break;
    case 'h':
      print_help();
    }
  }
  /* Parse cont_para to find which parameter we need to continue in */
  if (strcmp(cont_para,"eps")==0)  
    para=0;
  else if (strcmp(cont_para,"stiff1")==0)
    para=1;
  else if (strcmp(cont_para,"stiff2")==0)
    para=2;
  else if (strcmp(cont_para,"beta")==0)
    para=3;
  else if (strcmp(cont_para,"p")==0)
    para=4;
  else if (strcmp(cont_para,"mu")==0)
    para=5;
  
  int Tot_len= 6*NOE;

  Continuer c;
  c.ptr_F=firstvariation;
  c.setNoOfUnknownVar(Tot_len+4, Tot_len+4); //var:=Tot_len + no of Lag mult, eqs:=Tot_len(normal eq)+Total_len(tgt_equ) -phi_len + 4(constraints) 
  c.setNoOfInArgs(Tot_len+10);   // Tot_len + no.Of Lag. mult + other para
  //c.noOfEqs = 2*Tot_len + 4; 
  c.noOfIndVar=2;  //Corresponds to theta and phi 
  c.noOfPara=6;   // l(~1/eps), k1, k2, beta, pressure, mu
  c.correctorType="quasi-newton"; //"quasi-newton"
  c.solver="LU";
  c.Quad.order=65;
  //----------- GENERATING ELEMENT NODES -----------------
  nodes = new double[NOE+1];
  for (int n=0; n<= NOE; n++){
    //nodes[n] = n*PI/NOE; // UNIFORM
    //nodes[n] =n*n*n*PI/(NOE*NOE*NOE);
    nodes[n] = PI*pow(sin(n*PI/(2.*NOE)),2);
    //cout << n << " " << nodes[n]*180/PI << endl;
  }
 
  if(flag_trivial==true){
    c.posOfPara=Tot_len+4+0;  //0 Corresponds to l (~1/eps), 4: p
    /* Setting up the trivial solution*/
    double k1=1;
    double k2=0;
    double beta=1;
    double pressure=1;//-20.2975505;
    double mu=0.0;
    double w=W(mu,beta,-1,1);
    double wprime=Wprime(mu,beta,-1,1);
    double D=stiffness(mu,-1,1,k1,k1+k2);
    double D_prime=stiffness_prime(mu,-1,1,k1,k1+k2);
    double Lag_phi=-(wprime+D_prime*4.);
    double gamma=-beta*(mu*mu-1)*(mu*mu-1)+pressure/2.0;
    double Lag_z=0;
    double Lag_gz=0;
    /* -----------Starting Pt for Continuation -----------------*/
    double l= 0.98;
    
    /* **************Step Size ***************/
    c.stepSize=strt_step;
    c.BPFileName="primary_BP_eps";
    c.DataFileName=out_filename;
    
    double* trivial_sol=new double[Tot_len+10];
    for(int i=0;i<Tot_len;i++) trivial_sol[i]= 0.0;
    //for (int i=0; i<2*(NOE-1); i+=2) trivial_sol[i]=0.0;
    //trivial_sol[r_len]=2*mu*sqrt(PI); //trivial solution mu=constant
    trivial_sol[Tot_len]=  gamma;   
    trivial_sol[Tot_len+1]= Lag_phi;
    trivial_sol[Tot_len+2]= Lag_z;
    trivial_sol[Tot_len+3] = Lag_gz;  
    trivial_sol[Tot_len+4]= l;
    trivial_sol[Tot_len+5]= k1;
    trivial_sol[Tot_len+6]= k2;
    trivial_sol[Tot_len+7]= beta;
    trivial_sol[Tot_len+8]= pressure;
    trivial_sol[Tot_len+9]= mu;
    
    //trivial_sol[0]=0;
    //trivial_sol[5]=1e-5;
    /*
    
      gsl_vector* vec =gsl_vector_calloc(Tot_len+5);
      //c.F(trivial_sol,vec);
      
      //printf("b=[");
      //gsl_vector_fprintf(stdout,vec,"%le");
      //printf("];");
    
    
      c.setInitialSolution(trivial_sol);
      c.Jac(trivial_sol,c.jacobian);
      cout<<"import numpy"<<endl;
      cout<<"A=numpy.matrix("<<endl;
      pprintMat(c.jacobian);
      printf(");");
      gsl_matrix* V=gsl_matrix_alloc(c.jacobian->size2,c.jacobian->size2);
      gsl_vector* S=gsl_vector_alloc(c.jacobian->size2);
      gsl_vector* work=gsl_vector_alloc(c.jacobian->size2);
      //gsl_linalg_SV_decomp(c.jacobian,V,S,work);
      
      //printf("c=[");
      //gsl_vector_fprintf(stdout,S,"%le");
      //printf("];");   
    
    
      gsl_vector_free(work);
      gsl_vector_free(S);
      gsl_matrix_free(V);
      
      /*
      c.setInitialSolution(trivial_sol);
      //Matrix and Vector definitions needed for SVD
      gsl_matrix* jac_b=gsl_matrix_calloc(c.jacobian->size2,c.jacobian->size2);
      gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,c.jacobian,c.jacobian,0.0,jac_b);
   
      gsl_matrix* V=gsl_matrix_alloc(jac_b->size1,jac_b->size2);
      gsl_vector* S=gsl_vector_alloc(jac_b->size2);
      gsl_vector* work=gsl_vector_alloc(jac_b->size2);
      gsl_vector* tgt = gsl_vector_alloc(V->size1);

    
      gsl_linalg_SV_decomp(jac_b,V,S,work);	
      gsl_matrix_transpose(V);
      gsl_matrix_get_row(tgt,V,V->size1-1);
      gsl_vector_fprintf(stdout,S,"%le");	
      cout << endl;
      gsl_vector_fprintf(stdout,tgt, "%le");
      pprintMat(c.jacobian);
      gsl_vector_free(work);
      gsl_vector_free(S);
      gsl_matrix_free(V);
    */
    
    
    c.setInitialSolution(trivial_sol);
    c.Continue();
    return 0;
  }
  if(flag_switch==true){
    c.posOfPara=Tot_len+4+0;  //0 Corresponds to l (~1/eps), 4: pressure
    
    c.stepSize=strt_step;
    c.BPFileName="secondary_BP_eps";
    c.DataFileName=out_filename;
    c.load_BP("primary_BP_eps",1);
    c.Continue();
  }
  if(flag_cont_para==true){
    printf("Continuing in %s (%d) \n",cont_para,para);
    printf("Using data from %s \n",in_filename);
    c.posOfPara=Tot_len+4+para;    //3 - for three lagrange multipliers
    /* **************Step Size ***************/
    c.stepSize=strt_step;
    c.BPFileName="DecK_BP";
    c.DataFileName=out_filename;
    
    double* sol=new double[c.noOfInArgs];
    
    
    
    FILE* pfile;
    pfile=fopen(in_filename,"r");
    if (pfile!=NULL){
      for(int j=0;j<DataPt;j++){
	for(int i=0;i<c.noOfInArgs;i++){
	  double temp;
	  int retval=fscanf(pfile,"%le ", &temp);
	  sol[i]=temp;
	}
      }
    }
    fclose(pfile);
    
    c.setInitialSolution(sol);
    
    c.Continue();
    
  }

  // ----------------- HESSIAN -----------------------
  // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  
  if(flag_hessian==true){
    
    printf("Using data from %s \n",in_filename);
    printf("Along the parameter %s (%d) \n",cont_para,para);
    c.posOfPara=Tot_len+4+para;    //3 - for three lagrange multipliers
    /* **************Step Size ***************/
    c.DataFileName=out_filename;
    
    double* sol=new double[c.noOfInArgs];
    
    FILE* pfile;
    pfile=fopen(in_filename,"r");
    if (pfile!=NULL){
      for(int j=0;j<DataPt;j++){
	for(int i=0;i<c.noOfInArgs;i++){
	  double temp;
	  int retval=fscanf(pfile,"%le ", &temp);
	  sol[i]=temp;
	}
      }
    }
    else {printf("File does not exist\n");}
    printf("l = %f\n", sol[Tot_len+4]);
    fclose(pfile);
    
    int u_len = 2*(NOE);
    int phi_len = 2*(NOE);
    
    int constr=1;
    int no_of_constr=3;
    if (constr==1){
      //-----computing the hessian and constraints at the given data point 
      //--------------- AXISYMMETRIC ---------------------------
      // The basic hessian for nonaxi is u_len+phi_len+2, since the boundary
      // terms suppressed in axisym should be included. The number of rows
      // are increased by a further 2 to account for the constraints
      gsl_matrix* hess_constr = gsl_matrix_calloc(u_len+phi_len+no_of_constr,u_len+phi_len);
      gsl_matrix* mass_constr = gsl_matrix_calloc(u_len+phi_len+no_of_constr,u_len+phi_len);


      c.Quad.gaussLegendreQuad_hess(Hessian,sol,hess_constr);
      //c.Quad.adaptiveQuad_hess(Hessian,sol,0, PI, 1e-6, hess_constr);

      sol[Tot_len+7]=-1; // flag for mass computation (corresponds to beta)
      c.Quad.gaussLegendreQuad_hess(Hessian,sol,mass_constr);
      
      gsl_matrix_view hess = gsl_matrix_submatrix(hess_constr,0,0,u_len+phi_len,u_len+phi_len);
      gsl_matrix_view mass = gsl_matrix_submatrix(mass_constr,0,0,u_len+phi_len,u_len+phi_len);
      gsl_vector_view constr_1 = gsl_matrix_row(hess_constr,u_len+phi_len);
      gsl_vector_view constr_2 = gsl_matrix_row(hess_constr,u_len+phi_len+1);
      gsl_vector_view constr_3 = gsl_matrix_row(hess_constr,u_len+phi_len+2);
      
      /*
      //---------------- NON AXI SYMMETRIC ---------------------------
      // The basic hessian for nonaxi is u_len+phi_len+4, since the boundary
      // terms suppressed in axisym should be included. The number of rows
      // are increased by a further 2+2(w,psi) to account for the constraints
      gsl_matrix* hess_constr = gsl_matrix_calloc(u_len+phi_len+4+2,u_len+phi_len+4);
      gsl_matrix* mass_constr = gsl_matrix_calloc(u_len+phi_len+4+2,u_len+phi_len+4);
    
      //---------------------Two possible quadratures ----------------------------
      c.Quad.gaussLegendreQuad_hess(Hessian_nonaxi,sol,hess_constr);
    
      sol[Tot_len+7]=-1; // flag for mass computation (corresponds to beta)
    
      c.Quad.gaussLegendreQuad_hess(Hessian_nonaxi,sol,mass_constr);
    
      // Constructing the (projected) hessian matrix on the constraint space----
      gsl_matrix_view hess = gsl_matrix_submatrix(hess_constr,0,0,u_len+phi_len+4,u_len+phi_len+4);
      gsl_matrix_view mass = gsl_matrix_submatrix(mass_constr,0,0,u_len+phi_len+4,u_len+phi_len+4);
      gsl_vector_view constr_1 = gsl_matrix_row(hess_constr,u_len+phi_len+4);
      gsl_vector_view constr_2 = gsl_matrix_row(hess_constr,u_len+phi_len+5);
      gsl_vector_view constr_3 = gsl_matrix_row(hess_constr,u_len+phi_len+6);
    
      //----------------------------------------------------------------------------
      */
      //EigenValues using Arnoldi Iteration (Gives fast and accurate values)
      gsl_vector* eig=gsl_vector_calloc((&hess.matrix)->size1);
      gsl_matrix* evec=gsl_matrix_calloc((&hess.matrix)->size1,(&hess.matrix)->size2);
      arnoldi_eig(&hess.matrix,&mass.matrix, NO_OF_EIGS, eig, evec, 0);
      
      int no_of_neg_eigs=0;
      for(int i=0; i<NO_OF_EIGS;i++){
	//the relevant eigen values that determine stability
	if(eig->data[i]<0 && fabs(eig->data[i])>1e-6) no_of_neg_eigs++;
	printf("%le\n",eig->data[i]);
      
      }
      cout<< "No of significant negative eigen values: " <<no_of_neg_eigs<<endl;
      printf("--------Inner product of eigen vectors with constraints-----------------\n");
      double overlap;
      //gsl_vector_fprintf(stdout,&eigen_vector.vector,"%e"); cout << "------- iti : EVector --------"<<endl;

      for(int choice=0; choice< no_of_neg_eigs; choice++){
	cout<<"~~~~~~~ For Eigen vector "<<choice+1<<" ~~~~~~~~~~~"<<endl;
	gsl_vector_view eigen_vector = gsl_matrix_column(evec, choice);
	gsl_blas_ddot(&constr_1.vector,&eigen_vector.vector,&overlap);
	printf("<C1,ev> = %le\n",overlap);
	gsl_blas_ddot(&constr_2.vector,&eigen_vector.vector,&overlap);
	printf("<C2,ev> = %le\n",overlap);
	gsl_blas_ddot(&constr_3.vector,&eigen_vector.vector,&overlap);
	printf("<C3,ev> = %le\n",overlap);
    
      }

      //pprintMat(evec);
      //gsl_vector_fprintf(stdout,&eigen_vector.vector,"%f");
    
      /*
      // Check if Volume Constraint is satisfied 
      gsl_matrix* vol_const = gsl_matrix_calloc(1,1);
      for(int i=0;i<u_len-1;i++){
      // insert the mode into 'sol' at phi's location
      sol[i+u_len+v_len]=(&eigen_vector.vector)->data[i];
      
      }
    
      // Integrate 'vol_const' (1x1 matrix) over the domain
      c.Quad.gaussLegendreQuad_hess(constraint,sol,vol_const);
      //c.Quad.adaptiveQuad_hess(constraint,sol,0, PI, 1e-6, vol_const);		
      pprintMat(vol_const);
      */
      gsl_vector_free(eig);

      gsl_matrix_free(evec);
      gsl_matrix_free(hess_constr);
    
    }
    else{
    
      /*  
      //----------------------- FOURIER ------------------
      gsl_matrix* hess_constr = gsl_matrix_calloc(w_len+psi_len,w_len+psi_len);
      gsl_matrix* mass = gsl_matrix_calloc(w_len+psi_len,w_len+psi_len);

      //-----no of cols of proj =  dim of projected vector (hence smaller)
    
      c.Quad.gaussLegendreQuad_hess_fou(Hessian_fou,sol,hess_constr);
      //c.Quad.adaptiveQuad_hess(Hessian_fou,sol,0, PI, 1e-6, hess_constr);
      cout<<"import numpy"<<endl;
      cout<<"A=numpy.matrix("<<endl;
      pprintMat(hess_constr);
      cout<<")\n";
      sol[Tot_len+7]=-1; // flag for mass computation (corresponds to beta)
    
      c.Quad.gaussLegendreQuad_hess_fou(Hessian_fou,sol,mass);
      // c.Quad.adaptiveQuad_hess(Hessian_fou,sol,0, PI, 1e-6, mass);
      cout<<"M=numpy.matrix("<<endl;
      pprintMat(mass);
      cout<<")";
      //gsl_matrix_set_identity(mass);
      */
    
      
      //---------------------  Axi symmetric ------------------------
    
      gsl_matrix* hess_constr = gsl_matrix_calloc(u_len+phi_len,u_len+phi_len);
      gsl_matrix* mass = gsl_matrix_calloc(u_len+phi_len,u_len+phi_len);

      //-----no of cols of proj =  dim of projected vector (hence smaller)
    
   
      c.Quad.gaussLegendreQuad_hess(Hessian,sol,hess_constr);
    
      cout<<"import numpy"<<endl;
      cout<<"A=numpy.matrix("<<endl;
      //pprintMat(hess_constr);
      cout<<")\n";
      sol[Tot_len+7]=-1; // flag for mass computation (corresponds to beta)
    
      c.Quad.gaussLegendreQuad_hess(Hessian,sol,mass);
    
      cout<<"M=numpy.matrix("<<endl;
      //pprintMat(mass);
      cout<<")";
      //gsl_matrix_set_identity(mass);
    
    
      /*
      //--------------------- Non axisymmetric ---------------------
      gsl_matrix* hess_constr = gsl_matrix_calloc(u_len+phi_len+2,u_len+phi_len+2);
      gsl_matrix* mass = gsl_matrix_calloc(u_len+phi_len+2,u_len+phi_len+2);
    
      c.Quad.gaussLegendreQuad_hess(Hessian_nonaxi,sol,hess_constr);
      cout<<"import numpy"<<endl;
      cout<<"A=numpy.matrix("<<endl;
      //pprintMat(hess_constr);
      cout<<")";
      sol[Tot_len+7]=-1; // flag for mass computation (corresponds to beta)
     
      c.Quad.gaussLegendreQuad_hess(Hessian_nonaxi,sol,mass);
      //gsl_matrix_set_identity(mass);
      */
      gsl_vector* eig=gsl_vector_calloc(hess_constr->size1);
      gsl_matrix* evec=gsl_matrix_calloc(hess_constr->size1,hess_constr->size2);
      gsl_eigen_gensymmv_workspace* work=gsl_eigen_gensymmv_alloc(hess_constr->size1);
    
      /*
      // ---------------- Atha: GSL_eig solver----------------
      gsl_eigen_gensymmv(hess_constr,mass,eig,evec,work);
      gsl_eigen_gensymmv_sort(eig,evec,GSL_EIGEN_SORT_VAL_DESC);
      gsl_vector_view eigen_vector = gsl_matrix_column(evec, u_len+phi_len-1-vec_pos);

      //gsl_vector_fprintf(stdout,eig,"%e");
    
    
      //cout << "----"<<endl; gsl_vector_fprintf(stdout,&eigen_vector.vector,"%f");
      //-------------------Iti : GSL eig solver--------
      */
   
      /*
      //--------------- Atha: LAPACK eigen value ----------
      int m;
      lapack_gsl_sgeig_wrapper(hess_constr,mass, 20,  eig,  evec, m);
      gsl_vector_fprintf(stdout,eig,"%le");
      //--------------- Iti : LAPACK eigen value
      */
      arnoldi_eig(hess_constr,mass, NO_OF_EIGS, eig, evec, 0);
      gsl_vector_view eigen_vector = gsl_matrix_column(evec, 0);

      for(int i=0; i<NO_OF_EIGS;i++) printf("%le\n",eig->data[i]);
      //cout << "----"<<endl; gsl_vector_fprintf(stdout,&eigen_vector.vector,"%f");
      gsl_eigen_gensymmv_free(work);//gsl_eigen_symmv_free(work);
      gsl_vector_free(eig);
    
      gsl_matrix_free(evec);
      gsl_matrix_free(hess_constr);
    }
    return 0;
    
    
  }

  if(flag_profile==true){
    c.posOfPara=Tot_len+4;  //Corresponds to l (~1/eps)
    /* Setting up the trivial solution*/
    double k1=0.8;
    double k2=0.8;
    double beta=1;
    double pressure=1.0;
    double mu=0.0;
    double w=W(mu,beta,-1,1);
    double wprime=Wprime(mu,beta,-1,1);
    double D=stiffness(mu,-1,1,k1,k2);
    double D_prime=stiffness_prime(mu,-1,1,k1,k2);
    double Lag_phi=-(wprime+D_prime*4.);
    double gamma=-beta*(mu*mu-1)*(mu*mu-1)-Lag_phi*mu+pressure/2.0;
    double Lag_z=0;
    double Lag_gz=0.0;
    /* -----------Starting Pt for Continuation -----------------*/
    double l=0.9;
    //c.posOfPara=Tot_len+3;  //Corresponds to l (~1/eps)
    
    /* **************Step Size ***************/
    c.stepSize=strt_step;
    c.BPFileName="primary_BP_eps";
    c.DataFileName="temp_out";
    
    double* trivial_sol=new double[Tot_len+10];
    for(int i=0;i<Tot_len;i++) trivial_sol[i]= 0;
    trivial_sol[Tot_len]=  gamma;   
    trivial_sol[Tot_len+1]= Lag_phi;
    trivial_sol[Tot_len+2]= Lag_z;  
    trivial_sol[Tot_len+3]=Lag_gz;
    trivial_sol[Tot_len+4]= l;
    trivial_sol[Tot_len+5]= k1;
    trivial_sol[Tot_len+6]= k2;
    trivial_sol[Tot_len+7]= beta;
    trivial_sol[Tot_len+8]= pressure;
    trivial_sol[Tot_len+9]= mu;
    gsl_vector* vec=gsl_vector_alloc(Tot_len+5);
    
    //c.F(trivial_sol,vec);
    c.Jac(trivial_sol,c.jacobian);
    //gsl_vector_fprintf(stdout, vec, "%le");
    
    gsl_vector* tgt=gsl_vector_calloc(Tot_len+5);
    tgt->data[Tot_len+4]=1.0;
    c.tangent->data=tgt->data;
    //c.setInitialSolution(trivial_sol);
    //pprintMat(c.jacobian);
    gsl_vector_free(tgt);
    
    gsl_vector_free(vec);
    
    /*
    //Matrix and Vector definitions needed for SVD
    gsl_matrix* V=gsl_matrix_alloc(c.jacobian->size2,c.jacobian->size2);
    gsl_vector* S=gsl_vector_alloc(c.jacobian->size2);
    gsl_vector* work=gsl_vector_alloc(c.jacobian->size2);
    
    gsl_linalg_SV_decomp(c.jacobian,V,S,work);	
    gsl_vector_fprintf(stdout,S,"%le");	
    pprintMat(V);
    gsl_vector_free(work);
    gsl_vector_free(S);
    gsl_matrix_free(V);
    */
    
    //EigenValues
    /*
      gsl_eigen_symmv_workspace* work=gsl_eigen_symmv_alloc(c.jacobian->size1);
      gsl_vector* eig=gsl_vector_calloc(c.jacobian->size1);
      gsl_matrix* evec=gsl_matrix_calloc(c.jacobian->size1,c.jacobian->size2);
      gsl_eigen_symmv(c.jacobian,eig,evec,work);
      gsl_eigen_symmv_sort(eig,evec,GSL_EIGEN_SORT_VAL_DESC);
      gsl_vector_fprintf(stdout,eig,"%le");
      gsl_vector_free(eig);
      gsl_eigen_symmv_free(work);
      gsl_matrix_free(evec);
    */
    delete[] trivial_sol;
    return 0;
  }
}
void print_help(){
  printf("\033[1;33m t\033[m : Continue in 1/epsilon land along the trivial solution.\n");
  printf("\033[1;33m s\033[m : Switch branch.\n");
  printf("\033[1;33m f\033[m : Takes an argument \033[4;36mfilename\033[m and sets the input filename.\n");
  printf("\033[1;33m S\033[m : Sets the \033[4;36mstep length\033[m.\n");
  printf("\033[1;33m D\033[m : Sets the Data Pt with an integer.\n");
  printf("\033[1;33m c\033[m : Continues along a \033[4;36parameter\033[m.\n");
  printf("\033[1;33m o\033[m : Sets the \033[4;36moutput filename\033[m.\n");
  printf("\033[1;33m p\033[m : Profiler\n");
  printf("\033[1;33m H\033[m : Hessian\n");
  printf("\033[1;33m V\033[m : Select eigen vector postition (from the end) for vol constraint\n");
  printf("\033[1;33m h\033[m : Prints this help.\n");
  
}
int post_process_sol(Continuer& c){
  //printf("%f\n", c.solution[33]);
}
