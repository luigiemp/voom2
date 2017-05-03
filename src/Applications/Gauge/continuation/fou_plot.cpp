#include<iostream>
#include<math.h>
#include<string.h>
#include<sstream>
#include<gsl/gsl_sf.h>
#include"plstream.h"
#include"plplot.h"
#define PI (3.141592653589793)
int NOE=300;

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

double Hermite(double t, int index, int derv, double* nodes, int element){
  /* xi : lies in [0,1]
     index: 0, 1, 2 or 3 
            0,1 (u_0, u_0') and 2,3 (u_1, u_1')
     derv: derivatives 0,1,2 or 3
  */
  double h = nodes[element]-nodes[element-1];
  double x1 = nodes[element-1];
  double x2 = nodes[element];
  double xi = (t-x1)/(x2-x1);
  if(derv==0){
    if (index == 0)
      return (2*xi*xi*xi - 3*xi*xi +1);
    else if (index==1)
      return (xi*xi*xi - 2*xi*xi + xi)*h;
    else if (index==2)
      return (-2*xi*xi*xi + 3*xi*xi);
    else if (index==3)
      return (xi*xi*xi - xi*xi)*h;
  }
  else if(derv==1){
    if (index == 0)
      return 6*(xi*xi - xi)*1/h;
    else if (index==1)
      return (3*xi*xi - 4*xi +1)*1/h*h;
    else if (index==2)
      return 6*(-xi*xi + xi)*1/h;
    else if (index==3)
      return (3*xi*xi - 2*xi)*1/h*h;
  }
  else if(derv==2){
    if (index == 0)
      return (12*xi - 6)*pow(1/h,2);
    else if (index==1)
      return (6*xi - 4)*pow(1/h,2)*h;
    else if (index==2)
      return (-12*xi + 6)*pow(1/h,2);
    else if (index==3)
      return (6*xi - 2)*pow(1/h,2)*h;
  }
  else if(derv==3){
    if (index == 0)
      return 12.*pow(1/h,3);
    else if (index==1)
      return 6.*pow(1/h,3)*h;
    else if (index==2)
      return -12.*pow(1/h,3);
    else if (index==3)
      return 6.*pow(1/h,3)*h;
  }
}
double Lagrange(double t, int index, int derv, double* nodes, int element){
  /* xi: lies in [0,1]
     index: is 0 or 1
     derv: is 0 or 1
  */
  double h = nodes[element]-nodes[element-1];
  double xi = (t-nodes[element-1])/h;
  if (derv==0){
    if (index==1)
      return xi;
    else if(index==0)
      return (1-xi);
  }
  else if(derv==1){
    if (index==1)
      return 1.*(1/h);
    else if(index==0)
      return -1.*(1/h);

  }

}
int assign_data(double* data_in, int e, char var, double* data_out, double* nodes){
  /* data_out[0] = u0
     data_out[1] = u0_x
     data_out[2] = u1
     data_out[3] = u1_x
  */
  int v_start = 2*(NOE);
  int phi_start = 4*(NOE);
  int phi_end = 6*NOE;
  switch (var){
  case 'u':
    if (e == 1){
      /* Left Boundary */
      data_out[0] = data_in[0];
      data_out[1] = 0; //u0' = 0
      data_out[2] = data_in[1]; //u1
      data_out[3] = data_in[2]; //u1'

    }
    else if (e == NOE){
      /* Right Boundary*/
      data_out[0]=data_in[v_start-3]; //(uN) [....| u(N-1) u(N-1)' uN uN' || v ....]
      data_out[1]=data_in[v_start-2]; // uN'
      data_out[2]=data_in[v_start-1];
      data_out[3]=0; //uN+1'= 0
    }
    else{
      data_out[0]=data_in[(e-2)*2 +1];
      data_out[1]=data_in[(e-2)*2 +2];
      data_out[2]=data_in[(e-2)*2 +3];
      data_out[3]=data_in[(e-2)*2 +4];
    }
    break;
  case 'v':
    if (e == 1){
      /* Left Boundary */
      data_out[0] = 0; //v0=0
      data_out[1] = data_in[v_start];
      data_out[2] = data_in[v_start+1]; //v1
      data_out[3] = data_in[v_start + 2]; //v1'
      
    }
    else if (e == NOE){
      /* Right Boundary*/
      data_out[0]= data_in[phi_start-3]; //vN
      data_out[1]= data_in[phi_start-2];//vN'
      data_out[2]= 0; //vN+1=0
      data_out[3]= data_in[phi_start-1];

    }
    else{
      data_out[0]=data_in[v_start + (e-2)*2+1];
      data_out[1]=data_in[v_start + (e-2)*2 +2];
      data_out[2]=data_in[v_start + (e-2)*2 +3];
      data_out[3]=data_in[v_start + (e-2)*2 +4];
    }
    break;
  case 'c':
    if (e == 1){
      /* Left Boundary */
      data_out[0] = data_in[phi_start + 0];
      data_out[1] = 0; //u0' = 0
      data_out[2] = data_in[phi_start + 1]; //u1
      data_out[3] = data_in[phi_start + 2]; //u1'

    }
    else if (e == NOE){
      /* Right Boundary*/
      data_out[0]=data_in[6*NOE-3]; //(uN) [....| u(N-1) u(N-1)' uN uN' || v ....]
      data_out[1]=data_in[6*NOE-2]; // uN'
      data_out[2]= data_in[6*NOE-1];
      data_out[3]=0; //uN+1'= 0
    }
    else{
      data_out[0]=data_in[phi_start + (e-2)*2 +1];
      data_out[1]=data_in[phi_start + (e-2)*2 +2];
      data_out[2]=data_in[phi_start + (e-2)*2 +3];
      data_out[3]=data_in[phi_start + (e-2)*2 +4];
    }
    break;
  }
  return 0;
}
using namespace std;
class Plot{
public:
  plstream pl;
public:
  void plot_conc(double* data, double* nodes);
  void plot3d_vesicle(double* data, double* nodes);
  void plot_bif_dia(FILE* in_file, int para, int data_pt, int no_of_datapts, double* nodes);
  void plot_lag_mult(FILE* in_file, int para, int data_pt, int no_of_datapts,int lag_pos, double* nodes);
  void plot_r(double* data, double* nodes);
		
};
void file_seek(FILE* in_file, double* data, int data_pt);

void Plot::plot_conc(double* data, double* nodes){
  int N=10*NOE;
  PLFLT* s=new PLFLT[N+1]; //current arc length
  PLFLT* Conc=new PLFLT[N+1]; //assuming r_len=phi_len
  PLFLT max_x, min_x, max_y,min_y;
  PLFLT r,r_theta;
	
  int j=0;
  s[0]=0;
  int start=1; int end=NOE;
  cout<<N<<endl;
  for(int element=start; element<= end; element++){
    for(int i=1; i<=10; i++){
      double x=(i)*1./11;
      double t = nodes[element-1] + (nodes[element]-nodes[element-1])*x;  
      double h =nodes[element]-nodes[element-1];

      double Phi[4];
    //int element = (int) ceil(t*NOE/PI);
    //int element = (int) ceil(asin(sqrt(t/PI))*2*NOE/PI);
    
    assign_data(data, element, 'c', Phi,nodes);
    
    Conc[j] =data[6*NOE+9]+ Phi[0]*Hermite(t,0,0,nodes,element)+Phi[1]*Hermite(t,1,0,nodes,element)+Phi[2]*Hermite(t,2,0,nodes,element)+Phi[3]*Hermite(t,3,0,nodes,element);
    s[j] = t;
    
    if(element==start) { min_x=s[0]; max_x=s[0]; max_y=min_y=Conc[0]; j++;}
    else{
     if(s[j]>=max_x) max_x=s[j];
     if(s[j]<min_x) min_x=s[j];
     if(Conc[j]>=max_y) max_y=Conc[j];
     if(Conc[j]<min_y) min_y=Conc[j];
     j++;
     if(j<=N) s[j]=s[j-1];
    }
    }}

  cout<<min_x<<" "<<max_x<<" "<<min_y<<" "<<max_y<<endl;
  //plstream pl;
  //pl.start("xwin",1,1);
  //pl.adv(1);

  pl.vsta();
  pl.wind(0,PI,-1.5,1.5);
  pl.box( "bcnst", 0, 2, "bcnstv", 0, 2 );
	
  pl.wid(2); //thick lines
  pl.col0(15);
  //pl.line(N,s,Conc);
  
  // change color for each element
  for(int j=0; j<NOE; j++){
  pl.col0(1+(j%2)*(14));
  pl.line(N/NOE,s+j*N/NOE,Conc+j*N/NOE);
  }
  //pl.col0(15);
  //pl.poin(N,s,Phi,1);
  //pl.col0(2); //yellow
  pl.wid(1); //normal line - reset
  pl.lab("s","phi","");
  pl.col0(1); //reset red

}
void Plot::plot_r(double* arg, double* nodes){
  int N=10*NOE;
  
  PLFLT* s=new PLFLT[N+1]; //current arc length
  PLFLT* R=new PLFLT[N+1]; //assuming r_len=phi_len
  PLFLT max_x, min_x, max_y,min_y;
  int j=0;
  double temp;
  int start=2; int end=NOE;
  int plot_option=0;
  for(int element=start; element<= end; element++){
    for(int i=0; i<=9; i++){
      double x;
      if(element==1 && i==0)
	x=0.001;
      else if(element==NOE && i==9)
	x=0.999;
      else
	x=(i)*1./9;
    double t = nodes[element-1] + (nodes[element]-nodes[element-1])*x;  
    double h =nodes[element]-nodes[element-1];
  double p = 0;
  //element = (int) ceil(t*NOE/PI);
  //element = (int) ceil(sqrt(t*NOE*NOE/PI));
  //element = (int) ceil(asin(t/PI)*2*NOE/PI);
  //element = (int) ceil(asin(sqrt(t/PI))*2*NOE/PI);
  
    int v_start = 2*(NOE); 
  int phi_start = 4*(NOE);
  int Total_len = 6*NOE;
  
  double ct=cos(t);
  double st=sin(t);
  double cp=cos(p);
  double sp=sin(p);

  double gamma=arg[Total_len];
  double Lag_phi=arg[Total_len+1];
  double Lag_z=arg[Total_len+2];
  double Lag_gz = arg[Total_len+3];
  double l=arg[Total_len+4];
  double k1=arg[Total_len+5];
	
  double k2=k1*(1+arg[Total_len+6]/100.0);    //<<<------------------
  double beta=arg[Total_len+7];
  double pressure=arg[Total_len+8];
  double mu=arg[Total_len+9];
  // ---------------------Everything that needs plm--------------------------
  double U[4];
  double V[4];
  double Phi[2];
  
  assign_data(arg, element, 'c', Phi, nodes);
  assign_data(arg, element, 'u', U, nodes);
  assign_data(arg, element, 'v', V, nodes);
  
  double u = U[0]*Hermite(t,0,0,nodes,element)+U[1]*Hermite(t,1,0,nodes,element)+U[2]*Hermite(t,2,0,nodes,element)+U[3]*Hermite(t,3,0,nodes,element);
  double u_t = U[0]*Hermite(t,0,1,nodes,element)+U[1]*Hermite(t,1,1,nodes,element)+U[2]*Hermite(t,2,1,nodes,element)+U[3]*Hermite(t,3,1,nodes,element);
  double u_tt = U[0]*Hermite(t,0,2,nodes,element)+U[1]*Hermite(t,1,2,nodes,element)+U[2]*Hermite(t,2,2,nodes,element)+U[3]*Hermite(t,3,2,nodes,element);
  double u_ttt = U[0]*Hermite(t,0,3,nodes,element)+U[1]*Hermite(t,1,3,nodes,element)+U[2]*Hermite(t,2,3,nodes,element)+U[3]*Hermite(t,3,3,nodes,element);
  
  double v = V[0]*Hermite(t,0,0,nodes,element)+V[1]*Hermite(t,1,0,nodes,element)+V[2]*Hermite(t,2,0,nodes,element)+V[3]*Hermite(t,3,0,nodes,element);
  double v_t = V[0]*Hermite(t,0,1,nodes,element)+V[1]*Hermite(t,1,1,nodes,element)+V[2]*Hermite(t,2,1,nodes,element)+V[3]*Hermite(t,3,1,nodes,element);
  double v_tt = V[0]*Hermite(t,0,2,nodes,element)+V[1]*Hermite(t,1,2,nodes,element)+V[2]*Hermite(t,2,2,nodes,element)+V[3]*Hermite(t,3,2,nodes,element);
  double v_ttt = V[0]*Hermite(t,0,3,nodes,element)+V[1]*Hermite(t,1,3,nodes,element)+V[2]*Hermite(t,2,3,nodes,element)+V[3]*Hermite(t,3,3,nodes,element);
  
  double phi =Phi[0]*Hermite(t,0,0,nodes,element)+Phi[1]*Hermite(t,1,0,nodes,element)+Phi[2]*Hermite(t,2,0,nodes,element)+Phi[3]*Hermite(t,3,0,nodes,element);
  double phi_t = Phi[0]*Hermite(t,0,1,nodes,element)+Phi[1]*Hermite(t,1,1,nodes,element)+Phi[2]*Hermite(t,2,1,nodes,element)+Phi[3]*Hermite(t,3,1,nodes,element);
  
  //v= 100000*v; v_t=100000*v_t; v_tt=100000*v_tt;v_ttt=100000*v_ttt;

  double u_p=0.;
  double u_tp=0.;
  double u_pp=0.;
  double u_ttp=0; double u_tpp=0; double u_ppp=0;

  double v_p=0.;
  double v_tp=0.;
  double v_pp=0.;
  double v_ttp=0; double v_tpp=0; double v_ppp=0;
  
  double w=0;
  double w_t= 0;
  double w_p=0.;
  double w_tt= 0;
  double w_tp=0.;
  double w_pp=0.;
  double w_ttt= 0; double w_ttp=0; double w_tpp=0; double w_ppp=0;
  
  double phi_p=0.;
	
  phi = phi+mu;
  u=u+1.; // In the following (u,v,w) represents current configuration
  

  // Metric Tensor
  double At = (u_t-v); double Bt = (u+v_t); double Ct = w_t;
  double At_T = (u_tt-v_t); double Bt_T = (u_t+v_tt); double Ct_T = w_tt;
  double At_P = (u_tp-v_p); double Bt_P = (u_p+v_tp); double Ct_P = w_tp;
	
  double At_TT = (u_ttt-v_tt); double Bt_TT = (u_tt+v_ttt); double Ct_TT = w_ttt;
  double At_TP = (u_ttp-v_tp); double Bt_TP = (u_tp+v_ttp); double Ct_TP = w_ttp;
  double At_PP = (u_tpp-v_pp); double Bt_PP = (u_pp+v_tpp); double Ct_PP = w_tpp;
	
	
  double Ap = (u_p-w*st); double Bp = (v_p-w*ct); double Cp = (u*st+v*ct+w_p);
  double Ap_T = (u_tp-w_t*st-w*ct); double Bp_T = (v_tp-w_t*ct+w*st); double Cp_T = (u_t*st+u*ct+v_t*ct-v*st+w_tp);
  double Ap_P = (u_pp-w_p*st); double Bp_P = (v_pp-w_p*ct); double Cp_P = (u_p*st+v_p*ct+w_pp);
	
  double Ap_TT = (u_ttp-w_tt*st-2*w_t*ct+w*st); 
  double Bp_TT = (v_ttp-w_tt*ct+2*w_t*st+w*ct); 
  double Cp_TT = (u_tt*st+2*u_t*ct-u*st+v_tt*ct-2*v_t*st-v*ct+w_ttp);
  double Ap_TP = (u_tpp-w_tp*st-w_p*st); double Bp_TP = (v_tpp-w_tp*ct+w_p*st); double Cp_TP = (u_tp*st+u_p*ct+v_tp*ct-v_p*st+w_tpp);
  double Ap_PP = (u_ppp-w_pp*st); double Bp_PP = (v_ppp-w_pp*ct); double Cp_PP = (u_pp*st+v_pp*ct+w_ppp);
	
  double g_tt = At*At+Bt*Bt+Ct*Ct;
  double g_tp = At*Ap+Bt*Bp+Ct*Cp;
  double g_pp = Ap*Ap+Bp*Bp+Cp*Cp;
  double g = g_tt*g_pp-g_tp*g_tp;
  double Sg = pow(g,.5); double SSg = pow(g_tt,.5)*(u*st+v*ct);
  double gtt = g_pp/g; double gtp = -g_tp/g; double gpp = g_tt/g;
	
  double g_tt_T = 2*(At*At_T+Bt*Bt_T+Ct*Ct_T);
  double g_tp_T = At_T*Ap+At*Ap_T + Bt_T*Bp+Bt*Bp_T + Ct_T*Cp+Ct*Cp_T;
  double g_pp_T = 2*(Ap*Ap_T+Bp*Bp_T+Cp*Cp_T);
  //double g_tt_T = (2*g_pp*st*ct -st*st*g_pp_T)/(g_pp*g_pp);
		
  double g_tt_P = 2*(At*At_P+Bt*Bt_P+Ct*Ct_P);
  double g_tp_P = At_P*Ap+At*Ap_P + Bt_P*Bp+Bt*Bp_P + Ct_P*Cp+Ct*Cp_P;
  double g_pp_P = 2*(Ap*Ap_P+Bp*Bp_P+Cp*Cp_P);
  //double g_tt_P = 0;
	
  double g_T = g_tt_T*g_pp+g_tt*g_pp_T -2*g_tp*g_tp_T;
  double g_P = g_tt_P*g_pp+g_tt*g_pp_P -2*g_tp*g_tp_P;
  double Sg_T = (g_T)/(2*Sg);
  //double J_T =(st*Sg_T-Sg*ct)/(st*st);
	
  double gtt_T = g_pp_T/g - g_pp*g_T/(g*g);
  double gtt_P = g_pp_P/g - g_pp*g_P/(g*g);
  double gtp_T = -(g_tp_T/g - g_tp*g_T/(g*g));
  double gtp_P = -(g_tp_P/g - g_tp*g_P/(g*g));
  double gpp_T = g_tt_T/g - g_tt*g_T/(g*g);
  double gpp_P = g_tt_P/g - g_tt*g_P/(g*g);
	
  // Second Derivative of Position Vector
  double Rr_tt = At_T-Bt; 
  double Tr_tt = Bt_T+At; 
  double Pr_tt = Ct_T;
	
	
  double Rr_tp = At_P-st*Ct; 
  double Tr_tp = Bt_P-ct*Ct; 
  double Pr_tp = At*st+Bt*ct+Ct_P;
	
  double Rr_pp = Ap_P-st*Cp;
  double Tr_pp = Bp_P-ct*Cp;
  double Pr_pp = Ap*st+Bp*ct+Cp_P;
  
  double Rr_tt_T = At_TT-Bt_T; 
  double Tr_tt_T = Bt_TT+At_T; 
  double Pr_tt_T = Ct_TT;
  double Rr_tt_P = At_TP-Bt_P; 
  double Tr_tt_P = Bt_TP+At_P; 
  double Pr_tt_P = Ct_TP;
	
  double Rr_tp_T = At_TP-st*Ct_T-ct*Ct; 
  double Tr_tp_T = Bt_TP-ct*Ct_T+st*Ct; 
  double Pr_tp_T = At_T*st+At*ct+Bt_T*ct-Bt*st+Ct_TP;
  double Rr_tp_P = At_PP-st*Ct_P; 
  double Tr_tp_P = Bt_PP-ct*Ct_P; 
  double Pr_tp_P = At_P*st+Bt_P*ct+Ct_PP;
	
  double Rr_pp_T = Ap_TP-st*Cp_T-ct*Cp;
  double Tr_pp_T = Bp_TP-ct*Cp_T+st*Cp;
  double Pr_pp_T = Ap_T*st+Ap*ct+Bp_T*ct-Bp*st+Cp_TP;
  double Rr_pp_P = Ap_PP-st*Cp_P;
  double Tr_pp_P = Bp_PP-ct*Cp_P;
  double Pr_pp_P = Ap_P*st+Bp_P*ct+Cp_PP;
  
  // Normal Vector 
  double Rn = (Bt*Cp-Bp*Ct)/Sg;
  double Tn = (Ct*Ap-Cp*At)/Sg;
  double Pn = (At*Bp-Ap*Bt)/Sg;
	
  
  double Rn_T = (Bt_T*Cp+Bt*Cp_T-Bp_T*Ct-Bp*Ct_T)/Sg-0.5*Rn*g_T/g;
  double Tn_T = (Ct_T*Ap+Ct*Ap_T-Cp_T*At-Cp*At_T)/Sg-0.5*Tn*g_T/g;
  double Pn_T = (At_T*Bp+At*Bp_T-Ap_T*Bt-Ap*Bt_T)/Sg-0.5*Pn*g_T/g;
	
  double Rn_P = (Bt_P*Cp+Bt*Cp_P-Bp_P*Ct-Bp*Ct_P)/Sg-0.5*Rn*g_P/g;
  double Tn_P = (Ct_P*Ap+Ct*Ap_P-Cp_P*At-Cp*At_P)/Sg-0.5*Tn*g_P/g;
  double Pn_P = (At_P*Bp+At*Bp_P-Ap_P*Bt-Ap*Bt_P)/Sg-0.5*Pn*g_P/g;
  
  // Second Fundamental Form
  double b_tt = Rn*Rr_tt+Tn*Tr_tt+Pn*Pr_tt;
  double b_tp = Rn*Rr_tp+Tn*Tr_tp+Pn*Pr_tp;
  double b_pp = Rn*Rr_pp+Tn*Tr_pp+Pn*Pr_pp;
	
  //Curvature Tensor
  double bt_t = gtt*b_tt+gtp*b_tp;
  double bt_p = gtt*b_tp+gtp*b_pp;
  double bp_t = gtp*b_tt+gpp*b_tp;
  double bp_p = gtp*b_tp+gpp*b_pp;
	
  double btt = gtt*bt_t+gtp*bt_p; 
  double btp = gtt*bp_t+gtp*bp_p; 
  double bpp = gtp*bp_t+gpp*bp_p;
  //-----------------DO NOT DELETE------------------------------
  
  double b_tt_T = Rn_T*Rr_tt+Rn*Rr_tt_T + Tn_T*Tr_tt+Tn*Tr_tt_T + Pn_T*Pr_tt+Pn*Pr_tt_T;
  double b_tt_P = Rn_P*Rr_tt+Rn*Rr_tt_P + Tn_P*Tr_tt+Tn*Tr_tt_P + Pn_P*Pr_tt+Pn*Pr_tt_P;
  double b_tp_T = Rn_T*Rr_tp+Rn*Rr_tp_T + Tn_T*Tr_tp+Tn*Tr_tp_T + Pn_T*Pr_tp+Pn*Pr_tp_T;
  double b_tp_P = Rn_P*Rr_tp+Rn*Rr_tp_P + Tn_P*Tr_tp+Tn*Tr_tp_P + Pn_P*Pr_tp+Pn*Pr_tp_P;
  double b_pp_T = Rn_T*Rr_pp+Rn*Rr_pp_T + Tn_T*Tr_pp+Tn*Tr_pp_T + Pn_T*Pr_pp+Pn*Pr_pp_T;
  double b_pp_P = Rn_P*Rr_pp+Rn*Rr_pp_P + Tn_P*Tr_pp+Tn*Tr_pp_P + Pn_P*Pr_pp+Pn*Pr_pp_P;
	
  double bt_t_T = gtt_T*b_tt+gtt*b_tt_T + gtp_T*b_tp+gtp*b_tp_T;
  double bt_t_P = gtt_P*b_tt+gtt*b_tt_P + gtp_P*b_tp+gtp*b_tp_P;
  double bt_p_T = gtt_T*b_tp+gtt*b_tp_T + gtp_T*b_pp+gtp*b_pp_T;
  double bt_p_P = gtt_P*b_tp+gtt*b_tp_P + gtp_P*b_pp+gtp*b_pp_P;
  double bp_t_T = gtp_T*b_tt+gtp*b_tt_T + gpp_T*b_tp+gpp*b_tp_T;
  double bp_t_P = gtp_P*b_tt+gtp*b_tt_P + gpp_P*b_tp+gpp*b_tp_P;
  double bp_p_T = gtp_T*b_tp+gtp*b_tp_T + gpp_T*b_pp+gpp*b_pp_T;
  double bp_p_P = gtp_P*b_tp+gtp*b_tp_P + gpp_P*b_pp+gpp*b_pp_P;
  
  //----------------------------------------------------------------
  //Curvatures
  double H = 0.5*(bt_t+bp_p); double K = (bt_t*bp_p - bt_p*bp_t);
  double H_T = 0.5*(bt_t_T+bp_p_T); double H_P = 0.5*(bt_t_P+bp_p_P);
  double HT = (gtt*H_T+gtp*H_P); double HP = gtp*H_T+gpp*H_P;

  double Lap_herm_0 = gtt*Hermite(t,0,2,nodes,element)+(gtt_T+gtt*g_T/(2*g))*Hermite(t,0,1,nodes,element);
  double Lap_herm_1 = gtt*Hermite(t,1,2,nodes,element)+(gtt_T+gtt*g_T/(2*g))*Hermite(t,1,1,nodes,element);
  double Lap_herm_2 = gtt*Hermite(t,2,2,nodes,element)+(gtt_T+gtt*g_T/(2*g))*Hermite(t,2,1,nodes,element);
  double Lap_herm_3 = gtt*Hermite(t,3,2,nodes,element)+(gtt_T+gtt*g_T/(2*g))*Hermite(t,3,1,nodes,element);

  double Norm_eq_part =0;
  double Tgt_eq =0;
  double m1=-1.;
  double m2=1.;
  double epsilon=-W2prime(mu,1.0,m1,m2)/(l*(l+1));
  //cout<< "epsilon="<<epsilon<<endl;
  Norm_eq_part = ( 2*H*k1*(H*H-K)-2*H*gamma -pressure+Lag_z*(ct*(Rn-2*H*u)-st*(Tn-2*H*v)));
  double kappa=stiffness(phi,m1,m2,k1,k2); //<<<<------------------------
  cout <<kappa<< " "<<t<<endl;
//Current arc length 
  if (j == 0)
    s[j] = sqrt(g_tt)*PI/N;
  else
    s[j] = s[j-1]+sqrt(g_tt)*PI/N;
  
  s[j]= t;element;
  
  temp  = (g_tt*st*st)-g_pp;nodes[element];1*(k1*H*Lap_herm_1 + 0*Norm_eq_part*Hermite(t,1,0,nodes,element))*Sg*(element%2==0)+
         1*(k1*H*Lap_herm_3 + 0*Norm_eq_part*Hermite(t,3,0,nodes,element))*Sg*(element%2==1);

  ( epsilon*gtt*phi_t*Lagrange(t,0,1,nodes,element)+0*(Wprime(phi,beta,m1,m2)+Lag_phi)*Lagrange(t,0,0,nodes,element) )*Sg*(element%2==0)+
         ( epsilon*gtt*phi_t*Lagrange(t,1,1,nodes,element)+0*(Wprime(phi,beta,m1,m2)+Lag_phi)*Lagrange(t,1,0,nodes,element) )*Sg*(element%2==1);


;1*(k1*H*Lap_herm_3 + Norm_eq_part*Hermite(t,3,0,nodes,element))*Sg; 
gtt*Hermite(t,3,2,nodes,element)+(gtt_T+gtt*g_T/(2*g))*Hermite(t,3,1,nodes,element); 
(g_tt*st*st)/g_pp-1;(1-sqrt(1-g_tt*st*st))/(2*sqrt(g_tt)*sin(s[j]/2)*sin(s[j]/2))   ;//theta*PI/180+atan2(v,u); (u*st+v*ct);
  v_t; ct*u_tt+st*v_tt-v*ct-u*st-2*st*v_t+2*ct*u_t;

  R[j] = temp;
  //if (element==1) printf("%1.15le\n",H);
  if(temp==(g_tt*st*st)-g_pp) plot_option=1;
  
  if(element==start) { min_x=s[0]; max_x=s[0]; max_y=min_y=R[0]; j++;}
  else{
  if(s[j]>=max_x) max_x=s[j];
  if(s[j]<min_x) min_x=s[j];
  if(R[j]>=max_y) max_y=R[j];
  if(R[j]<min_y) min_y=R[j];
  j++;
  if(j<=N) s[j]=s[j-1];
  }
  
    }
  }
  cout<<min_x<<" "<<max_x<<" "<<min_y<<" "<<max_y<<"...."<<1/(nodes[1]-nodes[0])<<endl;
  //plstream pl;
  //pl.start("xwin",1,1);
  //pl.adv(1);
	
  pl.vsta();
  //min_x=2.00; 
  //max_x=3.3;
  //min_y=-2.10;
  //max_y=-2.08;
  if(plot_option==1)
    pl.wind(min_x,max_x, 0-1e-5,0+1e-5);
  else
    pl.wind(min_x,max_x, min_y-1e-8,max_y+1e-8);
  pl.box( "bcnst", 0, 2, "bcnstv", 0, 2 );
	
  pl.wid(2); //thick lines
  //------ uniform line
  pl.col0(15);
  //pl.line(N,s,R);
  // change color for each element
  
  //pl.col0(3);
  //pl.line(N/NOE,s+0*N/NOE,R+0*N/NOE);
  for(int j=0; j<NOE; j++){
  pl.col0(1+(j%2)*(14));
  pl.line(N/NOE,s+j*N/NOE,R+j*N/NOE);
  }
  
  pl.wid(1); //normal line - reset
  pl.lab("s","H","");
  pl.col0(1); //reset red

}

void Plot::plot3d_vesicle(double* data, double* nodes){
  int N=500;
  PLFLT* R=new PLFLT[N+1];
  PLFLT* Conc=new PLFLT[N+1]; //assuming r_len=phi_len
  PLFLT* Z=new PLFLT[N+1];
  PLFLT max_x, min_x, max_y,min_y;
  int j=0;
  for(int i=0; i<=N; i++){
    PLFLT t=((179.9999/N)*i+.000001)*PI/180;
    PLFLT ct,st;
    ct=cos(t);
    st=sin(t);
    R[j]=st;
    Z[j]=ct;
    Conc[j]=0;
    PLFLT u=0.0; PLFLT v=0;
    double U[4];
    double V[4];
    double Phi[2];
    //int element = (int) ceil(t*NOE/PI);
    int element = (int) ceil(asin(sqrt(t/PI))*2*NOE/PI);
    //double x=(t*NOE/PI +1 -element);
    assign_data(data, element, 'c', Phi,nodes);
    assign_data(data, element, 'u', U,nodes);
    assign_data(data, element, 'v', V,nodes);
    u = U[0]*Hermite(t,0,0,nodes,element)+U[1]*Hermite(t,1,0,nodes,element)+U[2]*Hermite(t,2,0,nodes,element)+U[3]*Hermite(t,3,0,nodes,element);
    v = V[0]*Hermite(t,0,0,nodes,element)+V[1]*Hermite(t,1,0,nodes,element)+V[2]*Hermite(t,2,0,nodes,element)+V[3]*Hermite(t,3,0,nodes,element);
    Conc[j] = data[6*NOE+9]+ Phi[0]*Hermite(t,0,0,nodes,element)+Phi[1]*Hermite(t,1,0,nodes,element)+Phi[2]*Hermite(t,2,0,nodes,element)+Phi[3]*Hermite(t,3,0,nodes,element);
    R[j]+=u*st+v*ct;
    Z[j]+=u*ct-v*st;
    if(i==0) { min_x=R[0]; max_x=R[0]; max_y=min_y=Z[0]; }
    if(R[j]>=max_x) max_x=R[j];
    if(R[j]<min_x) min_x=R[j];
    if(Z[j]>=max_y) max_y=Z[j];
    if(Z[j]<min_y) min_y=Z[j];
    j++;
		
  }
  //plstream pl;
  //pl.start("xwin",1,1);
  //pl.adv(1);
  pl.col0(0);  //Black for the border
  pl.vasp((max_y-min_y)/(2*max_x));
  pl.wind(-max_x,max_x,min_y,max_y);
	
  double R_prev=R[0];
  double Z_prev=Z[0];
  double X[2];
  double Y[2];
  pl.wid(2);
  for(int i=1;i <=N; i++){
    X[0]=R_prev;
    Y[0]=Z_prev;
    X[1]=R[i];
    Y[1]=Z[i];
    if(Conc[i-1]>=0.5){
      pl.col0(1);
      pl.line(2,X,Y);
      //pl.poin(2,X,Y,1);
      X[0]=-R_prev;
      X[1]=-R[i];
      pl.line(2,X,Y);
    }
    else if(Conc[i-1]<-0.5){
      pl.col0(12);
      pl.line(2,X,Y);
      //pl.poin(2,X,Y,1);
      X[0]=-R_prev;
      X[1]=-R[i];
      pl.line(2,X,Y);
    }
    else {
      pl.col0(7);
      pl.line(2,X,Y);
      //pl.poin(2,X,Y,1);
      X[0]=-R_prev;
      X[1]=-R[i];
      pl.line(2,X,Y);
    }
    R_prev=R[i];
    Z_prev=Z[i];
    pl.col0(1); //red
		
  }//pl.line(mesh_size,R,Z);
  double origin_x[1]={0};
  double origin_y[1]={0};
  pl.poin(1,origin_x,origin_y,1);
  pl.wid(1);
 
}

void Plot::plot_bif_dia(FILE* in_file, int para, int data_pt, int no_of_datapts, double* nodes){

  double* data = new double[6*(NOE)+ 10];
  PLFLT* nrm= new PLFLT[no_of_datapts];
  PLFLT* par = new PLFLT[no_of_datapts]; 
  PLFLT nrm_datapt, par_datapt;
  PLFLT max_x, min_x, max_y,min_y;
  double temp=0;
  rewind(in_file);
  for(int i=0; i<no_of_datapts; i++){
    for(int j=0;j<(6*(NOE)+10);j++){
      fscanf(in_file,"%le", &temp);
      data[j]=temp;
    }
    for(int k=0;k< 6*(NOE) + 4; k++) nrm[i] += data[k]*data[k]; //squared nrm
    nrm[i] = sqrt(nrm[i]); //L2 nrm
    par[i] = data[ 6*(NOE) + 4 + para];
    if(i==0) {max_x=min_x=par[0]; max_y=min_y=nrm[0];}
    if(par[i]>=max_x) max_x=par[i];
    if(par[i]<min_x) min_x=par[i];
    if(nrm[i]>=max_y) max_y=nrm[i];
    if(nrm[i]<min_y) min_y=nrm[i];
		
    if(i==data_pt-1) {nrm_datapt=nrm[i]; par_datapt=par[i];}
  }
  //plstream pl;
  //pl.start("xwin", 1,1);
  //pl.adv(1);
  cout<<"Par = "<<par_datapt<<endl;

  pl.vsta();
  pl.wind(min_x,max_x,min_y,max_y);
  pl.box( "bcnst", 0, 2, "bcnstv", 0, 2 );
	
  pl.wid(1); //thicker lines
  pl.col0(7);
  pl.line(no_of_datapts,par,nrm);
  //pl.poin(no_of_datapts,par,nrm,'*');
  pl.col0(15);
  pl.poin(1,&par_datapt,&nrm_datapt,'*');
  pl.wid(1); //normal line - reset
  pl.col0(1);
  stringstream ss;
  string title;
  ss<<par_datapt;
  title=ss.str();
  pl.lab("para","",title.c_str());
  pl.col0(1); //reset red
  //pl.vasp(1);

}

void Plot::plot_lag_mult(FILE* in_file, int para, int data_pt, int no_of_datapts,int lag_pos, double* nodes){
  /*
    int u_len = r_len[0];
    int v_len =r_len[1];
    int w_len =r_len[2];
    double* data = new double[r_len + phi_len + 9];
    PLFLT* lag_mult= new PLFLT[no_of_datapts];
    PLFLT* par = new PLFLT[no_of_datapts]; 
    PLFLT lag_mult_datapt, par_datapt;
    PLFLT max_x, min_x, max_y,min_y;
	
    for(int i=0; i<no_of_datapts; i++){
    file_seek(in_file,data, i+1, r_len, phi_len);
    lag_mult[i]=data[r_len + phi_len+ lag_pos];
    par[i] = data[ r_len + phi_len+ 3 + para];
    //cout<<i<<endl;
    if(i==0) {max_x=min_x=par[0]; max_y=min_y=lag_mult[0];}
    if(par[i]>=max_x) max_x=par[i];
    if(par[i]<min_x) min_x=par[i];
    if(lag_mult[i]>=max_y) max_y=lag_mult[i];
    if(lag_mult[i]<min_y) min_y=lag_mult[i];
		
    if(i==data_pt) {lag_mult_datapt=lag_mult[i]; par_datapt=par[i];}
    }
    //plstream pl;
    //pl.start("xwin", 1,1);
    //pl.adv(1);
    pl.vsta();
    pl.wind(min_x,max_x,min_y,max_y);
    pl.box( "bcnst", 0, 2, "bcnstv", 0, 2 );
	
    pl.wid(2); //thicker lines
    pl.col0(7);
    pl.line(no_of_datapts,par,lag_mult);
    //pl.poin(no_of_datapts,par,lag_mult,'*');
    pl.col0(9);
    pl.poin(1,&par_datapt,&lag_mult_datapt,'*');
    pl.wid(1); //normal line - reset
    pl.col0(2);
    pl.lab("para","","");
    pl.col0(1); //reset red

  */
}

void file_seek(FILE* in_file, double* data, int data_pt){
  
  rewind(in_file);
  double temp;
  for(int i=0;i<(data_pt-1)*(6*(NOE)+10);i++){
    fscanf(in_file,"%le", &temp);
  }
  for(int i=0;i<(6*(NOE)+10);i++){
    fscanf(in_file,"%le", &temp);
    if (i<= 6*(NOE)){
      data[i]=temp;
      //cout<<data[i];
    }
    else
      data[i]=temp;
  }
}
int main(int argc, char** argv){
  /*------------Parsing command line arguments --------------------*/
  char* in_filename, *out_filename, *cont_para;
  int data_pt, para, opt;
  bool flag_vesicle, flag_conc_plot, flag_bif_dia_vesicle, flag_bif_dia,flag_3d,flag_lag_mult;
  int lag_pos=0;
  double* nodes = new double[NOE+1];
  
  
  flag_vesicle=flag_conc_plot=flag_bif_dia_vesicle=flag_bif_dia=flag_3d=flag_lag_mult=false;
  in_filename = out_filename = cont_para = NULL;
  while ( (opt=getopt(argc,argv,"fDobVcBpuvwMh3l"))!=-1){
    switch (opt){
    case 'f': //set input filename
      in_filename=argv[optind];
			
      break;
    case 'D': //sets the data pt
      data_pt=atoi(argv[optind]);
      break;
    case 'o': //Output filename - a ps file
      out_filename=argv[optind];
      break;
    case 'b': //plots bifurcation diagram
      flag_bif_dia=true;
      break;
    case 'V': //plots vesicle
      flag_vesicle=true;
      break;
    case 'c':  // plots conc vs arclength
      flag_conc_plot=true;
      break;
    case 'B': //plots bifurcation dia and vesicle
      flag_bif_dia_vesicle=true;
      break;
    case 'p':
      cont_para=argv[optind];
      break;
    case 'l':
      flag_lag_mult = true;
      lag_pos=atoi(argv[optind]);
      break;
    case '3':
      flag_3d=true;
      break;
    case 'M':
      NOE=atoi(argv[optind]);
      //case 'h': //prints help
      //print_help();
    }
  }

	
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
	
  int no_of_elements;
	
  string str("wc ");
  string str2(in_filename);
  str=str+str2+"> wc_tmp.txt"; //wc "filename" > wc_tmp.dat" 
  system(str.c_str());
  FILE* pfile_tmp=fopen("wc_tmp.txt","r");
	
  fscanf(pfile_tmp,"%d", &no_of_elements);
  fclose(pfile_tmp);
  system("rm wc_tmp.txt");
  //plstream* pl;
  //pl->start("xwin", 1,1);
  FILE* in_file=fopen(in_filename,"r");
	
  Plot p;
	
  for (int n=0; n<= NOE; n++){
    //nodes[n] = n*PI/NOE; // UNIFORM
    //nodes[n] =n*n*PI/(NOE*NOE);
    nodes[n] = PI*pow(sin(n*PI/(2.*NOE)),2);
    //cout << n << " " << nodes[n]*180/PI << endl;
  }
	
  double* data=new double[6*(NOE)+10]; //modes + para
  file_seek(in_file, data,data_pt);	
  //for (int i=0; i<no_of_modes*2+9; i++)
  //	cout<<data[i]<<endl;
  
  // Load colour palettes
  //p.pl.spal0( "cmap0_black_on_white.pal" );
  //p.pl.spal1( "cmap1_blue_red.pal", true );
  //p.pl.spal0( "cmap0_black_on_white.pal" );
  //p.pl.spal1( "cmap1_blue_red.pal", true );
  //p.pl.scmap0n( 3 );

  if(flag_vesicle==true){
    p.pl.start("xwin",1,1);
    p.pl.adv(1);
    p.plot3d_vesicle(data,nodes);
  }
  if(flag_conc_plot == true){
    p.pl.start("xwin",1,1);
    p.pl.adv(1);
    //p.plot_conc(data,r_len,phi_len);
    p.plot_r(data,nodes); 
	
    //p.plot_spectrum(data,no_of_modes);
  }
  if(flag_bif_dia ==true){
    p.pl.start("xwin",1,1);
    p.pl.adv(1);
    p.plot_bif_dia(in_file, para, data_pt, no_of_elements,nodes);
  }
  if(flag_lag_mult ==true){
    p.pl.start("xwin",1,1);
    p.pl.adv(1);
    p.plot_lag_mult(in_file, para, data_pt, no_of_elements, lag_pos, nodes);
  }
  if(flag_bif_dia_vesicle==true){
    p.pl.start("xwin",2,2);
    //Vesicle 
    p.pl.adv(1);
    p.plot3d_vesicle(data, nodes);
		
    //Concentration profile
    p.pl.adv(2);
    p.plot_conc(data, nodes);
    //p.plot_r(data, nodes);
    //p.plot_spectrum(data,r_len, phi_len);
		
    // Lagrange multiplier plot or r plot
    // Last argument is either 0 - gamma, 1 -Lag_phi, 2 - Lag_z
    p.pl.adv(3);
    //p.plot_lag_mult(in_file, para, data_pt, r_len, phi_len, no_of_elements, 0); 
    p.plot_r(data, nodes); 
	
    //Bifurcation diagram
    p.pl.adv(4);
    p.plot_bif_dia(in_file, para, data_pt, no_of_elements, nodes);
    if(out_filename!=NULL)
      {
	p.pl.scolbg(255,255,255);
	//p.pl.scol0(0,0,0,0);
	printf("The current plot will be saved as a color Postscript with filename %s\n", out_filename);
	string str1(out_filename);
	str1=str1+".fig";
	plstream *pls2;
	pls2 = new plstream();       /* create a new one */
	pls2->sfnam( str1.c_str());       /* file name a fig file*/
	pls2->sdev( "xfig" );         /* device type */
	pls2->cpstrm( p.pl, false ); /* copy old stream parameters to new stream */
	pls2->replot();              /* do the save by replaying the plot buffer */
	delete pls2;
	str1="fig2mpdf -e "+str1; //command: fig2mpdf <filename.fig> outputs a eps file
	system(str1.c_str());
	str1="rm "+string(out_filename)+".fig";
	system(str1.c_str());
      }
		
	
	
  }
  cout<<no_of_elements<<endl;
	
}


