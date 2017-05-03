#include<iostream>
#include<gsl/gsl_math.h>
#include<gsl/gsl_sf.h>
#include"uvw_model2_FV.h"
#include<math.h>
using namespace std;
#define PI (3.141592653589793)

int local_firstvariation(double* indvar,double* arg,int element, double* out, double* nodes){

  /* Inputs: indvar - the independent variable, theta (colatitude)
             arg - All the unknowns and para
	     element - 1,2, ..., NOE (element number)
	             - if element = 0, we process constraints
	     out - is the return value. It contains the integrand
	           of the weak formulation (ie, normal_eq* H_j, tangt_eq*H_j, conc_eq*N_i)
		   for j = 0,1,2,3 and i=0,1 - that is, 'out' holds 8 doubles.
   */
  //PLM_COUNT++;
  double x = indvar[0];
  double t = nodes[element-1] + (nodes[element]-nodes[element-1])*x;
  double p = indvar[1];
  double dtdx = nodes[element]-nodes[element-1];
  double h = dtdx;
  int v_start = 2*NOE; 
  int phi_start = 4*NOE;
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
  double Phi[4];
  
  assign_data(arg, element, 'c', Phi, nodes);
  assign_data(arg, element, 'u', U, nodes);
  assign_data(arg, element, 'v', V, nodes);
    
  double Herm0 = Hermite(t,0,0,nodes,element);
  double Herm1 = Hermite(t,1,0,nodes,element);
  double Herm2 = Hermite(t,2,0,nodes,element);
  double Herm3 = Hermite(t,3,0,nodes,element);

  double Herm0_t = Hermite(t,0,1,nodes,element);
  double Herm1_t = Hermite(t,1,1,nodes,element);
  double Herm2_t = Hermite(t,2,1,nodes,element);
  double Herm3_t = Hermite(t,3,1,nodes,element);

  double Herm0_tt = Hermite(t,0,2,nodes,element);
  double Herm1_tt = Hermite(t,1,2,nodes,element);
  double Herm2_tt = Hermite(t,2,2,nodes,element);
  double Herm3_tt = Hermite(t,3,2,nodes,element);

  double Herm0_ttt = Hermite(t,0,3,nodes,element);
  double Herm1_ttt = Hermite(t,1,3,nodes,element);
  double Herm2_ttt = Hermite(t,2,3,nodes,element);
  double Herm3_ttt = Hermite(t,3,3,nodes,element);

  double u = U[0]*Herm0+U[1]*Herm1+U[2]*Herm2+U[3]*Herm3;
  double u_t = U[0]*Herm0_t+U[1]*Herm1_t+U[2]*Herm2_t+U[3]*Herm3_t;
  double u_tt = U[0]*Herm0_tt+U[1]*Herm1_tt+U[2]*Herm2_tt+U[3]*Herm3_tt;
  double u_ttt = U[0]*Herm0_ttt+U[1]*Herm1_ttt+U[2]*Herm2_ttt+U[3]*Herm3_ttt;
  
  double v = V[0]*Herm0+V[1]*Herm1+V[2]*Herm2+V[3]*Herm3;
  double v_t = V[0]*Herm0_t+V[1]*Herm1_t+V[2]*Herm2_t+V[3]*Herm3_t;
  double v_tt = V[0]*Herm0_tt+V[1]*Herm1_tt+V[2]*Herm2_tt+V[3]*Herm3_tt;
  double v_ttt = V[0]*Herm0_ttt+V[1]*Herm1_ttt+V[2]*Herm2_ttt+V[3]*Herm3_ttt;
  
  double phi = Phi[0]*Herm0+Phi[1]*Herm1+Phi[2]*Herm2+Phi[3]*Herm3;
  double phi_t = Phi[0]*Herm0_t+Phi[1]*Herm1_t+Phi[2]*Herm2_t+Phi[3]*Herm3_t;
  double phi_p = 0;
  phi = phi+mu;
  u=u+1.; // In the following (u,v,w) represents current configuration
  
  /*
  double u = U[0]*Hermite(t,0,0,nodes,element)+U[1]*Hermite(t,1,0,nodes,element)+U[2]*Hermite(t,2,0,nodes,element)+U[3]*Hermite(t,3,0,nodes,element);
  double u_t = U[0]*Hermite(t,0,1,nodes,element)+U[1]*Hermite(t,1,1,nodes,element)+U[2]*Hermite(t,2,1,nodes,element)+U[3]*Hermite(t,3,1,nodes,element);
  double u_tt = U[0]*Hermite(t,0,2,nodes,element)+U[1]*Hermite(t,1,2,nodes,element)+U[2]*Hermite(t,2,2,nodes,element)+U[3]*Hermite(t,3,2,nodes,element);
  double u_ttt = U[0]*Hermite(t,0,3,nodes,element)+U[1]*Hermite(t,1,3,nodes,element)+U[2]*Hermite(t,2,3,nodes,element)+U[3]*Hermite(t,3,3,nodes,element);
  
  double v = V[0]*Hermite(t,0,0,nodes,element)+V[1]*Hermite(t,1,0,nodes,element)+V[2]*Hermite(t,2,0,nodes,element)+V[3]*Hermite(t,3,0,nodes,element);
  double v_t = V[0]*Hermite(t,0,1,nodes,element)+V[1]*Hermite(t,1,1,nodes,element)+V[2]*Hermite(t,2,1,nodes,element)+V[3]*Hermite(t,3,1,nodes,element);
  double v_tt = V[0]*Hermite(t,0,2,nodes,element)+V[1]*Hermite(t,1,2,nodes,element)+V[2]*Hermite(t,2,2,nodes,element)+V[3]*Hermite(t,3,2,nodes,element);
  double v_ttt = V[0]*Hermite(t,0,3,nodes,element)+V[1]*Hermite(t,1,3,nodes,element)+V[2]*Hermite(t,2,3,nodes,element)+V[3]*Hermite(t,3,3,nodes,element);
  
  double phi = Phi[0]*Hermite(t,0,0,nodes,element)+Phi[1]*Hermite(t,1,0,nodes,element)+Phi[2]*Hermite(t,2,0,nodes,element)+Phi[3]*Hermite(t,3,0,nodes,element);
  double phi_t = Phi[0]*Hermite(t,0,1,nodes,element)+Phi[1]*Hermite(t,1,1,nodes,element)+Phi[2]*Hermite(t,2,1,nodes,element)+Phi[3]*Hermite(t,3,1,nodes,element);
  */

  /*
  //-------Atha: Slower --------

  double u_p=0.;
  double u_tp=0.;
  double u_pp=0.;
  double u_ttp=0; double u_tpp=0; double u_ppp=0;

  double v_p=0.;
  double v_tp=0.;
  double v_pp=0.;
  double v_ttp=0; double v_tpp=0; double v_ppp=0;
  
  //v= 100000*v; v_t=100000*v_t; v_tt=100000*v_tt;v_ttt=100000*v_ttt;
  double w=0;
  double w_t= 0;
  double w_p=0.;
  double w_tt= 0;
  double w_tp=0.;
  double w_pp=0.;
  double w_ttt= 0; double w_ttp=0; double w_tpp=0; double w_ppp=0;
  
  double phi_p=0.;
	
  
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
  */
  
    //-------------------  atha: FASTER ? --------------------------
  double g_tt = pow(u_t-v,2)+pow(u+v_t,2);
  double g_tt_T = 2*((u_t-v)*(u_tt-v_t)+(u+v_t)*(u_t+v_tt));
  double g_tp=0;
  double g_pp = pow(u*st+v*ct,2);
  double g_pp_T = 2*(u*st+v*ct)*((u_t-v)*st+(u+v_t)*ct);
  double Sg = sqrt(g_tt*g_pp);
  double g = Sg*Sg; double g_T = g_tt_T*g_pp+g_tt*g_pp_T;
  double gtt = 1/g_tt; double gtt_T=-g_tt_T/pow(g_tt,2);
  double gtp=0;
  double gpp = 1/g_pp;
  double denom = Sg/(u*st+v*ct);
  double Rn = (u+v_t)/denom; 
  double Tn = -(u_t-v)/denom;
  double b_tt = ((u+v_t)*(u_tt-2*v_t-u)-(u_t-v)*(v_tt+2*u_t-v))/denom;
  double b_tp=0;
  double b_pp = (-(u+v_t)*(u*st*st+v*ct*st)+(u_t-v)*(u*st*ct+v*ct*ct))/denom;
  double bt_t=gtt*b_tt; double bp_p=gpp*b_pp;
  double btt=gtt*bt_t; double bpp=gpp*bp_p; double btp=0;
  double H = 0.5*(bt_t+bp_p);
  double K = (bt_t*bp_p);
  //---------------Iti :  faster -----------------
  
  double m1=-1.;
  double m2=1.;
  double epsilon=-W2prime(mu,1.0,m1,m2)/(l*(l+1));
  //   Changing Bending Stiffness and epsilon together - YES 
  double scal=1;
    
  double nolispe=1/epsilon;
  ///%%%%%%%%%%%%%%%%%%%%%%%%%%% H-Bending stiffness equal %%%%%%%%%%%%%%%%%%%%%
  double kappa=stiffness(phi,m1,m2,k1,k1); 
  double kappa_prime=stiffness_prime(phi,m1,m2,k1,k1);
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  double c_gp=stiffness_prime(phi,m1,m2,k1,k2);
  
  double gamma_tilde = gamma + W(phi,beta,m1,m2)+epsilon/2*(gtt*phi_t*phi_t+\
							    2*gtp*phi_t*phi_p + gpp*phi_p*phi_p)+Lag_phi*(phi-mu);
  double Norm_eq_part =0;
  double Tgt_eq =0;
  double eps = 1;
  
  Norm_eq_part = ( -2*H*ct*Lag_gz + epsilon*(btt*phi_t*phi_t+2*btp*phi_t*phi_p+bpp*phi_p*phi_p)+ \
		   2*H*kappa*(H*H-K)-2*H*gamma_tilde-pressure+Lag_z*(ct*(Rn-2*H*u)-st*(Tn-2*H*v)));

  Tgt_eq = (gtt_T+0.5*(gtt/g)*g_T-st*ct*gpp+st*Lag_gz);
  //Tgt_eq = (-st*ct*gpp+st*Lag_gz);
  double Lap_herm_0 = gtt*Herm0_tt+(gtt_T+gtt*g_T/(2*g))*Herm0_t;
  double Lap_herm_1 =( gtt*Herm1_tt+(gtt_T+gtt*g_T/(2*g))*Herm1_t );
  double Lap_herm_2 = gtt*Herm2_tt+(gtt_T+gtt*g_T/(2*g))*Herm2_t;
  double Lap_herm_3 =( gtt*Herm3_tt+(gtt_T+gtt*g_T/(2*g))*Herm3_t );

  
  out[0]= (kappa*H*Lap_herm_0 -c_gp*(2*H*gtt-btt)*phi_t*Herm0_t + Norm_eq_part*Herm0)*Sg*dtdx*eps;
  out[1]= (kappa*H*Lap_herm_1 -c_gp*(2*H*gtt-btt)*phi_t*Herm1_t + Norm_eq_part*Herm1)*Sg*dtdx*eps; 
  out[2]= (kappa*H*Lap_herm_2 -c_gp*(2*H*gtt-btt)*phi_t*Herm2_t + Norm_eq_part*Herm2)*Sg*dtdx*eps; 
  out[3]= (kappa*H*Lap_herm_3 -c_gp*(2*H*gtt-btt)*phi_t*Herm3_t +Norm_eq_part*Herm3)*Sg*dtdx*eps; 
  
  int scale = 1;
  out[4]= scale*(-gtt*Herm0_t*0+(Tgt_eq)*Herm0)*Sg*dtdx;
  out[5]= scale*(-gtt*Herm1_t*0+(Tgt_eq)*Herm1)*Sg*dtdx;
  out[6]= scale*(-gtt*Herm2_t*0+(Tgt_eq)*Herm2)*Sg*dtdx;
  out[7]= scale*(-gtt*Herm3_t*0+(Tgt_eq)*Herm3)*Sg*dtdx;
  
  out[8]= nolispe*( epsilon*gtt*phi_t*Herm0_t+(Wprime(phi,beta,m1,m2)+Lag_phi+kappa_prime*H*H+c_gp*K)*Herm0 )*Sg*dtdx;
  out[9]= nolispe*( epsilon*gtt*phi_t*Herm1_t+(Wprime(phi,beta,m1,m2)+Lag_phi+ kappa_prime*H*H +c_gp*K)*Herm1)*Sg*dtdx;
  out[10]= nolispe*( epsilon*gtt*phi_t*Herm2_t+(Wprime(phi,beta,m1,m2)+Lag_phi+ kappa_prime*H*H +c_gp*K)*Herm2)*Sg*dtdx;
  out[11]= nolispe*( epsilon*gtt*phi_t*Herm3_t+(Wprime(phi,beta,m1,m2)+Lag_phi+ kappa_prime*H*H +c_gp*K)*Herm3)*Sg*dtdx;
  
  
  out[12]= (Sg-st)*dtdx*eps; //Area constraint 
  out[13]= (phi-mu)*Sg*dtdx*eps; //Lagrange mult phi constraint 
  out[14]= (u*(ct)-v*st)*Sg*dtdx*eps;//Z- constraint 
  out[15] = (ct)*Sg*dtdx*eps; //Gauge_const

  /*
  //------ Atha : individual hermite: slower
  double Lap_herm_0 = gtt*Hermite(t,0,2,nodes,element)+(gtt_T+gtt*g_T/(2*g))*Hermite(t,0,1,nodes,element);
  double Lap_herm_1 =( gtt*Hermite(t,1,2,nodes,element)+(gtt_T+gtt*g_T/(2*g))*Hermite(t,1,1,nodes,element) );
  double Lap_herm_2 = gtt*Hermite(t,2,2,nodes,element)+(gtt_T+gtt*g_T/(2*g))*Hermite(t,2,1,nodes,element);
  double Lap_herm_3 =( gtt*Hermite(t,3,2,nodes,element)+(gtt_T+gtt*g_T/(2*g))*Hermite(t,3,1,nodes,element) );

  out[0]= (k1*H*Lap_herm_0 + Norm_eq_part*Hermite(t,0,0,nodes,element))*Sg*dtdx*eps;
  out[1]= (k1*H*Lap_herm_1 + Norm_eq_part*Hermite(t,1,0,nodes,element))*Sg*dtdx*eps; 
  out[2]= (k1*H*Lap_herm_2 + Norm_eq_part*Hermite(t,2,0,nodes,element))*Sg*dtdx*eps; 
  out[3]= (k1*H*Lap_herm_3 + Norm_eq_part*Hermite(t,3,0,nodes,element))*Sg*dtdx*eps; 
  
  int scale = 1e5;
  out[4]= scale *(Tgt_eq)*Hermite(t,0,0,nodes,element)*Sg*dtdx;
  out[5]=  scale* (Tgt_eq)*Hermite(t,1,0,nodes,element)*Sg*dtdx;
  out[6]= scale* (Tgt_eq)*Hermite(t,2,0,nodes,element)*Sg*dtdx;
  out[7]=  scale* (Tgt_eq)*Hermite(t,3,0,nodes,element)*Sg*dtdx;
  
  out[8]= ( epsilon*gtt*phi_t*Hermite(t,0,1,nodes,element)+(Wprime(phi,beta,m1,m2)+Lag_phi)*Hermite(t,0,0,nodes,element) )*Sg*dtdx;
  out[9]= ( epsilon*gtt*phi_t*Hermite(t,1,1,nodes,element)+(Wprime(phi,beta,m1,m2)+Lag_phi)*Hermite(t,1,0,nodes,element) )*Sg*dtdx;
  out[10]= ( epsilon*gtt*phi_t*Hermite(t,2,1,nodes,element)+(Wprime(phi,beta,m1,m2)+Lag_phi)*Hermite(t,2,0,nodes,element) )*Sg*dtdx;
  out[11]= ( epsilon*gtt*phi_t*Hermite(t,3,1,nodes,element)+(Wprime(phi,beta,m1,m2)+Lag_phi)*Hermite(t,3,0,nodes,element) )*Sg*dtdx;
  
  
  out[12]= (Sg-st)*dtdx*eps; //Area constraint 
  out[13]= (phi-mu)*Sg*dtdx*eps; //Lagrange mult phi constraint 
  out[14]= (u*(ct)-v*st)*Sg*dtdx*eps;//Z- constraint 
  out[15] = (ct)*Sg*dtdx*eps; //Gauge_const
  //-------iti : individual hermite: slower-------
  */
  return 0; 


}

int firstvariation(double* indvar,double* arg,gsl_vector* out){
  /* This function will be handed over to the continuatino code.
     The arguments cannot be changed
  */
  
  double local[16];
  gsl_vector_set_zero(out);

  int i_v = 2*(NOE);
  int i_phi = 4*(NOE);
  for (int i=0; i<=NOE-1; i++){
    
    local_firstvariation(indvar,arg, i+1 ,local, nodes);
     
    //---------------- Normal Equation -------------------------
    if(i==0){
      out->data[0] += local[0];
      out->data[1] += local[2];
      out->data[2] += local[3];
    }
    else if(i==NOE-1){
      out->data[2*NOE-3] +=local[0];
      out->data[2*NOE-2] +=local[1];
      out->data[2*NOE-1] +=local[2];
      
    }
    else{
      out->data[2*i-1] += local[0] ;
      out->data[2*i] += local[1] ;
      out->data[2*i+1] += local[2] ;
      out->data[2*i+2] += local[3] ;
    }
    
    //--------------------- Tangential ---------------------------
    if(i==0){
      out->data[i_v + 0] += local[5];
      out->data[i_v + 1] += local[6];
      out->data[i_v + 2] += local[7];
    }
    else if(i==NOE-1){
      out->data[i_v + 2*NOE-3] +=local[4];
      out->data[i_v + 2*NOE-2] +=local[5];
      out->data[i_v + 2*NOE-1] +=local[7];
    }
    else{
      out->data[i_v + 2*i-1] += local[4] ;
      out->data[i_v + 2*i] += local[5] ;
      out->data[i_v + 2*i+1] += local[6] ;
      out->data[i_v + 2*i+2] += local[7] ;
    }
    
    //-------------------- Conc ----------------------------
    if(i==0){
      out->data[i_phi + 0] += local[8];
      out->data[i_phi + 1] += local[10];
      out->data[i_phi + 2] += local[11];
    }
    else if(i==NOE-1){
      out->data[i_phi + 2*NOE-3] +=local[8];
      out->data[i_phi + 2*NOE-2] +=local[9];
      out->data[i_phi + 2*NOE-1] +=local[11];
    }
    out->data[i_phi + 2*i-1] += local[8] ;
    out->data[i_phi + 2*i] += local[9] ;
    out->data[i_phi + 2*i+1] += local[10] ;
    out->data[i_phi + 2*i+2] += local[11] ;
    
    
    //Constraints
    out->data[6*NOE] += local[12];
    out->data[6*NOE+1] += local[13];
    out->data[6*NOE+2] += local[14];
    out->data[6*NOE+3] += local[15];

    
  }
 
  return 0;

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
  int phi_end = 6*(NOE);
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
      data_out[2]= data_in[v_start-1];
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
