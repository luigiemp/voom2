#include<iostream>
#include<gsl/gsl_math.h>
#include<gsl/gsl_sf.h>
#include"uvw_model2_FV.h"
#include"uvw_model2_SV.h"
#include<math.h>

using namespace std;
#define PI (3.141592653589793)
#define OUT(i,j) out->data[(i)*out->tda+ (j)]

int local_hessian(double* indvar,double* arg,int element, double hessian[11][11], double* nodes){

  /* Inputs: indvar - the independent variable, theta (colatitude)
     arg - All the unknowns and para
     element - 1,2, ..., NOE (element number)
     - if element = 0, we process constraints
     hessian- is a matrix(double**) that contains the return value. It contains the integrand
     of the weak formulation (ie, normal_eq* H_j, tangt_eq*H_j, conc_eq*N_i)
     for j = 0,1,2,3 and i=0,1 - that is, 'out' holds 8 doubles.
  */
  double x = indvar[0];
  double t = nodes[element-1] + (nodes[element]-nodes[element-1])*x;
  double p = indvar[1];
  double dtdx = nodes[element]-nodes[element-1];
  double h = dtdx;
  int v_start = 2*(NOE); 
  int phi_start = 4*(NOE);
  int Total_len = 6*(NOE);
  
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
	
  double k2= ((arg[Total_len+6]==1)? k1 : k1*(1+arg[Total_len+6]/100.0));    //<<<------------------
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

  double phi_tt = Phi[0]*Hermite(t,0,2,nodes,element)+Phi[1]*Hermite(t,1,2,nodes,element)+Phi[2]*Hermite(t,2,2,nodes,element)+Phi[3]*Hermite(t,3,2,nodes,element);

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
  double m1=-1.;
  double m2=1.;
  double Lap_phi = gtt*phi_tt+(gtt_T+gtt*g_T/(2*g))*phi_t;
  double epsilon=-W2prime(mu,1.0,m1,m2)/(l*(l+1));
  //   Changing Bending Stiffness and epsilon together - YES 
  double scal=1;
  //double r2=sqrt(u*u+v*v+w*w);
  
  double nolispe=1/epsilon;
  double kappa=stiffness(phi,m1,m2,scal*k1,scal*k2); //<<<<------------------------
  double kappa_prime=stiffness_prime(phi,m1,m2,scal*k1,scal*k2);
  
  double gamma_tilde = gamma + W(phi,beta,m1,m2)+epsilon/2*(gtt*phi_t*phi_t+\
							    2*gtp*phi_t*phi_p + gpp*phi_p*phi_p)+Lag_phi*(phi-mu);
  double Norm_eq_part =0;
  double Tgt_eq =0;
  double eps = 1;
  
  Norm_eq_part = ( -2*H*ct*Lag_gz*0 + epsilon*(btt*phi_t*phi_t+2*btp*phi_t*phi_p+bpp*phi_p*phi_p)+ \
		   2*H*k1*(H*H-K)-2*H*gamma_tilde-pressure+Lag_z*(ct*(Rn-2*H*u)-st*(Tn-2*H*v)));
  Tgt_eq = (gtt_T+0.5*(gtt/g)*g_T-st*ct*gpp+st*Lag_gz);
  //Tgt_eq = (v+st*Lag_gz);
  
  
  double Lap_herm_0 = gtt*Hermite(t,0,2,nodes,element)+(gtt_T+gtt*g_T/(2*g))*Hermite(t,0,1,nodes,element);
  double Lap_herm_1 =( gtt*Hermite(t,1,2,nodes,element)+(gtt_T+gtt*g_T/(2*g))*Hermite(t,1,1,nodes,element) );
  double Lap_herm_2 = gtt*Hermite(t,2,2,nodes,element)+(gtt_T+gtt*g_T/(2*g))*Hermite(t,2,1,nodes,element);
  double Lap_herm_3 =( gtt*Hermite(t,3,2,nodes,element)+(gtt_T+gtt*g_T/(2*g))*Hermite(t,3,1,nodes,element) );
  
  double w_i= 0;//r_i(i,t,p);
  double w_i_p=0; //rp_i(i,t,p);
  double w_i_t= 0;
  double w_i_pp=0; //rpp_i(i,t,p);
  double w_i_tp=0; //rzp_i(i,t,p);
  double w_i_tt= 0;
  double psi_i=0; //r modes and phi modes are same
  double psi_i_p=0; //rp_i(i-i0,t,p);
  double psi_i_t=0; //r_z modes and phi_z modes are same
	
  double w_j= 0;//r_i(i,t,p);
  double w_j_p=0; //rp_i(i,t,p);
  double w_j_t= 0;
  double w_j_pp=0; //rpp_i(i,t,p);
  double w_j_tp=0; //rzp_i(i,t,p);
  double w_j_tt= 0;
  double psi_j=0; //r modes and phi modes are same
  double psi_j_p=0; //rp_i(i-i0,t,p);
  double psi_j_t=0; //r_z modes and phi_z modes are same
  double Term_1, Term_2, Term_3;
  double Hess_ij=0;
  double Lap_w_i=0; double Lap_w_j=0;
  double Gammat_tt = 0.5*gtt*g_tt_T;
  double Gammat_pp = -0.5*gtt*g_pp_T;
  double w_i_TT = 0; 
  double w_j_TT = 0; 
  double w_i_PP = 0; 
  double w_j_PP = 0; 
  // ---- atha: may be deleted ------
  double w_i_ttt=0;
  double w_j_ttt=0;
  double Wa_i=0; double Wb_i=0;
  double Wa_i_t=0; double Wb_i_t=0;
  double Wa_i_tt=0; double Wb_i_tt=0;
  double Wa_i_ttt=0; double Wb_i_ttt=0;
  
  double Wa_j=0; double Wb_j=0;
  double Wa_j_t=0; double Wb_j_t=0;
  double Wa_j_tt=0; double Wb_j_tt=0;
  double Wa_j_ttt=0; double Wb_j_ttt=0;
  
  double fact_i = 1;
  double fact_j = 1;
  // may be deleted: iti ------------
  for(int i=0; i<8; i++){
    for(int j=0; j<=i; j++){
      
      fact_i = (i%2==0)? 1 : 1;
      fact_j = (j%2==0)? 1 : 1;
      
      if(i<4) { //------------- Normal Variation ---------------
	w_i =  fact_i*Hermite(t,i,0,nodes,element); 
	w_i_t = fact_i*Hermite(t,i,1,nodes,element); 
	w_i_tt =  fact_i*Hermite(t,i,2,nodes,element); 
	w_i_TT = fact_i*(w_i_tt - Gammat_tt*w_i_t);
	w_i_PP = fact_i*(-Gammat_pp*w_i_t);
	
	// atha: may be deleted 
	w_i_ttt =  Hermite(t,i,3,nodes,element);  //<<<<<<<
	Wa_i=Hermite(nodes[element-1],i,0,nodes,element);
	Wa_i_t=Hermite(nodes[element-1],i,1,nodes,element);
	Wa_i_tt=Hermite(nodes[element-1],i,2,nodes,element);
	Wa_i_ttt=Hermite(nodes[element-1],i,3,nodes,element);
	
	Wb_i=Hermite(nodes[element],i,0,nodes,element);
	Wb_i_t=Hermite(nodes[element],i,1,nodes,element);
	Wb_i_tt=Hermite(nodes[element],i,2,nodes,element);
	Wb_i_ttt=Hermite(nodes[element],i,3,nodes,element);
	// may be deleted: iti
	switch(i){
	case 0:
	  Lap_w_i = fact_i*Lap_herm_0;
	  break;
	case 1:
	  Lap_w_i = fact_i*Lap_herm_1;
	  break;
	case 2:
	  Lap_w_i = fact_i*Lap_herm_2;
	  break;
	case 3:
	  Lap_w_i = fact_i*Lap_herm_3;
	  break;
	}
	psi_i = 0; psi_i_t=0; 
      }
      if(i>=4 && i< 8) {//------- Concentration Variation ------------
	psi_i= Hermite(t,i-4,0,nodes,element);
	psi_i_t = Hermite(t,i-4,1,nodes,element);
	w_i=0; w_i_t=0; w_i_tt=0; w_i_TT=0; w_i_PP=0;
	Lap_w_i=0;
	// atha: may be deleted
	w_i_ttt=0;
	Wa_i=Wa_i_t=Wa_i_tt=Wa_i_ttt=0;
	Wb_i=Wb_i_t=Wb_i_tt=Wb_i_ttt=0;
	// iti: may be deleted
      }
      if(j<4) {//------------- Normal Variation ---------------
	w_j =  fact_j*Hermite(t,j,0,nodes,element); 
	w_j_t = fact_j*Hermite(t,j,1,nodes,element); 
	w_j_tt = fact_j*Hermite(t,j,2,nodes,element); 
	w_j_TT = fact_j*(w_j_tt - Gammat_tt*w_j_t);
	w_j_PP = fact_j*(-Gammat_pp*w_j_t);

	// atha: may be deleted 
	w_j_ttt =  Hermite(t,j,3,nodes,element); 
	Wa_j=Hermite(nodes[element-1],j,0,nodes,element);
	Wa_j_t=Hermite(nodes[element-1],j,1,nodes,element);
	Wa_j_tt=Hermite(nodes[element-1],j,2,nodes,element);
	Wa_j_ttt=Hermite(nodes[element-1],j,3,nodes,element);
	
	Wb_j=Hermite(nodes[element],j,0,nodes,element);
	Wb_j_t=Hermite(nodes[element],j,1,nodes,element);
	Wb_j_tt=Hermite(nodes[element],j,2,nodes,element);
	Wb_j_ttt=Hermite(nodes[element],j,3,nodes,element);
	// may be deleted: iti

	switch(j){
	case 0:
	  Lap_w_j = fact_j*Lap_herm_0;
	  break;
	case 1:
	  Lap_w_j = fact_j*Lap_herm_1;
	  break;
	case 2:
	  Lap_w_j = fact_j*Lap_herm_2;
	  break;
	case 3:
	  Lap_w_j = fact_j*Lap_herm_3;
	  break;
	}
	psi_j=0; psi_j_t=0;
      }
      if(j>=4 && j< 8) {//------- Concentration Variation ------------
	psi_j= Hermite(t,j-4,0,nodes,element);
	psi_j_t = Hermite(t,j-4,1,nodes,element);
	w_j=0; w_j_t=0;w_j_tt=0; w_j_TT=0;w_j_PP=0;
	Lap_w_j=0;

	// atha: may be deleted
	w_j_ttt=0;
	Wa_j=Wa_j_t=Wa_j_tt=Wa_j_ttt=0;
	Wb_j=Wb_j_t=Wb_j_tt=Wb_j_ttt=0;
	// iti: may be deleted
      } 
      			
      /*
      // --------------------- Terms involve H_T : Not accurate for both Helf, Phase field ---------
      Term_1 = 0.5*k1*Lap_w_i*Lap_w_j -((k1*H*H-gamma_tilde)*gtt+ epsilon*gtt*gtt*phi_t*phi_t+2*k1*H*btt)*w_i_t*w_j_t \ 
      + (2*K*gamma_tilde+k1*(8*pow(H,4)-10*pow(H,2)*K+2*pow(K,2))+ \
      2*epsilon*btt*pow(phi_t,2)*H-3*K*epsilon*gtt*pow(phi_t,2)+2*pressure*H)*w_i*w_j;
      Term_2 = 2*(epsilon*btt*phi_t*(w_j*psi_i_t+w_i*psi_j_t) - epsilon*H*gtt*phi_t*(w_j*psi_i_t+w_i*psi_j_t) - \
      H*(w_j*psi_i+w_i*psi_j)*(Lag_phi+Wprime(phi,beta,-1,1)));
      Term_3 = k1*( -2*btt*H_T - 2*H*HT)*(w_i*w_j_t+w_j*w_i_t)-k1*K*(Lap_w_i*w_j+Lap_w_j*w_i); 
      */
      
      
      // -------------------- Terms where H_T is removed by IBP: Accurate for Helf but not Phase field -----
      /*
      Term_1 = 0.5*k1*Lap_w_i*Lap_w_j -((k1*H*H-gamma_tilde)*gtt+ epsilon*gtt*gtt*phi_t*phi_t+2*k1*H*btt)*w_i_t*w_j_t \ 
	+ (2*K*gamma_tilde+k1*(8*pow(H,4)-10*pow(H,2)*K+2*pow(K,2))+ \
	   2*epsilon*btt*pow(phi_t,2)*H-3*K*epsilon*gtt*pow(phi_t,2)+2*pressure*H)*w_i*w_j;
     
      Term_2 = 2*(epsilon*btt*phi_t*(w_j*psi_i_t+w_i*psi_j_t) - epsilon*H*gtt*phi_t*(w_j*psi_i_t+w_i*psi_j_t) - \
		  H*(w_j*psi_i+w_i*psi_j)*(Lag_phi+Wprime(phi,beta,-1,1)));
      
      Term_3 = k1*(-H*H*(w_i*Lap_w_j+w_j*Lap_w_i)+(4*H*btt-2*H*H*gtt)*w_i_t*w_j_t+ \
		   2*H*( btt* (w_i*w_j_TT + w_j*w_i_TT)+bpp*(w_i*w_j_PP + w_j*w_i_PP)))+ \
	-k1*K*(Lap_w_i*w_j+Lap_w_j*w_i); 
      */
      // ------------------------------ No integration by parts --------------
      
	Term_1 = 0.5*kappa*Lap_w_i*Lap_w_j + ( kappa*(7*H*H-2*K) -gamma_tilde)*(w_i*Lap_w_j+w_j*Lap_w_i)/2 + \
	
	(epsilon*gtt*gtt*phi_t*phi_t-4*kappa*H*(2*H*gtt-btt))*(w_i*w_j_TT+w_j*w_i_TT)/2- 4*H*kappa*(2*H*gpp-bpp)*(w_i*w_j_PP+w_j*w_i_PP)/2+ \
	
	2*H*kappa*(btt*w_i_t*w_j_t-H*gtt*w_i_t*w_j_t) +

	  -kappa*H*H*gtt*w_i_t*w_j_t-kappa*H*H*(w_i*Lap_w_j+w_j*Lap_w_i)/2 +//2*H*kappa*HT*(w_i*w_j_t+w_j*w_i_t)/2 +	\

	(kappa*(8*pow(H,4)-10*pow(H,2)*K+2*pow(K,2))+2*K*gamma_tilde+2*H*pressure + \
	     
	epsilon*(2*H*btt*phi_t*phi_t - 3*K*gtt*phi_t*phi_t) )*w_i*w_j;

	Term_2 = 2*(epsilon*btt*phi_t*(w_j*psi_i_t+w_i*psi_j_t) - epsilon*H*gtt*phi_t*(w_j*psi_i_t+w_i*psi_j_t) - \
		    H*(w_j*psi_i+w_i*psi_j)*epsilon*Lap_phi);
		    //H*(w_j*psi_i+w_i*psi_j)*(Lag_phi+Wprime(phi,beta,-1,1)));
      
	Term_3 = 0;
      
	
	//  -H*H*k1*gtt*w_i_t*w_j*t-H*H*k1*(w_i*Lap_w_j+w_j*Lap_w_i)/2+

      //----------------------Second Variation----------------------------------------
		
      if(beta==-1){
	hessian[i][j]= //(psi_i*psi_j +w_i*w_j )*dtdx; 
	  (psi_i*psi_j +w_i*w_j )*Sg*dtdx;
      }
      else{
	hessian[i][j]= //(psi_i_t*psi_j_t +w_i_tt*w_j_tt + 0*(psi_i*w_j_t + psi_j*w_i_t))*dtdx;
	  (Term_1+Term_2 +Term_3 +						\
	   epsilon*gtt*psi_i_t*psi_j_t+W2prime(phi,beta,-1,1)*psi_i*psi_j )*Sg*dtdx;
      }
      hessian[j][i]=hessian[i][j];
      
      
      /*
      // ------------- Mass matrix -----------------------
      hessian[i][j] = (w_i*w_j + psi_i*psi_j)*dtdx;
      hessian[j][i] = hessian[i][j];
      */
      /*
      //----------------------Fist Variation---------------------------------------------
		
      Hess_ij=  (k1*H*Lap_w_i +				\
      ( epsilon*(btt*phi_t*phi_t+2*btp*phi_t*phi_p+bpp*phi_p*phi_p) \
      + 2*H*k1*(H*H-K)-2*H*gamma_tilde-pressure+Lag_z*(ct*(Rn-2*H*u)-st*(Tn-2*H*v)))*w_i\
      + epsilon*gtt*phi_t*psi_i_t+Wprime(phi,beta,-1,1)*psi_i+Lag_phi*psi_i)*Sg/(k1*epsilon);
      gsl_matrix_set(hessian,i,0,Hess_ij);
      */
      //--------Variations must satisfy the following constraints------------------
		
      if(i==j){
	//printf("%e\n",Lag_gz);
	//gsl_matrix_set(hessian,6,i, H*w_i*Sg ); // -------  Area constraint
	//gsl_matrix_set(hessian,7,i, (psi_i-2*H*(phi)*w_i)*Sg ); // Concentration constraint
	hessian[8][i]= H*w_i*Sg*dtdx;
	hessian[9][i]= (psi_i-2*H*(phi)*w_i)*Sg*dtdx;
	hessian[10][i] = w_i*Sg*dtdx;
      }
    }
  }
  return 0;
}
int local_hessian_2constr(double* indvar,double* arg,int element, double hessian[10][10], double* nodes){

  /* Inputs: indvar - the independent variable, theta (colatitude)
     arg - All the unknowns and para
     element - 1,2, ..., NOE (element number)
     - if element = 0, we process constraints
     hessian- is a matrix(double**) that contains the return value. It contains the integrand
     of the weak formulation (ie, normal_eq* H_j, tangt_eq*H_j, conc_eq*N_i)
     for j = 0,1,2,3 and i=0,1 - that is, 'out' holds 8 doubles.
  */
  double x = indvar[0];
  double t = nodes[element-1] + (nodes[element]-nodes[element-1])*x;
  double p = indvar[1];
  double dtdx = nodes[element]-nodes[element-1];
  double h = dtdx;
  int v_start = 2*(NOE); 
  int phi_start = 4*(NOE);
  int Total_len = 6*(NOE);
  
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
	
  double k2= ((arg[Total_len+6]==1)? k1 : k1*(1+arg[Total_len+6]/100.0));    //<<<------------------
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

  double phi_tt = Phi[0]*Hermite(t,0,2,nodes,element)+Phi[1]*Hermite(t,1,2,nodes,element)+Phi[2]*Hermite(t,2,2,nodes,element)+Phi[3]*Hermite(t,3,2,nodes,element);

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
  double m1=-1.;
  double m2=1.;
  double Lap_phi = gtt*phi_tt+(gtt_T+gtt*g_T/(2*g))*phi_t;
  double epsilon=-W2prime(mu,1.0,m1,m2)/(l*(l+1));
  //   Changing Bending Stiffness and epsilon together - YES 
  double scal=1;
  //double r2=sqrt(u*u+v*v+w*w);
  
  double nolispe=1/epsilon;
  double kappa=stiffness(phi,m1,m2,scal*k1,scal*k2); //<<<<------------------------
  double kappa_prime=stiffness_prime(phi,m1,m2,scal*k1,scal*k2);
  
  double gamma_tilde = gamma + W(phi,beta,m1,m2)+epsilon/2*(gtt*phi_t*phi_t+\
							    2*gtp*phi_t*phi_p + gpp*phi_p*phi_p)+Lag_phi*(phi-mu);
  double Norm_eq_part =0;
  double Tgt_eq =0;
  double eps = 1;
  
  Norm_eq_part = ( -2*H*ct*Lag_gz*0 + epsilon*(btt*phi_t*phi_t+2*btp*phi_t*phi_p+bpp*phi_p*phi_p)+ \
		   2*H*k1*(H*H-K)-2*H*gamma_tilde-pressure+Lag_z*(ct*(Rn-2*H*u)-st*(Tn-2*H*v)));
  Tgt_eq = (gtt_T+0.5*(gtt/g)*g_T-st*ct*gpp+st*Lag_gz);
  //Tgt_eq = (v+st*Lag_gz);
  
  
  double Lap_herm_0 = gtt*Hermite(t,0,2,nodes,element)+(gtt_T+gtt*g_T/(2*g))*Hermite(t,0,1,nodes,element);
  double Lap_herm_1 =( gtt*Hermite(t,1,2,nodes,element)+(gtt_T+gtt*g_T/(2*g))*Hermite(t,1,1,nodes,element) );
  double Lap_herm_2 = gtt*Hermite(t,2,2,nodes,element)+(gtt_T+gtt*g_T/(2*g))*Hermite(t,2,1,nodes,element);
  double Lap_herm_3 =( gtt*Hermite(t,3,2,nodes,element)+(gtt_T+gtt*g_T/(2*g))*Hermite(t,3,1,nodes,element) );
  
  double w_i= 0;//r_i(i,t,p);
  double w_i_p=0; //rp_i(i,t,p);
  double w_i_t= 0;
  double w_i_pp=0; //rpp_i(i,t,p);
  double w_i_tp=0; //rzp_i(i,t,p);
  double w_i_tt= 0;
  double psi_i=0; //r modes and phi modes are same
  double psi_i_p=0; //rp_i(i-i0,t,p);
  double psi_i_t=0; //r_z modes and phi_z modes are same
	
  double w_j= 0;//r_i(i,t,p);
  double w_j_p=0; //rp_i(i,t,p);
  double w_j_t= 0;
  double w_j_pp=0; //rpp_i(i,t,p);
  double w_j_tp=0; //rzp_i(i,t,p);
  double w_j_tt= 0;
  double psi_j=0; //r modes and phi modes are same
  double psi_j_p=0; //rp_i(i-i0,t,p);
  double psi_j_t=0; //r_z modes and phi_z modes are same
  double Term_1, Term_2, Term_3;
  double Hess_ij=0;
  double Lap_w_i=0; double Lap_w_j=0;
  double Gammat_tt = 0.5*gtt*g_tt_T;
  double Gammat_pp = -0.5*gtt*g_pp_T;
  double w_i_TT = 0; 
  double w_j_TT = 0; 
  double w_i_PP = 0; 
  double w_j_PP = 0; 
  // ---- atha: may be deleted ------
  double w_i_ttt=0;
  double w_j_ttt=0;
  double Wa_i=0; double Wb_i=0;
  double Wa_i_t=0; double Wb_i_t=0;
  double Wa_i_tt=0; double Wb_i_tt=0;
  double Wa_i_ttt=0; double Wb_i_ttt=0;
  
  double Wa_j=0; double Wb_j=0;
  double Wa_j_t=0; double Wb_j_t=0;
  double Wa_j_tt=0; double Wb_j_tt=0;
  double Wa_j_ttt=0; double Wb_j_ttt=0;
  
  double fact_i = 1;
  double fact_j = 1;
  // may be deleted: iti ------------
  for(int i=0; i<8; i++){
    for(int j=0; j<=i; j++){
      
      fact_i = (i%2==0)? 1 : 1;
      fact_j = (j%2==0)? 1 : 1;
      
      if(i<4) { //------------- Normal Variation ---------------
	w_i =  fact_i*Hermite(t,i,0,nodes,element); 
	w_i_t = fact_i*Hermite(t,i,1,nodes,element); 
	w_i_tt =  fact_i*Hermite(t,i,2,nodes,element); 
	w_i_TT = fact_i*(w_i_tt - Gammat_tt*w_i_t);
	w_i_PP = fact_i*(-Gammat_pp*w_i_t);
	
	// atha: may be deleted 
	w_i_ttt =  Hermite(t,i,3,nodes,element);  //<<<<<<<
	Wa_i=Hermite(nodes[element-1],i,0,nodes,element);
	Wa_i_t=Hermite(nodes[element-1],i,1,nodes,element);
	Wa_i_tt=Hermite(nodes[element-1],i,2,nodes,element);
	Wa_i_ttt=Hermite(nodes[element-1],i,3,nodes,element);
	
	Wb_i=Hermite(nodes[element],i,0,nodes,element);
	Wb_i_t=Hermite(nodes[element],i,1,nodes,element);
	Wb_i_tt=Hermite(nodes[element],i,2,nodes,element);
	Wb_i_ttt=Hermite(nodes[element],i,3,nodes,element);
	// may be deleted: iti
	switch(i){
	case 0:
	  Lap_w_i = fact_i*Lap_herm_0;
	  break;
	case 1:
	  Lap_w_i = fact_i*Lap_herm_1;
	  break;
	case 2:
	  Lap_w_i = fact_i*Lap_herm_2;
	  break;
	case 3:
	  Lap_w_i = fact_i*Lap_herm_3;
	  break;
	}
	psi_i = 0; psi_i_t=0; 
      }
      if(i>=4 && i< 8) {//------- Concentration Variation ------------
	psi_i= Hermite(t,i-4,0,nodes,element);
	psi_i_t = Hermite(t,i-4,1,nodes,element);
	w_i=0; w_i_t=0; w_i_tt=0; w_i_TT=0; w_i_PP=0;
	Lap_w_i=0;
	// atha: may be deleted
	w_i_ttt=0;
	Wa_i=Wa_i_t=Wa_i_tt=Wa_i_ttt=0;
	Wb_i=Wb_i_t=Wb_i_tt=Wb_i_ttt=0;
	// iti: may be deleted
      }
      if(j<4) {//------------- Normal Variation ---------------
	w_j =  fact_j*Hermite(t,j,0,nodes,element); 
	w_j_t = fact_j*Hermite(t,j,1,nodes,element); 
	w_j_tt = fact_j*Hermite(t,j,2,nodes,element); 
	w_j_TT = fact_j*(w_j_tt - Gammat_tt*w_j_t);
	w_j_PP = fact_j*(-Gammat_pp*w_j_t);

	// atha: may be deleted 
	w_j_ttt =  Hermite(t,j,3,nodes,element); 
	Wa_j=Hermite(nodes[element-1],j,0,nodes,element);
	Wa_j_t=Hermite(nodes[element-1],j,1,nodes,element);
	Wa_j_tt=Hermite(nodes[element-1],j,2,nodes,element);
	Wa_j_ttt=Hermite(nodes[element-1],j,3,nodes,element);
	
	Wb_j=Hermite(nodes[element],j,0,nodes,element);
	Wb_j_t=Hermite(nodes[element],j,1,nodes,element);
	Wb_j_tt=Hermite(nodes[element],j,2,nodes,element);
	Wb_j_ttt=Hermite(nodes[element],j,3,nodes,element);
	// may be deleted: iti

	switch(j){
	case 0:
	  Lap_w_j = fact_j*Lap_herm_0;
	  break;
	case 1:
	  Lap_w_j = fact_j*Lap_herm_1;
	  break;
	case 2:
	  Lap_w_j = fact_j*Lap_herm_2;
	  break;
	case 3:
	  Lap_w_j = fact_j*Lap_herm_3;
	  break;
	}
	psi_j=0; psi_j_t=0;
      }
      if(j>=4 && j< 8) {//------- Concentration Variation ------------
	psi_j= Hermite(t,j-4,0,nodes,element);
	psi_j_t = Hermite(t,j-4,1,nodes,element);
	w_j=0; w_j_t=0;w_j_tt=0; w_j_TT=0;w_j_PP=0;
	Lap_w_j=0;

	// atha: may be deleted
	w_j_ttt=0;
	Wa_j=Wa_j_t=Wa_j_tt=Wa_j_ttt=0;
	Wb_j=Wb_j_t=Wb_j_tt=Wb_j_ttt=0;
	// iti: may be deleted
      } 
      			
      /*
      // --------------------- Terms involve H_T : Not accurate for both Helf, Phase field ---------
      Term_1 = 0.5*k1*Lap_w_i*Lap_w_j -((k1*H*H-gamma_tilde)*gtt+ epsilon*gtt*gtt*phi_t*phi_t+2*k1*H*btt)*w_i_t*w_j_t \ 
      + (2*K*gamma_tilde+k1*(8*pow(H,4)-10*pow(H,2)*K+2*pow(K,2))+ \
      2*epsilon*btt*pow(phi_t,2)*H-3*K*epsilon*gtt*pow(phi_t,2)+2*pressure*H)*w_i*w_j;
      Term_2 = 2*(epsilon*btt*phi_t*(w_j*psi_i_t+w_i*psi_j_t) - epsilon*H*gtt*phi_t*(w_j*psi_i_t+w_i*psi_j_t) - \
      H*(w_j*psi_i+w_i*psi_j)*(Lag_phi+Wprime(phi,beta,-1,1)));
      Term_3 = k1*( -2*btt*H_T - 2*H*HT)*(w_i*w_j_t+w_j*w_i_t)-k1*K*(Lap_w_i*w_j+Lap_w_j*w_i); 
      */
      
      
      // -------------------- Terms where H_T is removed by IBP: Accurate for Helf but not Phase field -----
      /*
      Term_1 = 0.5*k1*Lap_w_i*Lap_w_j -((k1*H*H-gamma_tilde)*gtt+ epsilon*gtt*gtt*phi_t*phi_t+2*k1*H*btt)*w_i_t*w_j_t \ 
	+ (2*K*gamma_tilde+k1*(8*pow(H,4)-10*pow(H,2)*K+2*pow(K,2))+ \
	   2*epsilon*btt*pow(phi_t,2)*H-3*K*epsilon*gtt*pow(phi_t,2)+2*pressure*H)*w_i*w_j;
     
      Term_2 = 2*(epsilon*btt*phi_t*(w_j*psi_i_t+w_i*psi_j_t) - epsilon*H*gtt*phi_t*(w_j*psi_i_t+w_i*psi_j_t) - \
		  H*(w_j*psi_i+w_i*psi_j)*(Lag_phi+Wprime(phi,beta,-1,1)));
      
      Term_3 = k1*(-H*H*(w_i*Lap_w_j+w_j*Lap_w_i)+(4*H*btt-2*H*H*gtt)*w_i_t*w_j_t+ \
		   2*H*( btt* (w_i*w_j_TT + w_j*w_i_TT)+bpp*(w_i*w_j_PP + w_j*w_i_PP)))+ \
	-k1*K*(Lap_w_i*w_j+Lap_w_j*w_i); 
      */
      // ------------------------------ No integration by parts --------------
      
	Term_1 = 0.5*kappa*Lap_w_i*Lap_w_j + ( kappa*(7*H*H-2*K) -gamma_tilde)*(w_i*Lap_w_j+w_j*Lap_w_i)/2 + \
	
	(epsilon*gtt*gtt*phi_t*phi_t-4*kappa*H*(2*H*gtt-btt))*(w_i*w_j_TT+w_j*w_i_TT)/2- 4*H*kappa*(2*H*gpp-bpp)*(w_i*w_j_PP+w_j*w_i_PP)/2+ \
	
	2*H*kappa*(btt*w_i_t*w_j_t-H*gtt*w_i_t*w_j_t) +

	  -kappa*H*H*gtt*w_i_t*w_j_t-kappa*H*H*(w_i*Lap_w_j+w_j*Lap_w_i)/2 +//2*H*kappa*HT*(w_i*w_j_t+w_j*w_i_t)/2 +	\

	(kappa*(8*pow(H,4)-10*pow(H,2)*K+2*pow(K,2))+2*K*gamma_tilde+2*H*pressure + \
	     
	epsilon*(2*H*btt*phi_t*phi_t - 3*K*gtt*phi_t*phi_t) )*w_i*w_j;

	Term_2 = 2*(epsilon*btt*phi_t*(w_j*psi_i_t+w_i*psi_j_t) - epsilon*H*gtt*phi_t*(w_j*psi_i_t+w_i*psi_j_t) - \
		    H*(w_j*psi_i+w_i*psi_j)*epsilon*Lap_phi);
		    //H*(w_j*psi_i+w_i*psi_j)*(Lag_phi+Wprime(phi,beta,-1,1)));
      
	Term_3 = 0;
      
	
	//  -H*H*k1*gtt*w_i_t*w_j*t-H*H*k1*(w_i*Lap_w_j+w_j*Lap_w_i)/2+

      //----------------------Second Variation----------------------------------------
		
      if(beta==-1){
	hessian[i][j]= //(psi_i*psi_j +w_i*w_j )*dtdx; 
	  (psi_i*psi_j +w_i*w_j )*Sg*dtdx;
      }
      else{
	hessian[i][j]= //(psi_i_t*psi_j_t +w_i_tt*w_j_tt + 0*(psi_i*w_j_t + psi_j*w_i_t))*dtdx;
	  (Term_1+Term_2 +Term_3 +						\
	   epsilon*gtt*psi_i_t*psi_j_t+W2prime(phi,beta,-1,1)*psi_i*psi_j )*Sg*dtdx;
      }
      hessian[j][i]=hessian[i][j];
      
      
      /*
      // ------------- Mass matrix -----------------------
      hessian[i][j] = (w_i*w_j + psi_i*psi_j)*dtdx;
      hessian[j][i] = hessian[i][j];
      */
      /*
      //----------------------Fist Variation---------------------------------------------
		
      Hess_ij=  (k1*H*Lap_w_i +				\
      ( epsilon*(btt*phi_t*phi_t+2*btp*phi_t*phi_p+bpp*phi_p*phi_p) \
      + 2*H*k1*(H*H-K)-2*H*gamma_tilde-pressure+Lag_z*(ct*(Rn-2*H*u)-st*(Tn-2*H*v)))*w_i\
      + epsilon*gtt*phi_t*psi_i_t+Wprime(phi,beta,-1,1)*psi_i+Lag_phi*psi_i)*Sg/(k1*epsilon);
      gsl_matrix_set(hessian,i,0,Hess_ij);
      */
      //--------Variations must satisfy the following constraints------------------
		
      if(i==j){
	//printf("%e\n",Lag_gz);
	//gsl_matrix_set(hessian,6,i, H*w_i*Sg ); // -------  Area constraint
	//gsl_matrix_set(hessian,7,i, (psi_i-2*H*(phi)*w_i)*Sg ); // Concentration constraint
	hessian[8][i]= H*w_i*Sg*dtdx;
	hessian[9][i]= (psi_i-2*H*(phi)*w_i)*Sg*dtdx;
	//hessian[10][i] = w_i*Sg*dtdx;
      }
    }
  }
  return 0;
}
int local_hessian_nonaxi(double* indvar,double* arg,int element, double hessian[11][11], double* nodes){

  /* Inputs: indvar - the independent variable, theta (colatitude)
     arg - All the unknowns and para
     element - 1,2, ..., NOE (element number)
     - if element = 0, we process constraints
     hessian- is a matrix(double**) that contains the return value. It contains the integrand
     of the weak formulation (ie, normal_eq* H_j, tangt_eq*H_j, conc_eq*N_i)
     for j = 0,1,2,3 and i=0,1 - that is, 'out' holds 8 doubles.
  */

  float m = 1; // Non axisymmetric mode
  double x = indvar[0];
  double t = nodes[element-1] + (nodes[element]-nodes[element-1])*x;
  double p = indvar[1];
  double dtdx = nodes[element]-nodes[element-1];
  double h = dtdx;
  int v_start = 2*(NOE); 
  int phi_start = 4*(NOE);
  int Total_len = 6*NOE;;
  
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
	
  double k2=( (arg[Total_len+6]==1)? k1 : k1*(1+arg[Total_len+6]/100.0) ); 
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
  double m1=-1.;
  double m2=1.;
  double epsilon=-W2prime(mu,1.0,m1,m2)/(l*(l+1));
  //   Changing Bending Stiffness and epsilon together - YES 
  double scal=1;
  //double r2=sqrt(u*u+v*v+w*w);
  
  double nolispe=1/epsilon;
  double kappa=stiffness(phi,m1,m2,scal*k1,scal*k2); //<<<<------------------------
  double kappa_prime=stiffness_prime(phi,m1,m2,scal*k1,scal*k2);
  
  double gamma_tilde = gamma + W(phi,beta,m1,m2)+epsilon/2*(gtt*phi_t*phi_t+\
							    2*gtp*phi_t*phi_p + gpp*phi_p*phi_p)+Lag_phi*(phi-mu);
  double Norm_eq_part =0;
  double Tgt_eq =0;
  double eps = 1;
  
  Norm_eq_part = ( -2*H*ct*Lag_gz*0 + epsilon*(btt*phi_t*phi_t+2*btp*phi_t*phi_p+bpp*phi_p*phi_p)+ \
		   2*H*k1*(H*H-K)-2*H*gamma_tilde-pressure+Lag_z*(ct*(Rn-2*H*u)-st*(Tn-2*H*v)));
  Tgt_eq = (gtt_T+0.5*(gtt/g)*g_T-st*ct*gpp+st*Lag_gz);
  //Tgt_eq = (v+st*Lag_gz);
  
  
  double Lap_herm_0 = gtt*Hermite(t,0,2,nodes,element)+(gtt_T+gtt*g_T/(2*g))*Hermite(t,0,1,nodes,element);
  double Lap_herm_1 =( gtt*Hermite(t,1,2,nodes,element)+(gtt_T+gtt*g_T/(2*g))*Hermite(t,1,1,nodes,element) );
  double Lap_herm_2 = gtt*Hermite(t,2,2,nodes,element)+(gtt_T+gtt*g_T/(2*g))*Hermite(t,2,1,nodes,element);
  double Lap_herm_3 =( gtt*Hermite(t,3,2,nodes,element)+(gtt_T+gtt*g_T/(2*g))*Hermite(t,3,1,nodes,element) );
  
  double w_i= 0;//r_i(i,t,p);
  double w_i_p=0; //rp_i(i,t,p);
  double w_i_t= 0;
  double w_i_pp=0; //rpp_i(i,t,p);
  double w_i_tp=0; //rzp_i(i,t,p);
  double w_i_tt= 0;
  double psi_i=0; //r modes and phi modes are same
  double psi_i_p=0; //rp_i(i-i0,t,p);
  double psi_i_t=0; //r_z modes and phi_z modes are same
	
  double w_j= 0;//r_i(i,t,p);
  double w_j_p=0; //rp_i(i,t,p);
  double w_j_t= 0;
  double w_j_pp=0; //rpp_i(i,t,p);
  double w_j_tp=0; //rzp_i(i,t,p);
  double w_j_tt= 0;
  double psi_j=0; //r modes and phi modes are same
  double psi_j_p=0; //rp_i(i-i0,t,p);
  double psi_j_t=0; //r_z modes and phi_z modes are same
  double Term_1, Term_2, Term_3;
  double Hess_ij=0;
  double Lap_w_i=0; double Lap_w_j=0;
  double Gammat_tt = 0.5*gtt*g_tt_T;
  double Gammat_pp = -0.5*gtt*g_pp_T;
  double w_i_TT = 0; 
  double w_j_TT = 0; 
  double w_i_PP = 0; 
  double w_j_PP = 0; 

  // ---- atha: may be deleted ------
  double w_i_ttt=0;
  double w_j_ttt=0;
  double Wa_i=0; double Wb_i=0;
  double Wa_i_t=0; double Wb_i_t=0;
  double Wa_i_tt=0; double Wb_i_tt=0;
  double Wa_i_ttt=0; double Wb_i_ttt=0;
  
  double Wa_j=0; double Wb_j=0;
  double Wa_j_t=0; double Wb_j_t=0;
  double Wa_j_tt=0; double Wb_j_tt=0;
  double Wa_j_ttt=0; double Wb_j_ttt=0;
  
  // may be deleted: iti ------------

  for(int i=0; i<8; i++){
    for(int j=0; j<=i; j++){
      if(i<4) { //------------- Normal Variation ---------------
	w_i =  Hermite(t,i,0,nodes,element); 
	w_i_t = Hermite(t,i,1,nodes,element); 
	w_i_tt =  Hermite(t,i,2,nodes,element); 
	w_i_TT = w_i_tt - Gammat_tt*w_i_t;
	w_i_PP = -Gammat_pp*w_i_t;
	
	// atha: may be deleted 
	w_i_ttt =  Hermite(t,i,3,nodes,element);  //<<<<<<<
	Wa_i=Hermite(nodes[element-1],i,0,nodes,element);
	Wa_i_t=Hermite(nodes[element-1],i,1,nodes,element);
	Wa_i_tt=Hermite(nodes[element-1],i,2,nodes,element);
	Wa_i_ttt=Hermite(nodes[element-1],i,3,nodes,element);
	
	Wb_i=Hermite(nodes[element],i,0,nodes,element);
	Wb_i_t=Hermite(nodes[element],i,1,nodes,element);
	Wb_i_tt=Hermite(nodes[element],i,2,nodes,element);
	Wb_i_ttt=Hermite(nodes[element],i,3,nodes,element);
	// may be deleted: iti
	switch(i){
	case 0:
	  Lap_w_i = Lap_herm_0;
	  break;
	case 1:
	  Lap_w_i = Lap_herm_1;
	  break;
	case 2:
	  Lap_w_i = Lap_herm_2;
	  break;
	case 3:
	  Lap_w_i = Lap_herm_3;
	  break;
	}
	psi_i = 0; psi_i_t=0; 
      }
      if(i>=4 && i< 8) {//------- Concentration Variation ------------
	psi_i= Hermite(t,i-4,0,nodes,element);
	psi_i_t = Hermite(t,i-4,1,nodes,element);
	w_i=0; w_i_t=0; w_i_tt=0; w_i_TT=0; w_i_PP=0;
	Lap_w_i=0;
	// atha: may be deleted
	w_i_ttt=0;
	Wa_i=Wa_i_t=Wa_i_tt=Wa_i_ttt=0;
	Wb_i=Wb_i_t=Wb_i_tt=Wb_i_ttt=0;
	// iti: may be deleted
      }
      if(j<4) {//------------- Normal Variation ---------------
	w_j =  Hermite(t,j,0,nodes,element); 
	w_j_t = Hermite(t,j,1,nodes,element); 
	w_j_tt =  Hermite(t,j,2,nodes,element); 
	w_j_TT = w_j_tt - Gammat_tt*w_j_t;
	w_j_PP = -Gammat_pp*w_j_t;

	// atha: may be deleted 
	w_j_ttt =  Hermite(t,j,3,nodes,element); 
	Wa_j=Hermite(nodes[element-1],j,0,nodes,element);
	Wa_j_t=Hermite(nodes[element-1],j,1,nodes,element);
	Wa_j_tt=Hermite(nodes[element-1],j,2,nodes,element);
	Wa_j_ttt=Hermite(nodes[element-1],j,3,nodes,element);
	
	Wb_j=Hermite(nodes[element],j,0,nodes,element);
	Wb_j_t=Hermite(nodes[element],j,1,nodes,element);
	Wb_j_tt=Hermite(nodes[element],j,2,nodes,element);
	Wb_j_ttt=Hermite(nodes[element],j,3,nodes,element);
	// may be deleted: iti

	switch(j){
	case 0:
	  Lap_w_j = Lap_herm_0;
	  break;
	case 1:
	  Lap_w_j = Lap_herm_1;
	  break;
	case 2:
	  Lap_w_j = Lap_herm_2;
	  break;
	case 3:
	  Lap_w_j = Lap_herm_3;
	  break;
	}
	psi_j=0; psi_j_t=0;
      }
      if(j>=4 && j< 8) {//------- Concentration Variation ------------
	psi_j= Hermite(t,j-4,0,nodes,element);
	psi_j_t = Hermite(t,j-4,1,nodes,element);
	w_j=0; w_j_t=0;w_j_tt=0; w_j_TT=0;w_j_PP=0;
	Lap_w_j=0;

	// atha: may be deleted
	w_j_ttt=0;
	Wa_j=Wa_j_t=Wa_j_tt=Wa_j_ttt=0;
	Wb_j=Wb_j_t=Wb_j_tt=Wb_j_ttt=0;
	// iti: may be deleted
      } 
      			

      			
      /*
      // --------------------- Terms involve H_T : Not accurate for both Helf, Phase field ---------
      Term_1 = 0.5*k1*Lap_w_i*Lap_w_j -((k1*H*H-gamma_tilde)*gtt+ epsilon*gtt*gtt*phi_t*phi_t+2*k1*H*btt)*w_i_t*w_j_t \ 
      + (2*K*gamma_tilde+k1*(8*pow(H,4)-10*pow(H,2)*K+2*pow(K,2))+ \
      2*epsilon*btt*pow(phi_t,2)*H-3*K*epsilon*gtt*pow(phi_t,2)+2*pressure*H)*w_i*w_j;
      Term_2 = 2*(epsilon*btt*phi_t*(w_j*psi_i_t+w_i*psi_j_t) - epsilon*H*gtt*phi_t*(w_j*psi_i_t+w_i*psi_j_t) - \
      H*(w_j*psi_i+w_i*psi_j)*(Lag_phi+Wprime(phi,beta,-1,1)));
      Term_3 = k1*( -2*btt*H_T - 2*H*HT)*(w_i*w_j_t+w_j*w_i_t)-k1*K*(Lap_w_i*w_j+Lap_w_j*w_i); 
      
      */
      
      // -------------------- Terms where H_T is removed by IBP: Accurate for Helf but not Phase field -----
      /*
	Term_1 = 0.5*k1*Lap_w_i*Lap_w_j -((k1*H*H-gamma_tilde)*gtt+ epsilon*gtt*gtt*phi_t*phi_t+2*k1*H*btt)*w_i_t*w_j_t \ 
	+ (2*K*gamma_tilde+k1*(8*pow(H,4)-10*pow(H,2)*K+2*pow(K,2))+ \
	2*epsilon*btt*pow(phi_t,2)*H-3*K*epsilon*gtt*pow(phi_t,2)+2*pressure*H)*w_i*w_j;
     
	Term_2 = 2*(epsilon*btt*phi_t*(w_j*psi_i_t+w_i*psi_j_t) - epsilon*H*gtt*phi_t*(w_j*psi_i_t+w_i*psi_j_t) - \
	H*(w_j*psi_i+w_i*psi_j)*(Lag_phi+Wprime(phi,beta,-1,1)));
      
	Term_3 = k1*(-H*H*(w_i*Lap_w_j+w_j*Lap_w_i)+(4*H*btt-2*H*H*gtt)*w_i_t*w_j_t+ \
	2*H*( btt* (w_i*w_j_TT + w_j*w_i_TT)+bpp*(w_i*w_j_PP + w_j*w_i_PP)))+ \
	-k1*K*(Lap_w_i*w_j+Lap_w_j*w_i); 
      
      */
      // ------------------------------ No integration by parts --------------
      
      m=1;
      Term_1 = 0.5*kappa*Lap_w_i*Lap_w_j + ( kappa*(7*H*H-2*K) -gamma_tilde)*(w_i*Lap_w_j+w_j*Lap_w_i)/2 + \
	
	(epsilon*gtt*gtt*phi_t*phi_t-4*kappa*H*(2*H*gtt-btt))*(w_i*w_j_TT+w_j*w_i_TT)/2- 4*H*kappa*(2*H*gpp-bpp)*(w_i*w_j_PP+w_j*w_i_PP)/2+ \
	
	2*H*kappa*(btt*w_i_t*w_j_t-H*gtt*w_i_t*w_j_t) +

	-kappa*H*H*gtt*w_i_t*w_j_t-kappa*H*H*(w_i*Lap_w_j+w_j*Lap_w_i)/2 +//+2*H*kappa*HT*(w_i*w_j_t+w_j*w_i_t)/2 + \

	(kappa*(8*pow(H,4)-10*pow(H,2)*K+2*pow(K,2))+2*K*gamma_tilde+2*H*pressure + \
	     
	 epsilon*(2*H*btt*phi_t*phi_t - 3*K*gtt*phi_t*phi_t) )*w_i*w_j;

      Term_2 = 2*(epsilon*btt*phi_t*(w_j*psi_i_t+w_i*psi_j_t) - epsilon*H*gtt*phi_t*(w_j*psi_i_t+w_i*psi_j_t) - \
		  H*(w_j*psi_i+w_i*psi_j)*(Lag_phi+Wprime(phi,beta,-1,1)));
      
      Term_3 = (  kappa/2.*pow(m,4)*pow(gpp,2) - pow(m,2)*gpp*(kappa*(7*H*H-2*K)-gamma_tilde) + 4*H*kappa*pow(m,2)*(2*H*gpp-bpp) + \
		  2*H*kappa*(bpp-H*gpp)*pow(m,2) )*w_i*w_j - kappa*pow(m,2)*gpp*(w_i*Lap_w_j+w_j*Lap_w_i)/2 + \
	epsilon*gpp*pow(m,2)*psi_i*psi_j + pow(m,2)*kappa*H*H*gpp*w_i*w_j-kappa*H*H*pow(m,2)*gpp*w_i*w_j;
      
      //cout<<Lag_z<<" "<<(ct*(Rn-2*H*u)-st*(Tn-2*H*v))<<endl;
      //  -H*H*kappa*gtt*w_i_t*w_j*t-H*H*kappa*(w_i*Lap_w_j+w_j*Lap_w_i)/2+
      

      
      //----------------------Second Variation----------------------------------------
      if(beta==-1){
	hessian[i][j]= (psi_i*psi_j +w_i*w_j )*Sg*dtdx;
      }
      else{
	      
	hessian[i][j]= 
	  (Term_1+Term_2 +Term_3 +	\
	   (epsilon*gtt*psi_i_t*psi_j_t+W2prime(phi,beta,-1,1)*psi_i*psi_j) )*Sg*dtdx;
      }
      hessian[j][i]=hessian[i][j];

      
      /*
      // ------------- Mass matrix -----------------------
      hessian[i][j] = (w_i*w_j + psi_i*psi_j)*Sg*dtdx;
      hessian[j][i] = hessian[i][j];
      */
      /*
      //----------------------Fist Variation---------------------------------------------
		
      Hess_ij=  (k1*H*Lap_w_i +				\
      ( epsilon*(btt*phi_t*phi_t+2*btp*phi_t*phi_p+bpp*phi_p*phi_p) \
      + 2*H*k1*(H*H-K)-2*H*gamma_tilde-pressure+Lag_z*(ct*(Rn-2*H*u)-st*(Tn-2*H*v)))*w_i\
      + epsilon*gtt*phi_t*psi_i_t+Wprime(phi,beta,-1,1)*psi_i+Lag_phi*psi_i)*Sg/(k1*epsilon);
      gsl_matrix_set(hessian,i,0,Hess_ij);
      */
      //--------Variations must satisfy the following constraints------------------
		
      if(i==j){
	
	//gsl_matrix_set(hessian,6,i, H*w_i*Sg ); // -------  Area constraint
	//gsl_matrix_set(hessian,7,i, (psi_i-2*H*(phi)*w_i)*Sg ); // Concentration constraint
	//cout<<ct<<" "<<w_i<<"----"<<endl;
	hessian[8][i]= H*w_i*Sg;
	hessian[9][i]= (psi_i-2*H*(phi)*w_i)*Sg*dtdx;
	hessian[10][i] = w_i*Sg*dtdx;
      }
    }
  }
  return 0;
}


int Hessian(double* indvar,double* arg,gsl_matrix* out){
  /* This function will be handed over to the continuation code.
     The arguments cannot be changed
  */
  
  double local[11][11];
  
  gsl_matrix_set_zero(out); // initialize the matrix to zero
  int i_phi = 2*(NOE);
  // The first element: Special case, due to BC
  int constr=1;
  
  for (int i=0; i<=NOE-1; i++){
    
    local_hessian(indvar,arg, i+1 ,local, nodes);
    
    //------------------- SPARSE ORDERING ---------
    //---------------- Normal Equation ------------------------- 
    if (i==0){
      OUT(0,0) += local[0][0];
      OUT(0,2) += local[0][2];
      OUT(0,3) += local[0][3]; 
       
      OUT(2,0) += local[2][0];
      OUT(2,2) += local[2][2];
      OUT(2,3) += local[2][3];
       
      OUT(3,0) += local[3][0];
      OUT(3,2) += local[3][2];
      OUT(3,3) += local[3][3];
       
      //------- conc ----------
      OUT(1,1) += local[4][4];
      OUT(1,4) += local[4][6];
      OUT(1,5) += local[4][7];
       
      OUT(4,1) += local[6][4];
      OUT(4,4) += local[6][6];
      OUT(4,5) += local[6][7];

      OUT(5,1) += local[7][4];
      OUT(5,4) += local[7][6];
      OUT(5,5) += local[7][7];
       
    }
    else if (i==NOE-1){
      OUT(4*i-2,4*i-2) += local[0][0];
      OUT(4*i-2,4*i-1) += local[0][1];
      OUT(4*i-2,4*i+2) += local[0][2];

      OUT(4*i-1,4*i-2) += local[1][0];
      OUT(4*i-1,4*i-1) += local[1][1];
      OUT(4*i-1,4*i+2) += local[1][2];

      OUT(4*i+2,4*i-2) += local[2][0];
      OUT(4*i+2,4*i-1) += local[2][1];
      OUT(4*i+2,4*i+2) += local[2][2];

      
      //------- conc --------
      OUT(4*i,4*i) += local[4][4];
      OUT(4*i,4*i+1) += local[4][5];
      OUT(4*i,4*i+3) += local[4][6];
            
      OUT(4*i+1,4*i) += local[5][4];
      OUT(4*i+1,4*i+1) += local[5][5];
      OUT(4*i+1,4*i+3) += local[5][6];
            
      OUT(4*i+3,4*i) += local[6][4];
      OUT(4*i+3,4*i+1) += local[6][5];
      OUT(4*i+3,4*i+3) += local[6][6];
      
    }
    else{
      
      OUT(4*i-2,4*i-2) += local[0][0];
      OUT(4*i-2,4*i-1) += local[0][1];
      OUT(4*i-2,4*i+2) += local[0][2];
      OUT(4*i-2,4*i+3) += local[0][3];
      
      OUT(4*i-1,4*i-2) += local[1][0];
      OUT(4*i-1,4*i-1) += local[1][1];
      OUT(4*i-1,4*i+2) += local[1][2];
      OUT(4*i-1,4*i+3) += local[1][3];
      
      OUT(4*i+2,4*i-2) += local[2][0];
      OUT(4*i+2,4*i-1) += local[2][1];
      OUT(4*i+2,4*i+2) += local[2][2];
      OUT(4*i+2,4*i+3) += local[2][3];
      
      OUT(4*i+3,4*i-2) += local[3][0];
      OUT(4*i+3,4*i-1) += local[3][1];
      OUT(4*i+3,4*i+2) += local[3][2];
      OUT(4*i+3,4*i+3) += local[3][3];
      
      //-------------------- Conc ----------------------------
      OUT(4*i,4*i) += local[4][4];
      OUT(4*i,4*i+1) += local[4][5];
      OUT(4*i,4*i+4) += local[4][6];
      OUT(4*i,4*i+5) += local[4][7];
      
      OUT(4*i+1,4*i) += local[5][4];
      OUT(4*i+1,4*i+1) += local[5][5];
      OUT(4*i+1,4*i+4) += local[5][6];
      OUT(4*i+1,4*i+5) += local[5][7];
      
      OUT(4*i+4,4*i) += local[6][4];
      OUT(4*i+4,4*i+1) += local[6][5];
      OUT(4*i+4,4*i+4) += local[6][6];
      OUT(4*i+4,4*i+5) += local[6][7];
      
      OUT(4*i+5,4*i) += local[7][4];
      OUT(4*i+5,4*i+1) += local[7][5];
      OUT(4*i+5,4*i+4) += local[7][6];
      OUT(4*i+5,4*i+5) += local[7][7];
      
      
    }
    
    //---------------- Coupled ------------------------
    if (i==0){
      OUT(1,0) += local[4][0];
      OUT(1,2) += local[4][2];
      OUT(1,3) += local[4][3]; 
      
      OUT(4,0) += local[6][0];
      OUT(4,2) += local[6][2];
      OUT(4,3) += local[6][3];
      
      OUT(5,0) += local[7][0];
      OUT(5,2) += local[7][2];
      OUT(5,3) += local[7][3];
      
      //----------------
      OUT(0,1) += local[0][4];
      OUT(2,1) += local[2][4];
      OUT(3,1) += local[3][4]; 
      
      OUT(0,4) += local[0][6];
      OUT(2,4) += local[2][6];
      OUT(3,4) += local[3][6];
      
      OUT(0,5) += local[0][7];
      OUT(2,5) += local[2][7];
      OUT(3,5) += local[3][7];
       
    }
    else if (i==NOE-1){
      OUT(4*i,4*i-2) += local[4][0];
      OUT(4*i,4*i-1) += local[4][1];
      OUT(4*i,4*i+2) += local[4][2];

      OUT(4*i+1,4*i-2) += local[5][0];
      OUT(4*i+1,4*i-1) += local[5][1];
      OUT(4*i+1,4*i+2) += local[5][2];

      OUT(4*i+3,4*i-2) += local[6][0];
      OUT(4*i+3,4*i-1) += local[6][1];
      OUT(4*i+3,4*i+2) += local[6][2];

      //---------------------------
      OUT(4*i-2,4*i) += local[0][4];
      OUT(4*i-1,4*i) += local[1][4];
      OUT(4*i+2,4*i) += local[2][4];

      OUT(4*i-2,4*i+1) += local[0][5];
      OUT(4*i-1,4*i+1) += local[1][5];
      OUT(4*i+2,4*i+1) += local[2][5];

      OUT(4*i-2,4*i+3) += local[0][6];
      OUT(4*i-1,4*i+3) += local[1][6];
      OUT(4*i+2,4*i+3) += local[2][6];
      
    }
    else{
      
      OUT(4*i,4*i-2) += local[4][0];
      OUT(4*i,4*i-1) += local[4][1];
      OUT(4*i,4*i+2) += local[4][2];
      OUT(4*i,4*i+3) += local[4][3];
      
      OUT(4*i+1,4*i-2) += local[5][0];
      OUT(4*i+1,4*i-1) += local[5][1];
      OUT(4*i+1,4*i+2) += local[5][2];
      OUT(4*i+1,4*i+3) += local[5][3];
      
      OUT(4*i+4,4*i-2) += local[6][0];
      OUT(4*i+4,4*i-1) += local[6][1];
      OUT(4*i+4,4*i+2) += local[6][2];
      OUT(4*i+4,4*i+3) += local[6][3];
      
      OUT(4*i+5,4*i-2) += local[7][0];
      OUT(4*i+5,4*i-1) += local[7][1];
      OUT(4*i+5,4*i+2) += local[7][2];
      OUT(4*i+5,4*i+3) += local[7][3];
      //---------------------------------
      OUT(4*i-2,4*i) += local[0][4];
      OUT(4*i-1,4*i) += local[1][4];
      OUT(4*i+2,4*i) += local[2][4];
      OUT(4*i+3,4*i) += local[3][4];
      
      OUT(4*i-2,4*i+1) += local[0][5];
      OUT(4*i-1,4*i+1) += local[1][5];
      OUT(4*i+2,4*i+1) += local[2][5];
      OUT(4*i+3,4*i+1) += local[3][5];
      
      OUT(4*i-2,4*i+4) += local[0][6];
      OUT(4*i-1,4*i+4) += local[1][6];
      OUT(4*i+2,4*i+4) += local[2][6];
      OUT(4*i+3,4*i+4) += local[3][6];
      
      OUT(4*i-2,4*i+5) += local[0][7];
      OUT(4*i-1,4*i+5) += local[1][7];
      OUT(4*i+2,4*i+5) += local[2][7];
      OUT(4*i+3,4*i+5) += local[3][7];
      
      
      
    }
    if(constr==1){
      // The last two(three) rows are the constraints area/conc/volume
      //Constraints
      if(i==0){
	// Atha: constrain 1 :
	// ---- w ------
	OUT(4*NOE,0) += local[8][0];
	OUT(4*NOE,2) += local[8][2];
	OUT(4*NOE,3) += local[8][3];
	//----- psi -----
	OUT(4*NOE,1) += local[8][4];
	OUT(4*NOE,4) += local[8][6];
	OUT(4*NOE,5) += local[8][7];
	// Iti: constraint 1
	// Atha: constrain 2 :
	// ---- w ------
	OUT(4*NOE+1,0) += local[9][0];
	OUT(4*NOE+1,2) += local[9][2];
	OUT(4*NOE+1,3) += local[9][3];
	//----- psi -----
	OUT(4*NOE+1,1) += local[9][4];
	OUT(4*NOE+1,4) += local[9][6];
	OUT(4*NOE+1,5) += local[9][7];
	// Iti: constraint 2
	// Atha: constrain 3 :
	// ---- w ------
	OUT(4*NOE+2,0) += local[10][0];
	OUT(4*NOE+2,2) += local[10][2];
	OUT(4*NOE+2,3) += local[10][3];
	//----- psi -----
	OUT(4*NOE+2,1) += local[10][4];
	OUT(4*NOE+2,4) += local[10][6];
	OUT(4*NOE+2,5) += local[10][7];
	// Iti: constraint 3

      }
      else if (i==NOE-1){
	// Atha: constrain 1 :
	// ---- w ------
	OUT(4*NOE,4*i-2) += local[8][0];
	OUT(4*NOE,4*i-1) += local[8][1];
	OUT(4*NOE,4*i+2) += local[8][2];
	
	//----- psi -----
	OUT(4*NOE,4*i) += local[8][4];
	OUT(4*NOE,4*i+1) += local[8][5];
	OUT(4*NOE,4*i+3) += local[8][6];
	
	// Iti: constraint 1
	// Atha: constrain 2 :
	// ---- w ------
	OUT(4*NOE+1,4*i-2) += local[9][0];
	OUT(4*NOE+1,4*i-1) += local[9][1];
	OUT(4*NOE+1,4*i+2) += local[9][2];
	//----- psi -----
	OUT(4*NOE+1,4*i) += local[9][4];
	OUT(4*NOE+1,4*i+1) += local[9][5];
	OUT(4*NOE+1,4*i+3) += local[9][6];
	// Iti: constraint 2
	// Atha: constrain 3 :
	// ---- w ------
	OUT(4*NOE+2,4*i-2) += local[10][0];
	OUT(4*NOE+2,4*i-1) += local[10][1];
	OUT(4*NOE+2,4*i+2) += local[10][2];
	//----- psi -----
	OUT(4*NOE+2,4*i) += local[10][4];
	OUT(4*NOE+2,4*i+1) += local[10][5];
	OUT(4*NOE+2,4*i+3) += local[10][6];
	// Iti: constraint 3

      }
      else{
	// Atha: constrain 1 :
	// ---- w ------
	OUT(4*NOE,4*i-2) += local[8][0];
	OUT(4*NOE,4*i-1) += local[8][1];
	OUT(4*NOE,4*i+2) += local[8][2];
	OUT(4*NOE,4*i+3) += local[8][3];
	//----- psi -----
	OUT(4*NOE,4*i) += local[8][4];
	OUT(4*NOE,4*i+1) += local[8][5];
	OUT(4*NOE,4*i+4) += local[8][6];
	OUT(4*NOE,4*i+5) += local[8][7];
	// Iti: constraint 1
	// Atha: constrain 2 :
	// ---- w ------
	OUT(4*NOE+1,4*i-2) += local[9][0];
	OUT(4*NOE+1,4*i-1) += local[9][1];
	OUT(4*NOE+1,4*i+2) += local[9][2];
	OUT(4*NOE+1,4*i+3) += local[9][3];
	//----- psi -----
	OUT(4*NOE+1,4*i) += local[9][4];
	OUT(4*NOE+1,4*i+1) += local[9][5];
	OUT(4*NOE+1,4*i+4) += local[9][6];
	OUT(4*NOE+1,4*i+5) += local[9][7];
	// Iti: constraint 2
	// Atha: constrain 3 :
	// ---- w ------
	OUT(4*NOE+2,4*i-2) += local[10][0];
	OUT(4*NOE+2,4*i-1) += local[10][1];
	OUT(4*NOE+2,4*i+2) += local[10][2];
	OUT(4*NOE+2,4*i+3) += local[10][3];
	//----- psi -----
	OUT(4*NOE+2,4*i) += local[10][4];
	OUT(4*NOE+2,4*i+1) += local[10][5];
	OUT(4*NOE+2,4*i+4) += local[10][6];
	OUT(4*NOE+2,4*i+5) += local[10][7];
	// Iti: constraint 3

      }
    }
    
    
  }
  return 0;

}

int Hessian_2constr(double* indvar,double* arg,gsl_matrix* out){
  /* This function will be handed over to the continuation code.
     The arguments cannot be changed
  */
  
  double local[10][10];
  
  gsl_matrix_set_zero(out); // initialize the matrix to zero
  int i_phi = 2*(NOE);
  // The first element: Special case, due to BC
  int constr=1;
  
  for (int i=0; i<=NOE-1; i++){
    
    local_hessian_2constr(indvar,arg, i+1 ,local, nodes);
    
    //------------------- SPARSE ORDERING ---------
    //---------------- Normal Equation ------------------------- 
    if (i==0){
      OUT(0,0) += local[0][0];
      OUT(0,2) += local[0][2];
      OUT(0,3) += local[0][3]; 
       
      OUT(2,0) += local[2][0];
      OUT(2,2) += local[2][2];
      OUT(2,3) += local[2][3];
       
      OUT(3,0) += local[3][0];
      OUT(3,2) += local[3][2];
      OUT(3,3) += local[3][3];
       
      //------- conc ----------
      OUT(1,1) += local[4][4];
      OUT(1,4) += local[4][6];
      OUT(1,5) += local[4][7];
       
      OUT(4,1) += local[6][4];
      OUT(4,4) += local[6][6];
      OUT(4,5) += local[6][7];

      OUT(5,1) += local[7][4];
      OUT(5,4) += local[7][6];
      OUT(5,5) += local[7][7];
       
    }
    else if (i==NOE-1){
      OUT(4*i-2,4*i-2) += local[0][0];
      OUT(4*i-2,4*i-1) += local[0][1];
      OUT(4*i-2,4*i+2) += local[0][2];

      OUT(4*i-1,4*i-2) += local[1][0];
      OUT(4*i-1,4*i-1) += local[1][1];
      OUT(4*i-1,4*i+2) += local[1][2];

      OUT(4*i+2,4*i-2) += local[2][0];
      OUT(4*i+2,4*i-1) += local[2][1];
      OUT(4*i+2,4*i+2) += local[2][2];

      
      //------- conc --------
      OUT(4*i,4*i) += local[4][4];
      OUT(4*i,4*i+1) += local[4][5];
      OUT(4*i,4*i+3) += local[4][6];
            
      OUT(4*i+1,4*i) += local[5][4];
      OUT(4*i+1,4*i+1) += local[5][5];
      OUT(4*i+1,4*i+3) += local[5][6];
            
      OUT(4*i+3,4*i) += local[6][4];
      OUT(4*i+3,4*i+1) += local[6][5];
      OUT(4*i+3,4*i+3) += local[6][6];
      
    }
    else{
      
      OUT(4*i-2,4*i-2) += local[0][0];
      OUT(4*i-2,4*i-1) += local[0][1];
      OUT(4*i-2,4*i+2) += local[0][2];
      OUT(4*i-2,4*i+3) += local[0][3];
      
      OUT(4*i-1,4*i-2) += local[1][0];
      OUT(4*i-1,4*i-1) += local[1][1];
      OUT(4*i-1,4*i+2) += local[1][2];
      OUT(4*i-1,4*i+3) += local[1][3];
      
      OUT(4*i+2,4*i-2) += local[2][0];
      OUT(4*i+2,4*i-1) += local[2][1];
      OUT(4*i+2,4*i+2) += local[2][2];
      OUT(4*i+2,4*i+3) += local[2][3];
      
      OUT(4*i+3,4*i-2) += local[3][0];
      OUT(4*i+3,4*i-1) += local[3][1];
      OUT(4*i+3,4*i+2) += local[3][2];
      OUT(4*i+3,4*i+3) += local[3][3];
      
      //-------------------- Conc ----------------------------
      OUT(4*i,4*i) += local[4][4];
      OUT(4*i,4*i+1) += local[4][5];
      OUT(4*i,4*i+4) += local[4][6];
      OUT(4*i,4*i+5) += local[4][7];
      
      OUT(4*i+1,4*i) += local[5][4];
      OUT(4*i+1,4*i+1) += local[5][5];
      OUT(4*i+1,4*i+4) += local[5][6];
      OUT(4*i+1,4*i+5) += local[5][7];
      
      OUT(4*i+4,4*i) += local[6][4];
      OUT(4*i+4,4*i+1) += local[6][5];
      OUT(4*i+4,4*i+4) += local[6][6];
      OUT(4*i+4,4*i+5) += local[6][7];
      
      OUT(4*i+5,4*i) += local[7][4];
      OUT(4*i+5,4*i+1) += local[7][5];
      OUT(4*i+5,4*i+4) += local[7][6];
      OUT(4*i+5,4*i+5) += local[7][7];
      
      
    }
    
    //---------------- Coupled ------------------------
    if (i==0){
      OUT(1,0) += local[4][0];
      OUT(1,2) += local[4][2];
      OUT(1,3) += local[4][3]; 
      
      OUT(4,0) += local[6][0];
      OUT(4,2) += local[6][2];
      OUT(4,3) += local[6][3];
      
      OUT(5,0) += local[7][0];
      OUT(5,2) += local[7][2];
      OUT(5,3) += local[7][3];
      
      //----------------
      OUT(0,1) += local[0][4];
      OUT(2,1) += local[2][4];
      OUT(3,1) += local[3][4]; 
      
      OUT(0,4) += local[0][6];
      OUT(2,4) += local[2][6];
      OUT(3,4) += local[3][6];
      
      OUT(0,5) += local[0][7];
      OUT(2,5) += local[2][7];
      OUT(3,5) += local[3][7];
       
    }
    else if (i==NOE-1){
      OUT(4*i,4*i-2) += local[4][0];
      OUT(4*i,4*i-1) += local[4][1];
      OUT(4*i,4*i+2) += local[4][2];

      OUT(4*i+1,4*i-2) += local[5][0];
      OUT(4*i+1,4*i-1) += local[5][1];
      OUT(4*i+1,4*i+2) += local[5][2];

      OUT(4*i+3,4*i-2) += local[6][0];
      OUT(4*i+3,4*i-1) += local[6][1];
      OUT(4*i+3,4*i+2) += local[6][2];

      //---------------------------
      OUT(4*i-2,4*i) += local[0][4];
      OUT(4*i-1,4*i) += local[1][4];
      OUT(4*i+2,4*i) += local[2][4];

      OUT(4*i-2,4*i+1) += local[0][5];
      OUT(4*i-1,4*i+1) += local[1][5];
      OUT(4*i+2,4*i+1) += local[2][5];

      OUT(4*i-2,4*i+3) += local[0][6];
      OUT(4*i-1,4*i+3) += local[1][6];
      OUT(4*i+2,4*i+3) += local[2][6];
      
    }
    else{
      
      OUT(4*i,4*i-2) += local[4][0];
      OUT(4*i,4*i-1) += local[4][1];
      OUT(4*i,4*i+2) += local[4][2];
      OUT(4*i,4*i+3) += local[4][3];
      
      OUT(4*i+1,4*i-2) += local[5][0];
      OUT(4*i+1,4*i-1) += local[5][1];
      OUT(4*i+1,4*i+2) += local[5][2];
      OUT(4*i+1,4*i+3) += local[5][3];
      
      OUT(4*i+4,4*i-2) += local[6][0];
      OUT(4*i+4,4*i-1) += local[6][1];
      OUT(4*i+4,4*i+2) += local[6][2];
      OUT(4*i+4,4*i+3) += local[6][3];
      
      OUT(4*i+5,4*i-2) += local[7][0];
      OUT(4*i+5,4*i-1) += local[7][1];
      OUT(4*i+5,4*i+2) += local[7][2];
      OUT(4*i+5,4*i+3) += local[7][3];
      //---------------------------------
      OUT(4*i-2,4*i) += local[0][4];
      OUT(4*i-1,4*i) += local[1][4];
      OUT(4*i+2,4*i) += local[2][4];
      OUT(4*i+3,4*i) += local[3][4];
      
      OUT(4*i-2,4*i+1) += local[0][5];
      OUT(4*i-1,4*i+1) += local[1][5];
      OUT(4*i+2,4*i+1) += local[2][5];
      OUT(4*i+3,4*i+1) += local[3][5];
      
      OUT(4*i-2,4*i+4) += local[0][6];
      OUT(4*i-1,4*i+4) += local[1][6];
      OUT(4*i+2,4*i+4) += local[2][6];
      OUT(4*i+3,4*i+4) += local[3][6];
      
      OUT(4*i-2,4*i+5) += local[0][7];
      OUT(4*i-1,4*i+5) += local[1][7];
      OUT(4*i+2,4*i+5) += local[2][7];
      OUT(4*i+3,4*i+5) += local[3][7];
      
      
      
    }
    if(constr==1){
      // The last two(three) rows are the constraints area/conc/volume
      //Constraints
      if(i==0){
	// Atha: constrain 1 :
	// ---- w ------
	OUT(4*NOE,0) += local[8][0];
	OUT(4*NOE,2) += local[8][2];
	OUT(4*NOE,3) += local[8][3];
	//----- psi -----
	OUT(4*NOE,1) += local[8][4];
	OUT(4*NOE,4) += local[8][6];
	OUT(4*NOE,5) += local[8][7];
	// Iti: constraint 1
	// Atha: constrain 2 :
	// ---- w ------
	OUT(4*NOE+1,0) += local[9][0];
	OUT(4*NOE+1,2) += local[9][2];
	OUT(4*NOE+1,3) += local[9][3];
	//----- psi -----
	OUT(4*NOE+1,1) += local[9][4];
	OUT(4*NOE+1,4) += local[9][6];
	OUT(4*NOE+1,5) += local[9][7];
	// Iti: constraint 2
	
      }
      else if (i==NOE-1){
	// Atha: constrain 1 :
	// ---- w ------
	OUT(4*NOE,4*i-2) += local[8][0];
	OUT(4*NOE,4*i-1) += local[8][1];
	OUT(4*NOE,4*i+2) += local[8][2];
	
	//----- psi -----
	OUT(4*NOE,4*i) += local[8][4];
	OUT(4*NOE,4*i+1) += local[8][5];
	OUT(4*NOE,4*i+3) += local[8][6];
	
	// Iti: constraint 1
	// Atha: constrain 2 :
	// ---- w ------
	OUT(4*NOE+1,4*i-2) += local[9][0];
	OUT(4*NOE+1,4*i-1) += local[9][1];
	OUT(4*NOE+1,4*i+2) += local[9][2];
	//----- psi -----
	OUT(4*NOE+1,4*i) += local[9][4];
	OUT(4*NOE+1,4*i+1) += local[9][5];
	OUT(4*NOE+1,4*i+3) += local[9][6];
	// Iti: constraint 2
	
      }
      else{
	// Atha: constrain 1 :
	// ---- w ------
	OUT(4*NOE,4*i-2) += local[8][0];
	OUT(4*NOE,4*i-1) += local[8][1];
	OUT(4*NOE,4*i+2) += local[8][2];
	OUT(4*NOE,4*i+3) += local[8][3];
	//----- psi -----
	OUT(4*NOE,4*i) += local[8][4];
	OUT(4*NOE,4*i+1) += local[8][5];
	OUT(4*NOE,4*i+4) += local[8][6];
	OUT(4*NOE,4*i+5) += local[8][7];
	// Iti: constraint 1
	// Atha: constrain 2 :
	// ---- w ------
	OUT(4*NOE+1,4*i-2) += local[9][0];
	OUT(4*NOE+1,4*i-1) += local[9][1];
	OUT(4*NOE+1,4*i+2) += local[9][2];
	OUT(4*NOE+1,4*i+3) += local[9][3];
	//----- psi -----
	OUT(4*NOE+1,4*i) += local[9][4];
	OUT(4*NOE+1,4*i+1) += local[9][5];
	OUT(4*NOE+1,4*i+4) += local[9][6];
	OUT(4*NOE+1,4*i+5) += local[9][7];
	// Iti: constraint 2
	
      }
    }
    
    
  }
  return 0;

}


int Hessian_nonaxi(double* indvar,double* arg,gsl_matrix* out){
  /* This function will be handed over to the continuation code.
     The arguments cannot be changed
  */
  
  double local[11][11];
  
  gsl_matrix_set_zero(out); // initialize the matrix to zero
  int i_phi = 2*(NOE)+2;
  // The first element: Special case, due to BC
  int constr=1;
  
  for (int i=0; i<=NOE-1; i++){
    
    local_hessian_nonaxi(indvar,arg, i+1 ,local, nodes);
    //cout<<out->size1<<" "<<out->size2<<" "<<i_phi+i+1<<endl;
        
    //------------------- SPARSE ORDERING ---------
    //---------------- Normal Equation ------------------------- 
      
    OUT(4*i,4*i) += local[0][0];
    OUT(4*i,4*i+1) += local[0][1];
    OUT(4*i,4*i+4) += local[0][2];
    OUT(4*i,4*i+5) += local[0][3];
      
    OUT(4*i+1,4*i) += local[1][0];
    OUT(4*i+1,4*i+1) += local[1][1];
    OUT(4*i+1,4*i+4) += local[1][2];
    OUT(4*i+1,4*i+5) += local[1][3];
      
    OUT(4*i+4,4*i) += local[2][0];
    OUT(4*i+4,4*i+1) += local[2][1];
    OUT(4*i+4,4*i+4) += local[2][2];
    OUT(4*i+4,4*i+5) += local[2][3];
      
    OUT(4*i+5,4*i) += local[3][0];
    OUT(4*i+5,4*i+1) += local[3][1];
    OUT(4*i+5,4*i+4) += local[3][2];
    OUT(4*i+5,4*i+5) += local[3][3];
      
    //-------------------- Conc ----------------------------
    OUT(4*i+2,4*i+2) += local[4][4];
    OUT(4*i+2,4*i+3) += local[4][5];
    OUT(4*i+2,4*i+6) += local[4][6];
    OUT(4*i+2,4*i+7) += local[4][7];
      
    OUT(4*i+3,4*i+2) += local[5][4];
    OUT(4*i+3,4*i+3) += local[5][5];
    OUT(4*i+3,4*i+6) += local[5][6];
    OUT(4*i+3,4*i+7) += local[5][7];
      
    OUT(4*i+6,4*i+2) += local[6][4];
    OUT(4*i+6,4*i+3) += local[6][5];
    OUT(4*i+6,4*i+6) += local[6][6];
    OUT(4*i+6,4*i+7) += local[6][7];
      
    OUT(4*i+7,4*i+2) += local[7][4];
    OUT(4*i+7,4*i+3) += local[7][5];
    OUT(4*i+7,4*i+6) += local[7][6];
    OUT(4*i+7,4*i+7) += local[7][7];
      
    
    //---------------- Coupled ------------------------

      
    OUT(4*i+2,4*i) += local[4][0];
    OUT(4*i+2,4*i+1) += local[4][1];
    OUT(4*i+2,4*i+4) += local[4][2];
    OUT(4*i+2,4*i+5) += local[4][3];
      
    OUT(4*i+3,4*i) += local[5][0];
    OUT(4*i+3,4*i+1) += local[5][1];
    OUT(4*i+3,4*i+4) += local[5][2];
    OUT(4*i+3,4*i+5) += local[5][3];
      
    OUT(4*i+6,4*i) += local[6][0];
    OUT(4*i+6,4*i+1) += local[6][1];
    OUT(4*i+6,4*i+4) += local[6][2];
    OUT(4*i+6,4*i+5) += local[6][3];
      
    OUT(4*i+7,4*i) += local[7][0];
    OUT(4*i+7,4*i+1) += local[7][1];
    OUT(4*i+7,4*i+4) += local[7][2];
    OUT(4*i+7,4*i+5) += local[7][3];
    //---------------------------------
    OUT(4*i,4*i+2) += local[0][4];
    OUT(4*i+1,4*i+2) += local[1][4];
    OUT(4*i+4,4*i+2) += local[2][4];
    OUT(4*i+5,4*i+2) += local[3][4];
      
    OUT(4*i,4*i+3) += local[0][5];
    OUT(4*i+1,4*i+3) += local[1][5];
    OUT(4*i+4,4*i+3) += local[2][5];
    OUT(4*i+5,4*i+3) += local[3][5];
      
    OUT(4*i,4*i+6) += local[0][6];
    OUT(4*i+1,4*i+6) += local[1][6];
    OUT(4*i+4,4*i+6) += local[2][6];
    OUT(4*i+5,4*i+6) += local[3][6];
      
    OUT(4*i,4*i+7) += local[0][7];
    OUT(4*i+1,4*i+7) += local[1][7];
    OUT(4*i+4,4*i+7) += local[2][7];
    OUT(4*i+5,4*i+7) += local[3][7];
      
      
      

    if(constr==1){
      // The last two(three) rows are the constraints area/conc/volume
      //Constraints

      // Atha: constrain 1 :
      // ---- w ------
      OUT(4*NOE+4,4*i) += local[8][0];
      OUT(4*NOE+4,4*i+1) += local[8][1];
      OUT(4*NOE+4,4*i+4) += local[8][2];
      OUT(4*NOE+4,4*i+5) += local[8][3];
      //----- psi -----
      OUT(4*NOE+4,4*i+2) += local[8][4];
      OUT(4*NOE+4,4*i+3) += local[8][5];
      OUT(4*NOE+4,4*i+6) += local[8][6];
      OUT(4*NOE+4,4*i+7) += local[8][7];
      // Iti: constraint 1
      // Atha: constrain 2 :
      // ---- w ------
      OUT(4*NOE+5,4*i) += local[9][0];
      OUT(4*NOE+5,4*i+1) += local[9][1];
      OUT(4*NOE+5,4*i+4) += local[9][2];
      OUT(4*NOE+5,4*i+5) += local[9][3];
      //----- psi -----
      OUT(4*NOE+5,4*i+2) += local[9][4];
      OUT(4*NOE+5,4*i+3) += local[9][5];
      OUT(4*NOE+5,4*i+6) += local[9][6];
      OUT(4*NOE+5,4*i+7) += local[9][7];
      // Iti: constraint 2
      // Atha: constrain 3 :
      // ---- w ------
      OUT(4*NOE+6,4*i) += local[10][0];
      OUT(4*NOE+6,4*i+1) += local[10][1];
      OUT(4*NOE+6,4*i+4) += local[10][2];
      OUT(4*NOE+6,4*i+5) += local[10][3];
      //----- psi -----
      OUT(4*NOE+6,4*i+2) += local[10][4];
      OUT(4*NOE+6,4*i+3) += local[10][5];
      OUT(4*NOE+6,4*i+6) += local[10][6];
      OUT(4*NOE+6,4*i+7) += local[10][7];
      // Iti: constraint 2


    }


  }
  
  return 0;

}
