#include<iostream>
#include<gsl/gsl_math.h>
#include<gsl/gsl_sf.h>
#include"uvw_model2_const.h"
extern int PLM_COUNT;
using namespace std;
#define PI (3.141592653589793)

int constraint(double* indvar,double* arg,gsl_matrix* out){
	
  //PLM_COUNT++;
  double t=indvar[0];
  double p=indvar[1];
  int chi_len = phi_len;
  int v_start=u_len; int w_start=v_start+v_len;
  int chi_start=w_start+w_len; //chi represents normal variation
  int Total_len=chi_start+chi_len;
      
  double u_discrete[u_len];
  double v_discrete[v_len];
  double w_discrete[w_len];
  double chi_discrete[chi_len];
  int max_size=max(chi_len,u_len);  // NEEDS TO INCORPORATE U V W CHI
  max_size=max(max_size,v_len);

  double ct=cos(t);
  double st=sin(t);
  double cp=cos(p);
  double sp=sin(p);
	
  Plm_all(max_size+1, 0, ct, plm);    //Axisymmmetric - m=0; +2 needed for derivatives
  Plm_all(max_size+1, 1, ct, plm_1);    //Needed to compute first derivative
  Plm_all(max_size+1, 2, ct, plm_2);    //Needed to compute second derivative
  Plm_all(max_size+1, 3, ct, plm_3);    //Needed to compute third derivative
  Plm_all(max_size+1, 4, ct, plm_4);    //Needed to compute fourth derivative (3rd derv of v,w)

  double gamma=arg[Total_len];
  double Lag_phi=arg[Total_len+1];
  double Lag_z=arg[Total_len+2];
  double Lag_gz = arg[Total_len+3];
  double l=arg[Total_len+4];
  double k1=arg[Total_len+5];
	
  double k2=k1;
  //double k2=arg[Total_len+6];
	
  double beta=arg[Total_len+7];
  double pressure=arg[Total_len+8];
  double mu=arg[Total_len+9];
  // ---------------------Everything that needs plm--------------------------
  double chi=0; 
  double chi_t=0;
  double chi_p=0.;
	
  double u=0; 
  double u_t=0;
  double u_p=0.;
  double u_tt=0;
  double u_tp=0.;
  double u_pp=0.;
  double u_ttt=0; double u_ttp=0; double u_tpp=0; double u_ppp=0;
	
  double v=0; 
  double v_t=0;
  double v_p=0.;
  double v_tt=0;
  double v_tp=0.;
  double v_pp=0.;
  double v_ttt=0; double v_ttp=0; double v_tpp=0; double v_ppp=0;
	
  double w=0; 
  double w_t=0;
  double w_p=0.;
  double w_tt=0;
  double w_tp=0.;
  double w_pp=0.;
  double w_ttt=0; double w_ttp=0; double w_tpp=0; double w_ppp=0;
	
  double flm=0;
  double sqfac_i_1=0;
  double sqfac_i_2=0;
  double sqfac_i_3=0;
  double sqfac_i_4=0;	
  //------------Legendre Poly----------------------
	
  for(int i=0; i< max_size; i++){
    flm=1;
    sqfac_i_1=SQFAC_i_1[i];//sqfac(i,1);
    sqfac_i_2=SQFAC_i_2[i];//sqfac(i,2);
    sqfac_i_3=SQFAC_i_3[i];//sqfac(i,3);
	   
    if (i<u_len-1){
      //printf("%f\n",arg[i+chi_start]);
      /* The normal variation mode is stored at phi's location.
	 But by incorporating area constraint, we have reduced
	 the number of components of normal variation by one. 
	 More specifically, we have removed the zeroth mode.
	 Hence, we need to start from mode 1, plm[i+1]
      */
      chi += plm[i+1]*arg[i+chi_start]; 
      //chi += plm[i+1]*float(i==0);
      //printf("\033[1;33m %f\033[m",arg[i+chi_start]);
    }
    //printf("\n-----\n");
	  
    if(i<u_len){
      u += plm[i]*arg[i];
      u_t += sqfac_i_1* (plm_1[i]) *arg[i];
      u_tt += 0.5*(sqfac_i_2*plm_2[i]-i*(i+1)*plm[i] ) *arg[i];
      u_ttt += -0.25* ((3*i*i+3*i-2)*sqfac_i_1*plm_1[i]-sqfac_i_3*plm_3[i])*arg[i];
	     
    }
	   
    if(i<v_len){
      int k=i+1;
      //flm=sqrt((2*k+1.0)/(4*M_PI));
      sqfac_i_1=SQFAC_i_1[k];//sqfac(k,1);//
      sqfac_i_2=SQFAC_i_2[k];//sqfac(k,2);//
      sqfac_i_3=SQFAC_i_3[k];//sqfac(k,3);//
      sqfac_i_4=SQFAC_i_4[k];;//sqfac(k,4);//
	     
      v += sqfac_i_1* (plm_1[k]) *arg[i+v_start];
      v_t += 0.5*(sqfac_i_2*plm_2[k]-k*(k+1)*plm[k] ) *arg[i+v_start];
      v_tt += -0.25* ((3*k*k+3*k-2)*sqfac_i_1*plm_1[k]-sqfac_i_3*plm_3[k])*arg[i+v_start];
      //v_ttt += flm*(plm_4[k]-4*plm_2[k]+ct/st*(3*ct/st*plm_2[k]+6*plm_3[k]-plm_1[k]))*v_discrete[i];
      v_ttt += (1/8.0)*(sqfac_i_4*plm_4[k]-(4*k*k+4*k-8)*sqfac_i_2*plm_2[k]+k*(k+1)*(3*k*k+3*k-2)*plm[k])*arg[i+v_start];
	     
    }
	   
    if(i<w_len){
      int k=i+1;
      w += sqfac(k,1)* (plm_1[k]) *arg[i+w_start];
      w_t += 0.5*(sqfac(k,2)*plm_2[k]-k*(k+1)*plm[k] ) *arg[i+w_start];
      w_tt += -0.25* ((3*k*k+3*k-2)*sqfac(k,1)*plm_1[k]-sqfac(k,3)*plm_3[k])*arg[i+w_start];
      w_ttt += (1/8.0)*(sqfac(k,4)*plm_4[k]-(4*k*k+4*k-8)*sqfac(k,2)*plm_2[k]+k*(k+1)*(3*k*k+3*k-2)*plm[k])*arg[i+w_start];
    }
	   
  }
  //r=r+1.;
	
  //printf("-----\n");
  u=u+1; // In the following (u,v,w) represents current configuration
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
  /*
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
  */
  // Normal Vector 
  double Rn = (Bt*Cp-Bp*Ct)/Sg;
  double Tn = (Ct*Ap-Cp*At)/Sg;
  double Pn = (At*Bp-Ap*Bt)/Sg;
	
  /*
    double Rn_T = (Bt_T*Cp+Bt*Cp_T-Bp_T*Ct-Bp*Ct_T)/Sg-0.5*Rn*g_T/g;
    double Tn_T = (Ct_T*Ap+Ct*Ap_T-Cp_T*At-Cp*At_T)/Sg-0.5*Tn*g_T/g;
    double Pn_T = (At_T*Bp+At*Bp_T-Ap_T*Bt-Ap*Bt_T)/Sg-0.5*Pn*g_T/g;
	
    double Rn_P = (Bt_P*Cp+Bt*Cp_P-Bp_P*Ct-Bp*Ct_P)/Sg-0.5*Rn*g_P/g;
    double Tn_P = (Ct_P*Ap+Ct*Ap_P-Cp_P*At-Cp*At_P)/Sg-0.5*Tn*g_P/g;
    double Pn_P = (At_P*Bp+At*Bp_P-Ap_P*Bt-Ap*Bt_P)/Sg-0.5*Pn*g_P/g;
  */
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
  /*
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
  */
  //----------------------------------------------------------------
  double H = 0.5*(bt_t+bp_p); double K = (bt_t*bp_p - bt_p*bp_t);
  //printf("\n\n~~~%f~~~\n",Sg/st);
  gsl_matrix_set(out,0,0,chi*Sg);
  return 0; 

}

