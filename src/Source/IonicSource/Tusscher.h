//-*-C++-*-
/*! \brief
  Tusscher el ionic model for human ventricular tissue. Developed based on 
  the paper "Alternans and spiral breakup in a human ventricular tissue model"
  by K. H. W. J. ten Tusscher and A. V. Panfilov, Am J Physiol Heart Circ 
  Physiol 291:H1088-H1100, 2006. First published 24 March 2006
  - Implemented methods described in "Efficient Numerical Technique for the 
  solution of the monodomain and Bidomain Equations" by Whiteleley. We can
  take larger time steps for Ionic Solve
*/
#ifndef _Tusscher_h_
#define _Tusscher_h_

#include "IonicMaterial.h"

namespace voom{
  class Tusscher: public IonicMaterial {
  private:
    //! Constants for the model
    Real *_Constants;
    
    //! State Variables
    std::vector<Real>  _State;

    //! Rate Variables
    Real  _Rates[19];

    //! Algebraic Variables
    Real _Algebraic[70];

    //! Adaptive Update of State Variables
    void UpdateStateVariables(const Real dt);

  public:
    //! Constructor
    Tusscher(Real *constants){
      _Constants = constants;
      _State.resize( 19 );
      /* As in Tusscher Paper
      // State Variables
      _State[0] = -85.423;    _State[1] = 138.52;    _State[2] = 10.132;
      _State[3] = 0.000153;   _State[4] = 0.0165;    _State[5] = 0.473;
      _State[6] = 0.0174;     _State[7] = 0.00165;   _State[8] = 0.749;
      _State[9] = 0.6788;     _State[10] = 0.00042;  _State[11] = 3.288e-5;
      _State[12] = 0.7026;    _State[13] = 0.9526;   _State[14] = 0.9942;
      _State[15] = 0.999998;  _State[16] = 2.347e-8; _State[17] = 4.272;
      _State[18] = 0.8978;
      */

      // Modified State Variables as requested by KCL
      _State[0] = -85.423; // Resting Voltage
      _State[1] = 136.89;  // Intracellular_K
      _State[2] = 8.604;   // Intracellular_Ni
      _State[3] = 0.000126;// Intracellular Calcium
      _State[4] = 0.00621; // Xr1
      _State[5] = 0.4712;  // Xr2
      _State[6] = 0.0095;  // Xs
      _State[7] = 0.00172; // m
      _State[8] = 0.7444;  // h
      _State[9] = 0.7045;  // j
      _State[10] = 0.00036; // Subspace Calcium
      _State[11] = 3.373e-5; // L_Type d
      _State[12] = 0.7888;   // L_Type f
      _State[13] = 0.9755;   // L_type f2
      _State[14] = 0.9953;   // L_Type fCass
      _State[15] = 0.999998; // Transient outward current s gate
      _State[16] = 2.42e-8;  // Transient outwatf current r gate
      _State[17] = 3.64;    // Sarcoplasmic Reticulum Calcium
      _State[18] = 0.9073;   // R_prime
      
      //! Surface Area to Volume Ratio from the paper mentioned above
      _xi        = 2000.;
    }

    //! Destructor
    ~Tusscher() {;}

    //! Compute Rates
    void ComputeRates();

    //! Compute Ionic Current
    Real compute(Real Xi, Real C_m, Real dt, Real volt,
		     Real istim);

    //! Get gamma
    Real getGamma();

    // get internal variables
    const std::vector<Real>& getInternalParameters(int& nData) const { 
      nData = 19; return _State; }

    // Set internal variables
    void setInternalParameters(const std::vector<Real>& data) {
      for(int i = 0; i < 19; i++) _State[i] = data[i];
    }      
  };
}
#endif
/*
 * VOI is time in component environment (millisecond).
 * STATES[0] is V in component membrane (millivolt).
 * CONSTANTS[0] is R in component membrane (joule_per_mole_kelvin).
 * CONSTANTS[1] is T in component membrane (kelvin).
 * CONSTANTS[2] is F in component membrane (coulomb_per_millimole).
 * CONSTANTS[3] is Cm in component membrane (microF).
 * CONSTANTS[4] is V_c in component membrane (micrometre3).
 * CONSTANTS[5] is stim_start in component membrane (millisecond).
 * CONSTANTS[6] is stim_period in component membrane (millisecond).
 * CONSTANTS[7] is stim_duration in component membrane (millisecond).
 * CONSTANTS[8] is stim_amplitude in component membrane (picoA_per_picoF).
 * CONSTANTS[9] is P_kna in component reversal_potentials (dimensionless).
 * CONSTANTS[10] is K_o in component potassium_dynamics (millimolar).
 * CONSTANTS[11] is Na_o in component sodium_dynamics (millimolar).
 * CONSTANTS[12] is Ca_o in component calcium_dynamics (millimolar).
 * CONSTANTS[13] is g_K1 in component inward_rectifier_potassium_current (nanoS_per_picoF).
 * CONSTANTS[14] is g_Kr in component rapid_time_dependent_potassium_current (nanoS_per_picoF).
 * CONSTANTS[15] is g_Ks in component slow_time_dependent_potassium_current (nanoS_per_picoF).
 * CONSTANTS[16] is g_Na in component fast_sodium_current (nanoS_per_picoF).
 * CONSTANTS[17] is g_bna in component sodium_background_current (nanoS_per_picoF).
 * CONSTANTS[18] is g_CaL in component L_type_Ca_current (litre_per_farad_second).
 * CONSTANTS[19] is g_bca in component calcium_background_current (nanoS_per_picoF).
 * CONSTANTS[20] is g_to in component transient_outward_current (nanoS_per_picoF).
 * CONSTANTS[21] is P_NaK in component sodium_potassium_pump_current (picoA_per_picoF).
 * CONSTANTS[22] is K_mk in component sodium_potassium_pump_current (millimolar).
 * CONSTANTS[23] is K_mNa in component sodium_potassium_pump_current (millimolar).
 * CONSTANTS[24] is K_NaCa in component sodium_calcium_exchanger_current (picoA_per_picoF).
 * CONSTANTS[25] is K_sat in component sodium_calcium_exchanger_current (dimensionless).
 * CONSTANTS[26] is alpha in component sodium_calcium_exchanger_current (dimensionless).
 * CONSTANTS[27] is gamma in component sodium_calcium_exchanger_current (dimensionless).
 * CONSTANTS[28] is Km_Ca in component sodium_calcium_exchanger_current (millimolar).
 * CONSTANTS[29] is Km_Nai in component sodium_calcium_exchanger_current (millimolar).
 * CONSTANTS[30] is g_pCa in component calcium_pump_current (picoA_per_picoF).
 * CONSTANTS[31] is K_pCa in component calcium_pump_current (millimolar).
 * CONSTANTS[32] is g_pK in component potassium_pump_current (nanoS_per_picoF).
 * CONSTANTS[33] is k1_prime in component calcium_dynamics (per_millimolar2_per_millisecond).
 * CONSTANTS[34] is k2_prime in component calcium_dynamics (per_millimolar_per_millisecond).
 * CONSTANTS[35] is k3 in component calcium_dynamics (per_millisecond).
 * CONSTANTS[36] is k4 in component calcium_dynamics (per_millisecond).
 * CONSTANTS[37] is EC in component calcium_dynamics (millimolar).
 * CONSTANTS[38] is max_sr in component calcium_dynamics (dimensionless).
 * CONSTANTS[39] is min_sr in component calcium_dynamics (dimensionless).
 * CONSTANTS[40] is V_rel in component calcium_dynamics (per_millisecond).
 * CONSTANTS[41] is V_xfer in component calcium_dynamics (per_millisecond).
 * CONSTANTS[42] is K_up in component calcium_dynamics (millimolar).
 * CONSTANTS[43] is V_leak in component calcium_dynamics (per_millisecond).
 * CONSTANTS[44] is Vmax_up in component calcium_dynamics (millimolar_per_millisecond).
 * CONSTANTS[45] is Buf_c in component calcium_dynamics (millimolar).
 * CONSTANTS[46] is K_buf_c in component calcium_dynamics (millimolar).
 * CONSTANTS[47] is Buf_sr in component calcium_dynamics (millimolar).
 * CONSTANTS[48] is K_buf_sr in component calcium_dynamics (millimolar).
 * CONSTANTS[49] is Buf_ss in component calcium_dynamics (millimolar).
 * CONSTANTS[50] is K_buf_ss in component calcium_dynamics (millimolar).
 * CONSTANTS[51] is V_sr in component calcium_dynamics (micrometre3).
 * CONSTANTS[52] is V_ss in component calcium_dynamics (micrometre3).

 * ALGEBRAIC[47] is i_K1 in component inward_rectifier_potassium_current (picoA_per_picoF).
 * ALGEBRAIC[54] is i_to in component transient_outward_current (picoA_per_picoF).
 * ALGEBRAIC[48] is i_Kr in component rapid_time_dependent_potassium_current (picoA_per_picoF).
 * ALGEBRAIC[49] is i_Ks in component slow_time_dependent_potassium_current (picoA_per_picoF).
 * ALGEBRAIC[52] is i_CaL in component L_type_Ca_current (picoA_per_picoF).
 * ALGEBRAIC[55] is i_NaK in component sodium_potassium_pump_current (picoA_per_picoF).
 * ALGEBRAIC[50] is i_Na in component fast_sodium_current (picoA_per_picoF).
 * ALGEBRAIC[51] is i_b_Na in component sodium_background_current (picoA_per_picoF).
 * ALGEBRAIC[56] is i_NaCa in component sodium_calcium_exchanger_current (picoA_per_picoF).
 * ALGEBRAIC[53] is i_b_Ca in component calcium_background_current (picoA_per_picoF).
 * ALGEBRAIC[58] is i_p_K in component potassium_pump_current (picoA_per_picoF).
 * ALGEBRAIC[57] is i_p_Ca in component calcium_pump_current (picoA_per_picoF).
 * ALGEBRAIC[12] is i_Stim in component membrane (picoA_per_picoF).

 * ALGEBRAIC[25] is E_Na in component reversal_potentials (millivolt).
 * ALGEBRAIC[33] is E_K in component reversal_potentials (millivolt).
 * ALGEBRAIC[41] is E_Ks in component reversal_potentials (millivolt).
 * ALGEBRAIC[43] is E_Ca in component reversal_potentials (millivolt).


 * ALGEBRAIC[46] is xK1_inf in component inward_rectifier_potassium_current (dimensionless).
 * ALGEBRAIC[44] is alpha_K1 in component inward_rectifier_potassium_current (dimensionless).
 * ALGEBRAIC[45] is beta_K1 in component inward_rectifier_potassium_current (dimensionless).
 * ALGEBRAIC[0] is xr1_inf in component rapid_time_dependent_potassium_current_Xr1_gate (dimensionless).
 * ALGEBRAIC[13] is alpha_xr1 in component rapid_time_dependent_potassium_current_Xr1_gate (dimensionless).
 * ALGEBRAIC[26] is beta_xr1 in component rapid_time_dependent_potassium_current_Xr1_gate (dimensionless).
 * ALGEBRAIC[34] is tau_xr1 in component rapid_time_dependent_potassium_current_Xr1_gate (millisecond).
 * ALGEBRAIC[1] is xr2_inf in component rapid_time_dependent_potassium_current_Xr2_gate (dimensionless).
 * ALGEBRAIC[14] is alpha_xr2 in component rapid_time_dependent_potassium_current_Xr2_gate (dimensionless).
 * ALGEBRAIC[27] is beta_xr2 in component rapid_time_dependent_potassium_current_Xr2_gate (dimensionless).
 * ALGEBRAIC[35] is tau_xr2 in component rapid_time_dependent_potassium_current_Xr2_gate (millisecond).
 * ALGEBRAIC[2] is xs_inf in component slow_time_dependent_potassium_current_Xs_gate (dimensionless).
 * ALGEBRAIC[15] is alpha_xs in component slow_time_dependent_potassium_current_Xs_gate (dimensionless).
 * ALGEBRAIC[28] is beta_xs in component slow_time_dependent_potassium_current_Xs_gate (dimensionless).
 * ALGEBRAIC[36] is tau_xs in component slow_time_dependent_potassium_current_Xs_gate (millisecond).
 * ALGEBRAIC[3] is m_inf in component fast_sodium_current_m_gate (dimensionless).
 * ALGEBRAIC[16] is alpha_m in component fast_sodium_current_m_gate (dimensionless).
 * ALGEBRAIC[29] is beta_m in component fast_sodium_current_m_gate (dimensionless).
 * ALGEBRAIC[37] is tau_m in component fast_sodium_current_m_gate (millisecond).
 * ALGEBRAIC[4] is h_inf in component fast_sodium_current_h_gate (dimensionless).
 * ALGEBRAIC[17] is alpha_h in component fast_sodium_current_h_gate (per_millisecond).
 * ALGEBRAIC[30] is beta_h in component fast_sodium_current_h_gate (per_millisecond).
 * ALGEBRAIC[38] is tau_h in component fast_sodium_current_h_gate (millisecond).
 * ALGEBRAIC[5] is j_inf in component fast_sodium_current_j_gate (dimensionless).
 * ALGEBRAIC[18] is alpha_j in component fast_sodium_current_j_gate (per_millisecond).
 * ALGEBRAIC[31] is beta_j in component fast_sodium_current_j_gate (per_millisecond).
 * ALGEBRAIC[39] is tau_j in component fast_sodium_current_j_gate (millisecond).

 * STATES[1] is K_i in component potassium_dynamics (millimolar).
 * STATES[2] is Na_i in component sodium_dynamics (millimolar).
 * STATES[3] is Ca_i in component calcium_dynamics (millimolar).
 * STATES[4] is Xr1 in component rapid_time_dependent_potassium_current_Xr1_gate (dimensionless).
 * STATES[5] is Xr2 in component rapid_time_dependent_potassium_current_Xr2_gate (dimensionless).
 * STATES[6] is Xs in component slow_time_dependent_potassium_current_Xs_gate (dimensionless).
 * STATES[7] is m in component fast_sodium_current_m_gate (dimensionless).
 * STATES[8] is h in component fast_sodium_current_h_gate (dimensionless).
 * STATES[9] is j in component fast_sodium_current_j_gate (dimensionless).
 * STATES[10] is Ca_ss in component calcium_dynamics (millimolar).
 * STATES[11] is d in component L_type_Ca_current_d_gate (dimensionless).
 * STATES[12] is f in component L_type_Ca_current_f_gate (dimensionless).
 * STATES[13] is f2 in component L_type_Ca_current_f2_gate (dimensionless).
 * STATES[14] is fCass in component L_type_Ca_current_fCass_gate (dimensionless).
 * STATES[15] is s in component transient_outward_current_s_gate (dimensionless).
 * STATES[16] is r in component transient_outward_current_r_gate (dimensionless).
 * STATES[17] is Ca_SR in component calcium_dynamics (millimolar).


 * ALGEBRAIC[6] is d_inf in component L_type_Ca_current_d_gate (dimensionless).
 * ALGEBRAIC[19] is alpha_d in component L_type_Ca_current_d_gate (dimensionless).
 * ALGEBRAIC[32] is beta_d in component L_type_Ca_current_d_gate (dimensionless).
 * ALGEBRAIC[40] is gamma_d in component L_type_Ca_current_d_gate (millisecond).
 * ALGEBRAIC[42] is tau_d in component L_type_Ca_current_d_gate (millisecond).
 * ALGEBRAIC[7] is f_inf in component L_type_Ca_current_f_gate (dimensionless).
 * ALGEBRAIC[20] is tau_f in component L_type_Ca_current_f_gate (millisecond).
 * ALGEBRAIC[8] is f2_inf in component L_type_Ca_current_f2_gate (dimensionless).
 * ALGEBRAIC[21] is tau_f2 in component L_type_Ca_current_f2_gate (millisecond).
 * ALGEBRAIC[9] is fCass_inf in component L_type_Ca_current_fCass_gate (dimensionless).
 * ALGEBRAIC[22] is tau_fCass in component L_type_Ca_current_fCass_gate (millisecond).
 * ALGEBRAIC[10] is s_inf in component transient_outward_current_s_gate (dimensionless).
 * ALGEBRAIC[23] is tau_s in component transient_outward_current_s_gate (millisecond).
 * ALGEBRAIC[11] is r_inf in component transient_outward_current_r_gate (dimensionless).
 * ALGEBRAIC[24] is tau_r in component transient_outward_current_r_gate (millisecond).

 * ALGEBRAIC[67] is i_rel in component calcium_dynamics (millimolar_per_millisecond).
 * ALGEBRAIC[59] is i_up in component calcium_dynamics (millimolar_per_millisecond).
 * ALGEBRAIC[60] is i_leak in component calcium_dynamics (millimolar_per_millisecond).
 * ALGEBRAIC[61] is i_xfer in component calcium_dynamics (millimolar_per_millisecond).
 * ALGEBRAIC[62] is kcasr in component calcium_dynamics (dimensionless).
 * ALGEBRAIC[66] is O in component calcium_dynamics (dimensionless).
 * STATES[18] is R_prime in component calcium_dynamics (dimensionless).
 * ALGEBRAIC[64] is k1 in component calcium_dynamics (per_millimolar2_per_millisecond).
 * ALGEBRAIC[65] is k2 in component calcium_dynamics (per_millimolar_per_millisecond).
 * ALGEBRAIC[63] is Ca_i_bufc in component calcium_dynamics (dimensionless).
 * ALGEBRAIC[68] is Ca_sr_bufsr in component calcium_dynamics (dimensionless).
 * ALGEBRAIC[69] is Ca_ss_bufss in component calcium_dynamics (dimensionless).
 * RATES[0] is d/dt V in component membrane (millivolt).
 * RATES[4] is d/dt Xr1 in component rapid_time_dependent_potassium_current_Xr1_gate (dimensionless).
 * RATES[5] is d/dt Xr2 in component rapid_time_dependent_potassium_current_Xr2_gate (dimensionless).
 * RATES[6] is d/dt Xs in component slow_time_dependent_potassium_current_Xs_gate (dimensionless).
 * RATES[7] is d/dt m in component fast_sodium_current_m_gate (dimensionless).
 * RATES[8] is d/dt h in component fast_sodium_current_h_gate (dimensionless).
 * RATES[9] is d/dt j in component fast_sodium_current_j_gate (dimensionless).
 * RATES[11] is d/dt d in component L_type_Ca_current_d_gate (dimensionless).
 * RATES[12] is d/dt f in component L_type_Ca_current_f_gate (dimensionless).
 * RATES[13] is d/dt f2 in component L_type_Ca_current_f2_gate (dimensionless).
 * RATES[14] is d/dt fCass in component L_type_Ca_current_fCass_gate (dimensionless).
 * RATES[15] is d/dt s in component transient_outward_current_s_gate (dimensionless).
 * RATES[16] is d/dt r in component transient_outward_current_r_gate (dimensionless).
 * RATES[18] is d/dt R_prime in component calcium_dynamics (dimensionless).
 * RATES[3] is d/dt Ca_i in component calcium_dynamics (millimolar).
 * RATES[17] is d/dt Ca_SR in component calcium_dynamics (millimolar).
 * RATES[10] is d/dt Ca_ss in component calcium_dynamics (millimolar).
 * RATES[2] is d/dt Na_i in component sodium_dynamics (millimolar).
 * RATES[1] is d/dt K_i in component potassium_dynamics (millimolar)
 */

