//-*-C++-*-
/*! \brief
  Ionic model for Rabbit Purkinje structure. Developed based on the paper by
  Corrias et al titled "Ionic mechanisms of electrophysiological properties 
  and repolarization abnormalities in rabbit Purkinje fibers", American
  Journal of HeartCirculation Physiology, Feb 2011
  - Added Faster computation as suggested by Whiteley. We do Explicit 
  update on Linear ODE and Newton-Raphson on Non-Linear ODE's.

*/
#ifndef _Purkinje_h_
#define _Purkinje_h_

#include "IonicMaterial.h"

namespace voom {
  class Purkinje: public IonicMaterial {
  private:
    //! Constants for the model
    Real  *_Constants;

    //! State variables
    std::vector<Real>  _State;

    //! Rate Variables
    Real  _Rates[22];

    //! Algebraic variable array
    Real  _Algebraic[69];

    //! Surface Area of the Cell
    Real  _SurfaceArea;

    //! Forward Euler Adaptive Update of State Variables
    void UpdateStateVariables(const Real dt);

    //!
  public:
    //! Constructor
    Purkinje(Real *constants){
      _Constants = constants;
      _State.resize( 22 );
      // Initialize States
      _State[0] = -88.34;      _State[1] = 0.00001;      _State[2] = 0.000032;
      _State[3] = 0.17;        _State[4] = 6.7;          _State[5] = 140;    
      _State[6] = 0.001337;    _State[7] = 0.01;         _State[8] = 0.000003;
      _State[9] = 0.1;         _State[10] = 0.7;         _State[11] = 0.0;
      _State[12] = 0.7;        _State[13] = 0.000007;    _State[14] = 0.978861;
      _State[15] = 0.000012;   _State[16] = 0.864489;    _State[17] = 0.25;
      _State[18] = 1.0;        _State[19] = 0.0;         _State[20] = 0.011099;
      _State[21] = 0.011099;
      // Surface Area to Volume ratio from the paper cm^-1
      _xi        =  2500.;
      // Surface Area of the Cell from the paper 4612 um^2. Converted to cm^2
      _SurfaceArea = 4.612E-5 ;// cm^2
    }

    //! Destructor
    ~Purkinje() {;}
    
    //! Compute Rates
    void ComputeRates();    

    //! Compute Variables
    void ComputeVariables();

    //! Compute Ionic Current
    Real compute(Real Xi, Real C_m, Real dt, Real volt,
		     Real istim);

    //! Get gamma
    Real getGamma();

    // get internal variables
    const std::vector<Real>& getInternalParameters(int &nData) const { 
      nData = 22; return _State; }

    // Set internal variables
    void setInternalParameters(const std::vector<Real>& data) {
      for(int i = 0; i < 22; i++) _State[i] = data[i];
    }
  };  
}

/*
  DEFINITION OF PRIVATE VARIABLES IN THE MODEL
  ============================================
  There are a total of 69 entries in the algebraic variable array.
  There are a total of 22 entries in each of the rate and state variable 
  arrays. There are a total of 62 entries in the constant variable array.
  
  * VOI is time in component Time (time_units).
  * CONSTANTS[0] is T in component Environment (Temperature_units).
  * CONSTANTS[1] is F in component Environment (F_units).
  * CONSTANTS[2] is R in component Environment (R_units).
  * CONSTANTS[3] is Ca_o in component Environment (millimolar).
  * CONSTANTS[4] is Na_o in component Environment (millimolar).
  * CONSTANTS[5] is K_o in component Environment (millimolar).
  * CONSTANTS[6] is Cm in component membrane (capacitance_units).
  * CONSTANTS[7] is cell_volume in component membrane (volume_units).
  * CONSTANTS[8] is bulk_fraction in component membrane (dimensionless).
  * CONSTANTS[9] is periphery_fraction in component membrane (dimensionless).
  * CONSTANTS[10] is sr_fraction in component membrane (dimensionless).
  * CONSTANTS[11] is stim_start in component membrane (time_units).
  * CONSTANTS[12] is stim_end in component membrane (time_units).
  * CONSTANTS[13] is stim_period in component membrane (time_units).
  * CONSTANTS[14] is stim_duration in component membrane (time_units).
  * CONSTANTS[15] is stim_amplitude in component membrane (current_units).
  * CONSTANTS[16] is tau_x_Ttype in component x_Ttype (time_units).
  * CONSTANTS[17] is G_max_Ttype in component i_cat (conductance_units).
  * CONSTANTS[18] is tau_x_Ltype in component x_Ltype (time_units).
  * CONSTANTS[19] is G_max_Ltype in component i_cal (conductance_units).
  * CONSTANTS[20] is tau_x_to_fast in component x_to_fast (time_units).
  * CONSTANTS[21] is tau_y_to_fast in component y_to_fast (time_units).
  * CONSTANTS[22] is G_max_to_fast in component i_to_fast (conductance_units).
  * CONSTANTS[23] is G_max_to_sustained in component i_to_sustained 
  (conductance_units).
  * CONSTANTS[24] is tau_x_na_fast in component x_na_fast (time_units).
  * CONSTANTS[25] is tau_y_na_fast in component y_na_fast (time_units).
  * CONSTANTS[26] is G_max_na_fast in component i_na_fast (conductance_units).
  * CONSTANTS[27] is tau_x_na_late in component x_na_late (time_units).
  * CONSTANTS[28] is G_max_na_late in component i_na_late (conductance_units).
  * CONSTANTS[29] is G_max_k1 in component i_k1 (conductance_units).
  * CONSTANTS[30] is G_max_kr in component i_kr (conductance_units).
  * CONSTANTS[31] is G_max_ks in component i_ks (conductance_units).
  * CONSTANTS[32] is n_NaCa in component i_naca (dimensionless).
  * CONSTANTS[33] is g_NaCa in component i_naca (current_units).
  * CONSTANTS[34] is d_NaCa in component i_naca (dimensionless).
  * CONSTANTS[35] is gamma in component i_naca (dimensionless).
  * CONSTANTS[36] is g_nak in component i_nak (current_units).
  * CONSTANTS[37] is G_max_kb in component i_kb (conductance_units).
  * CONSTANTS[38] is G_max_nab in component i_nab (conductance_units).
  * CONSTANTS[39] is G_max_cab in component i_cab (conductance_units).
  * CONSTANTS[40] is PMCA_max in component i_pmca (current_units).
  * CONSTANTS[41] is Kpmca in component i_pmca (millimolar).
  * CONSTANTS[42] is Hpmca in component i_pmca (dimensionless).
  * CONSTANTS[43] is G_f_k in component i_f_k (conductance_units).
  * CONSTANTS[44] is G_f_na in component i_f_na (conductance_units).
  * CONSTANTS[45] is REL_max in component i_rel_per 
  (current_per_millimolar_units).
  * CONSTANTS[46] is Krel in component i_rel_per (millimolar).
  * CONSTANTS[47] is SERCA_max in component i_serca_per (current_units).
  * CONSTANTS[48] is Kmf in component i_serca_per (millimolar).
  * CONSTANTS[49] is Kmr in component i_serca_per (millimolar).
  * CONSTANTS[50] is H in component i_serca_per (dimensionless).
  * CONSTANTS[51] is Kmf in component i_serca_bulk (millimolar).
  * CONSTANTS[52] is Kmr in component i_serca_bulk (millimolar).
  * CONSTANTS[53] is H in component i_serca_bulk (dimensionless).
  * CONSTANTS[54] is LEAK_max in component i_leak(current_per_millimolar_units)
  * CONSTANTS[55] is DIFF_max in component i_diff(current_per_millimolar_units)
  * CONSTANTS[56] is FoRT in component Environment (inverse_voltage_units).
  * CONSTANTS[57] is RToF in component Environment (voltage_units).
  * CONSTANTS[58] is peripheral_volume in component membrane (volume_units).
  * CONSTANTS[59] is bulk_volume in component membrane (volume_units).
  * CONSTANTS[60] is diffusable_volume in component membrane (volume_units).
  * CONSTANTS[61] is sr_volume in component membrane (volume_units).

  * STATES[0] is Vm in component membrane (voltage_units).
  * STATES[1] is Ca_i_peripheral in component membrane (millimolar).
  * STATES[2] is Ca_i_bulk in component membrane (millimolar).
  * STATES[3] is Ca_sr in component membrane (millimolar).
  * STATES[4] is Na_i in component membrane (millimolar).
  * STATES[5] is K_i in component membrane (millimolar).
  * STATES[6] is x_Ttype in component x_Ttype (dimensionless).
  * STATES[7] is y_Ttype in component y_Ttype (dimensionless).
  * STATES[8] is x_Ltype in component x_Ltype (dimensionless).
  * STATES[9] is y_Ltype in component y_Ltype (dimensionless).
  * STATES[10] is y_ca_Ltype in component y_ca_Ltype (dimensionless).
  * STATES[11] is x_to_fast in component x_to_fast (dimensionless).
  * STATES[12] is y_to_fast in component y_to_fast (dimensionless).
  * STATES[13] is x_na_fast in component x_na_fast (dimensionless).
  * STATES[14] is y_na_fast in component y_na_fast (dimensionless).
  * STATES[15] is x_na_late in component x_na_late (dimensionless).
  * STATES[16] is y_na_late in component y_na_late (dimensionless).
  * STATES[17] is y_kr in component y_kr (dimensionless).
  * STATES[18] is x_ks in component x_ks (dimensionless).
  * STATES[19] is y_ks in component y_ks (dimensionless).
  * STATES[20] is y_gate_f_k in component y_gate_f_k (dimensionless).
  * STATES[21] is y_gate_f_na in component y_gate_f_na (dimensionless).

  * ALGEBRAIC[0] is i_stim in component membrane (current_units).
  * ALGEBRAIC[1] is x_inf_Ttype in component x_Ttype (dimensionless).
  * ALGEBRAIC[2] is y_inf_Ttype in component y_Ttype (dimensionless).
  * ALGEBRAIC[3] is x_inf_Ltype in component x_Ltype (dimensionless).
  * ALGEBRAIC[4] is y_inf_Ltype in component y_Ltype (dimensionless).
  * ALGEBRAIC[5] is y_ca_inf_Ltype in component y_ca_Ltype (dimensionless).
  * ALGEBRAIC[6] is x_inf_to_fast in component x_to_fast (dimensionless)
  * ALGEBRAIC[7] is y_inf_to_fast in component y_to_fast (dimensionless).
  * ALGEBRAIC[8] is x_inf_na_fast in component x_na_fast (dimensionless).
  * ALGEBRAIC[9] is y_inf_na_fast in component y_na_fast (dimensionless).
  * ALGEBRAIC[10] is x_inf_na_late in component x_na_late (dimensionless).
  * ALGEBRAIC[11] is y_inf_na_late in component y_na_late (dimensionless).
  * ALGEBRAIC[12] is y_inf_kr in component y_kr (dimensionless).
  * ALGEBRAIC[13] is x_inf_ks in component x_ks (dimensionless).
  * ALGEBRAIC[14] is y_inf_f_gate in component y_gate_f_k (dimensionless).
  * ALGEBRAIC[15] is y_inf_f_gate in component y_gate_f_na (dimensionless).
  * ALGEBRAIC[16] is E_Ca in component i_cat (voltage_units).
  * ALGEBRAIC[17] is tau_y_Ttype in component y_Ttype (time_units).
  * ALGEBRAIC[18] is tau_y_Ltype in component y_Ltype (time_units).
  * ALGEBRAIC[19] is tau_y_ca_Ltype in component y_ca_Ltype (time_units).
  * ALGEBRAIC[20] is tau_y_na_late in component y_na_late (time_units).
  * ALGEBRAIC[21] is ykrv1 in component y_kr (rate_constants_units).
  * ALGEBRAIC[22] is tau_x_ks in component x_ks (time_units).
  * ALGEBRAIC[23] is tau_y_f_gate in component y_gate_f_k (time_units).
  * ALGEBRAIC[24] is tau_y_f_gate in component y_gate_f_na (time_units).
  * ALGEBRAIC[25] is i_cat in component i_cat (current_units).
  * ALGEBRAIC[26] is ykrv2 in component y_kr (rate_constants_units).
  * ALGEBRAIC[27] is y_inf_ks in component y_ks (dimensionless).
  * ALGEBRAIC[28] is E_Ca in component i_cal (voltage_units).
  * ALGEBRAIC[29] is tau_y_kr in component y_kr (time_units).
  * ALGEBRAIC[30] is tau_y_ks in component y_ks (time_units).
  * ALGEBRAIC[31] is i_cal in component i_cal (current_units).
  * ALGEBRAIC[32] is E_k in component i_to_fast (voltage_units).
  * ALGEBRAIC[33] is i_to_fast in component i_to_fast (current_units).
  * ALGEBRAIC[34] is x_to_sustained in component x_to_sustained (dimensionless)
  * ALGEBRAIC[35] is E_k in component i_to_sustained (voltage_units).
  * ALGEBRAIC[36] is i_to_sustained in component i_to_sustained (current_units)
  * ALGEBRAIC[37] is E_na in component i_na_fast (voltage_units).
  * ALGEBRAIC[38] is i_na_fast in component i_na_fast (current_units).
  * ALGEBRAIC[39] is E_na in component i_na_late (voltage_units).
  * ALGEBRAIC[40] is i_na_late in component i_na_late (current_units).
  * ALGEBRAIC[41] is x_k1 in component x_k1 (dimensionless).
  * ALGEBRAIC[42] is E_k in component i_k1 (voltage_units).
  * ALGEBRAIC[43] is i_k1 in component i_k1 (current_units).
  * ALGEBRAIC[44] is x_kr in component x_kr (dimensionless).
  * ALGEBRAIC[45] is E_k in component i_kr (voltage_units).
  * ALGEBRAIC[46] is i_kr in component i_kr (current_units).
  * ALGEBRAIC[47] is E_k in component i_ks (voltage_units).
  * ALGEBRAIC[48] is i_ks in component i_ks (current_units).
  * ALGEBRAIC[49] is i_naca in component i_naca (current_units).
  * ALGEBRAIC[50] is x_nak in component x_nak (dimensionless).
  * ALGEBRAIC[51] is y_nak in component y_nak (dimensionless).
  * ALGEBRAIC[52] is i_nak in component i_nak (current_units).
  * ALGEBRAIC[53] is E_k in component i_kb (voltage_units).
  * ALGEBRAIC[54] is i_kb in component i_kb (current_units).
  * ALGEBRAIC[55] is E_na in component i_nab (voltage_units).
  * ALGEBRAIC[56] is i_nab in component i_nab (current_units).
  * ALGEBRAIC[57] is E_ca in component i_cab (voltage_units).
  * ALGEBRAIC[58] is i_cab in component i_cab (current_units).
  * ALGEBRAIC[59] is i_pmca in component i_pmca (current_units).
  * ALGEBRAIC[60] is E_k in component i_f_k (voltage_units).
  * ALGEBRAIC[61] is i_rel_per in component i_rel_per (current_units).
  * ALGEBRAIC[62] is i_f_k in component i_f_k (current_units).
  * ALGEBRAIC[63] is i_serca_per in component i_serca_per (current_units).
  * ALGEBRAIC[64] is E_na in component i_f_na (voltage_units).
  * ALGEBRAIC[65] is i_serca_bulk in component i_serca_bulk (current_units).
  * ALGEBRAIC[66] is i_f_na in component i_f_na (current_units).
  * ALGEBRAIC[67] is i_leak in component i_leak (current_units).
  * ALGEBRAIC[68] is i_diff in component i_diff (current_units).
  * RATES[0] is d/dt Vm in component membrane (voltage_units).
  * RATES[1] is d/dt Ca_i_peripheral in component membrane (millimolar).
  * RATES[2] is d/dt Ca_i_bulk in component membrane (millimolar).
  * RATES[3] is d/dt Ca_sr in component membrane (millimolar).
  * RATES[5] is d/dt K_i in component membrane (millimolar).
  * RATES[4] is d/dt Na_i in component membrane (millimolar).
  * RATES[6] is d/dt x_Ttype in component x_Ttype (dimensionless).
  * RATES[7] is d/dt y_Ttype in component y_Ttype (dimensionless).
  * RATES[8] is d/dt x_Ltype in component x_Ltype (dimensionless).
  * RATES[9] is d/dt y_Ltype in component y_Ltype (dimensionless).
  * RATES[10] is d/dt y_ca_Ltype in component y_ca_Ltype (dimensionless).
  * RATES[11] is d/dt x_to_fast in component x_to_fast (dimensionless).
  * RATES[12] is d/dt y_to_fast in component y_to_fast (dimensionless).
  * RATES[13] is d/dt x_na_fast in component x_na_fast (dimensionless).
  * RATES[14] is d/dt y_na_fast in component y_na_fast (dimensionless).
  * RATES[15] is d/dt x_na_late in component x_na_late (dimensionless).
  * RATES[16] is d/dt y_na_late in component y_na_late (dimensionless).
  * RATES[17] is d/dt y_kr in component y_kr (dimensionless).
  * RATES[18] is d/dt x_ks in component x_ks (dimensionless).
  * RATES[19] is d/dt y_ks in component y_ks (dimensionless).
  * RATES[20] is d/dt y_gate_f_k in component y_gate_f_k (dimensionless).
  * RATES[21] is d/dt y_gate_f_na in component y_gate_f_na (dimensionless).
  */

#endif

