/*! 
  \file PassMyoA.h
  \brief Interface for a passive myocardium finite deformation
  hyperelasticity model.
*/

#ifndef _PASSMYOA_H_
#define _PASSMYOA_H_

#include "MechanicsMaterial.h"

namespace voom {
  
  class PassMyoA : public MechanicsMaterial
  {
  public: 
    // Constructors/destructors:
  PassMyoA(int ID):  MechanicsMaterial(ID), _alpha1(1.0), _alpha2(1.0), _beta(1.0), _gamma(1.0), _a1(1.0), _a2(1.0) {};
  PassMyoA(int ID, Real Alpha1, Real Alpha2, Real Beta, Real Gamma, Real A1, Real A2, vector<Vector3d > Fibers): MechanicsMaterial(ID), _alpha1(Alpha1), _alpha2(Alpha2), _beta(Beta), _gamma(Gamma), _a1(A1), _a2(A2), _fibers(Fibers) {};
  PassMyoA(PassMyoA* BaseMaterial): 
    MechanicsMaterial(BaseMaterial->_matID), _alpha1(BaseMaterial->_alpha1),  _alpha2(BaseMaterial->_alpha2), _beta(BaseMaterial->_beta), _gamma(BaseMaterial->_gamma), _a1(BaseMaterial->_a1), _a2(BaseMaterial->_a2), _fibers(BaseMaterial->_fibers) {};

    // Clone
    virtual PassMyoA* clone() const {
      return new PassMyoA(*this);
    }
      
    // Default copy constructor (compiler should already provide exactly this)
    PassMyoA(const PassMyoA & Old): 
    MechanicsMaterial(Old._matID), _alpha1(Old._alpha1), _alpha2(Old._alpha2), _beta(Old._beta), _gamma(Old._gamma), _a1(Old._a1), _a2(Old._a2), _fibers(Old._fibers) {};

    void setMaterialParameters(const vector<Real > & MatPar) {
      _alpha1 = MatPar[0];
      _alpha2 = MatPar[1];
    }
    void setInternalParameters(const vector<Real > & IntPar) {
      _a1 = IntPar[0];
      _a2 = IntPar[1];
    }
    void setRegularizationParameters(const vector<Real > & RegPar) { 
      _beta   = RegPar[0];
      _gamma  = RegPar[1];
    };

    vector<Real > getMaterialParameters() {
      vector<Real > MatProp(2, 0.0);
      MatProp[0] = _alpha1;
      MatProp[1] = _alpha2;
      return MatProp;
    }
    vector<Real > getInternalParameters() {
      vector<Real > IntPar(2, 0.0);
      IntPar[0] = _a1;
      IntPar[1] = _a2;
      return IntPar;
    }
    vector<Real > getRegularizationParameters() {
      vector<Real > RegProp(2, 0.0);
      RegProp[0] = _beta;
      RegProp[1] = _gamma;
      return RegProp;      
    };

    // Operators
    //! Based on deformation gradient tensor F, calculates state of material
    void compute(FKresults & R, const Matrix3d & F);

    //! Tells if material has history variables and needs to be duplicated at each quadrature point
    // It is used in the Model derived classes
    bool HasHistoryVariables() { return false; };
    
  private:
    //! Data
    Real _alpha1;
    Real _alpha2;
    Real _beta;
    Real _gamma;
    Real _a1;
    Real _a2;

    vector<Vector3d > _fibers;
	
  }; // class PassMyoA
  
}
#endif // _PASSMYOA_H_
