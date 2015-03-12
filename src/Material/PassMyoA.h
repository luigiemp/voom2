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
    PassMyoA(): _alpha(1.0), _beta(1.0), _a1(1.0), _a2(1.0) {};
    PassMyoA(Real Alpha, Real Beta, Real A1, Real A2): _alpha(Alpha), _beta(Beta), _a1(A1), _a2(A2) {};
    PassMyoA(PassMyoA* BaseMaterial): 
    _alpha(BaseMaterial->_alpha), _beta(BaseMaterial->_beta), _a1(BaseMaterial->_a1), _a2(BaseMaterial->_a2) {};

    // Clone
    virtual PassMyoA* clone() const {
      return new PassMyoA(*this);
    }
      
    // Default copy constructor (compiler should already provide exactly this)
    PassMyoA(const PassMyoA & Old): 
    _alpha(Old._alpha), _beta(Old._beta), _a1(Old._a1), _a2(Old._a2) {};

    void setMaterialParameters(const vector<Real > & MatPar) {
      _alpha = MatPar[0];
      _beta  = MatPar[1];
    }
    void setInternalParameters(const vector<Real > & IntPar) {
      _a1 = IntPar[0];
      _a2 = IntPar[1];
    }
    void setN(const Vector3d & N) {
      _N = N;
    }

    vector<Real > getMaterialParameters() {
      vector<Real > MatProp(2, 0.0);
      MatProp[0] = _alpha;
      MatProp[1] = _beta;
      return MatProp;
    }
    vector<Real > getInternalParameters() {
      vector<Real > IntPar(2, 0.0);
      IntPar[0] = _a1;
      IntPar[1] = _a2;
      return IntPar;
    }
    Vector3d getN() {
      return _N;
    }

    // Operators
    //! Based on deformation gradient tensor F, calculates state of material
    void compute(FKresults & R, const Matrix3d & F);

    //! Tells if material has history variables and needs to be duplicated at each quadrature point
    // It is used in the Model derived classes
    bool HasHistoryVariables() { return false; };
    
  private:
    //! Data
    Real _alpha;
    Real _beta;
    Real _a1;
    Real _a2;

    Vector3d _N;
	
  }; // class PassMyoA
  
}
#endif // _PASSMYOA_H_
