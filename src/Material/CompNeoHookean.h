/*! 
  \file CompNeoHookean.h
  \brief Interface for a compressible Neo-Hookean finite deformation
  hyperelasticity model.
*/

#ifndef _COMPNEOHOOKEAN_H_
#define _COMPNEOHOOKEAN_H_

#include "MechanicsMaterial.h"

namespace voom {
  
  class CompNeoHookean : public MechanicsMaterial
  {
  public: 
    // Constructors/destructors:
    CompNeoHookean(int ID): MechanicsMaterial(ID), _lambda(1.0), _mu(1.0) {};
    CompNeoHookean(int ID, Real Lambda, Real Mu): MechanicsMaterial(ID), _lambda(Lambda), _mu(Mu) {};
    CompNeoHookean(CompNeoHookean* BaseMaterial): 
    MechanicsMaterial(BaseMaterial->_matID), _lambda(BaseMaterial->_lambda), _mu(BaseMaterial->_mu) {};

    // Clone
    virtual CompNeoHookean* clone() const {
      return new CompNeoHookean(*this);
    }
      
    // Default copy constructor (compiler should already provide exactly this)
    CompNeoHookean(const CompNeoHookean & Old): 
    MechanicsMaterial(Old._matID), _lambda(Old._lambda), _mu(Old._mu) {};

    Real getLambda() {return _lambda;};
    Real getMu() {return _mu;};
    void setLambda(Real Lambda) {_lambda = Lambda;};
    void setMu(Real Mu) {_mu = Mu;};
    void setLameConstFromEnu(Real E, Real nu) {
      _lambda = nu*E/((1.0+nu)*(1.0-2.0*nu)); 
      _mu = 0.5*E/(1.0+nu);
    }

    void setMaterialParameters(const vector<Real > & LambdaMu) {
      _lambda = LambdaMu[0];
      _mu = LambdaMu[1];
    }
    void setInternalParameters(const vector<Real > & IntParam) {}; // No internal parameters for CompNeoHookean
    void setRegularizationParameters(const vector<Real > &)    {}; // CompNeoHookean does not have regularization parameters

    vector<Real > getMaterialParameters() {
      vector<Real > MatProp(2, 0.0);
      MatProp[0] = _lambda;
      MatProp[1] = _mu;
      return MatProp;
    }
    vector<Real > getInternalParameters() { // No internal parameters for CompNeoHookean
      vector<Real > IntParam;
      return IntParam;
    }
    vector<Real > getRegularizationParameters() { // No regularization parameters for CompNeoHookean
      vector<Real > RegParam;
      return RegParam;
    }

    
    // Operators
    //! Based on deformation gradient tensor F, calculates state of material
    void compute(FKresults & R, const Matrix3d & F, Vector3d * Fiber = NULL);

    //! Tells if material has history variables and needs to be duplicated at each quadrature point
    // It is used in the Model derived classes
    bool HasHistoryVariables() { return false; };
    
  private:
    //! Data
    Real _lambda;
    Real _mu;
	
  }; // class CompNeoHookean 
  
}
#endif // _COMPNEOHOOKEAN_H_
