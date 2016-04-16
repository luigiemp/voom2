/*! 
  \file Jacobian.h
  \brief Interface to a material used to compute the ratio between the element reference and current volume
*/

#ifndef _ANGLESPRING_H_
#define _ANGLESPRING_H_

#include "FilamentMaterial.h"

namespace voom {
  
  class AngleSpring : public FilamentMaterial
  {
  public: 
    // Constructors/destructors:

  AngleSpring(int ID, Real kappa): FilamentMaterial(ID), _kappa(kappa) {}; 
  AngleSpring(AngleSpring* BaseMaterial): 
    FilamentMaterial(BaseMaterial->_matID), _kappa(BaseMaterial->_kappa) {};
    
    // Clone
    virtual AngleSpring* clone() const {
      return new AngleSpring(*this);
    }
    
    // Default copy constructor (compiler should already provide exactly this)
    AngleSpring(const AngleSpring & Old): 
    FilamentMaterial(Old._matID), _kappa(Old._kappa) {};

    void setMaterialParameters(const Real  & kappa) {
      _kappa = kappa;
    }; 
    void setInternalParameters(const Vector3d & ) {}; 
    void setRegularizationParameters(const vector<Real > &)    {}; // No regularization parameters for Spring

    Real getMaterialParameters() { 
      Real  MatProp;
      MatProp = _kappa;
      return MatProp;
    }
    Real getInternalParameters() {
      Real IntParam;
      return IntParam;
    }

    vector<Real> getRegularizationParameters() { // No regularization parameters for Spring
      vector<Real > RegParam;
      return RegParam;
    }

    
    // Operators
    //! Based on current length d, calculates state of the spring
    //    void compute(Filresults & R,  Vector3d nodeA, Vector3d nodeB, Vector3d nodeC);

    void compute(Filresults & R,  Vector3d & d,const Vector3d & d0);

    //! Tells if material has history variables 
    // It is used in the Model derived classes
    bool HasHistoryVariables() { return false; };
    
  private:
    Real _kappa;
    
  }; // class AngleSpring
  
}
#endif // _ANGLESPRING_H_
