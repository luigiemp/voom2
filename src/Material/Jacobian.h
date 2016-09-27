/*! 
  \file Jacobian.h
  \brief Interface to a material used to compute the ratio between the element reference and current volume
*/

#ifndef _JACOBIAN_H_
#define _JACOBIAN_H_

#include "MechanicsMaterial.h"

namespace voom {
  
  class Jacobian : public MechanicsMaterial
  {
  public: 
    // Constructors/destructors:
    Jacobian(int ID): MechanicsMaterial(ID) {}; 
    Jacobian(Jacobian* BaseMaterial): 
    MechanicsMaterial(BaseMaterial->_matID) {};

    // Clone
    virtual Jacobian* clone() const {
      return new Jacobian(*this);
    }
      
    // Default copy constructor (compiler should already provide exactly this)
    Jacobian(const Jacobian & Old): 
    MechanicsMaterial(Old._matID) {};

    void setMaterialParameters(const vector<Real > & MatParam) {}; // No material parameters for Jacobian
    void setInternalParameters(const vector<Real > & IntParam) {}; // No internal parameters for Jacobian
    void setRegularizationParameters(const vector<Real > &)    {}; // No regularization parameters for Jacobian

    vector<Real > getMaterialParameters() { // No material parameters for Jacobian
      vector<Real > MatProp;
      return MatProp;
    }
    vector<Real > getInternalParameters() { // No internal parameters for Jacobian
      vector<Real > IntParam;
      return IntParam;
    }
    vector<Real > getRegularizationParameters() { // No regularization parameters for Jacobian
      vector<Real > RegParam;
      return RegParam;
    }

    
    // Operators
    //! Based on deformation gradient tensor F, calculates state of material
    void compute(FKresults & R, const Matrix3d & F);

    //! Tells if material has history variables and needs to be duplicated at each quadrature point
    // It is used in the Model derived classes
    bool HasHistoryVariables() { return false; };
    
  private:
	
  }; // class Jacobian
  
}
#endif // _JACOBIAN_H_
