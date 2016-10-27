/*! 
  \file Holzapfel.h
  \brief Interface for Holzapfel passive myocardium finite deformation
  hyperelasticity model.
*/

#ifndef _HOLZAPFEL_H_
#define _HOLZAPFEL_H_

#include "MechanicsMaterial.h"

namespace voom {
  
  class Holzapfel : public MechanicsMaterial
  {
  public: 
    // Constructors/destructors:
  Holzapfel(int ID):  MechanicsMaterial(ID), _a1(1.0), _a2(1.0),  _a3(1.0), _a4(1.0), _b1(1.0), _b2(1.0), _b3(1.0), _b4(1.0) {};
  Holzapfel(int ID, Real a1, Real a2, Real a3, Real a4, Real b1, Real b2, Real b3, Real b4, vector<Vector3d > Fibers): 
    MechanicsMaterial(ID), _a1(a1), _a2(a2), _a3(a3), _a4(a4), _b1(b1), _b2(b2), _b3(b3), _b4(b4), _fibers(Fibers) {};
  Holzapfel(Holzapfel* BaseMaterial): 
    MechanicsMaterial(BaseMaterial->_matID), _a1(BaseMaterial->_a1), _a2(BaseMaterial->_a2), _a3(BaseMaterial->_a3), _a4(BaseMaterial->_a4),
      _b1(BaseMaterial->_b1), _b2(BaseMaterial->_b2), _b3(BaseMaterial->_b3), _b4(BaseMaterial->_b4), _fibers(BaseMaterial->_fibers) {};

    // Clone
    virtual Holzapfel* clone() const {
      return new Holzapfel(*this);
    }
      
    // Default copy constructor (compiler should already provide exactly this)
    Holzapfel(const Holzapfel & Old): 
    MechanicsMaterial(Old._matID), _a1(Old._a1), _a2(Old._a2), _a3(Old._a3), _a4(Old._a4), _b1(Old._b1), _b2(Old._b2), _b3(Old._b3), _b4(Old._b4), _fibers(Old._fibers) {};

    void setMaterialParameters(const vector<Real > & MatPar) {
      _a1 = MatPar[0];
      _a2 = MatPar[1];
      _a3 = MatPar[2];
      _a4 = MatPar[3];
    }
    void setInternalParameters(const vector<Real > & IntPar) {
      _b1 = IntPar[0];
      _b2 = IntPar[1];
      _b3 = IntPar[2];
      _b4 = IntPar[3];
    }
    void setRegularizationParameters(const vector<Real > &) {}; // No regularization parameters for Holzapfel


    vector<Real > getMaterialParameters() {
      vector<Real > MatProp(4, 0.0);
      MatProp[0] = _a1;
      MatProp[1] = _a2;
      MatProp[2] = _a3;
      MatProp[3] = _a4;
      return MatProp;
    }
    vector<Real > getInternalParameters() {
      vector<Real > IntPar(4, 0.0);
      IntPar[0] = _b1;
      IntPar[1] = _b2;
      IntPar[2] = _b3;
      IntPar[3] = _b4;
      return IntPar;
    }
    vector<Real > getRegularizationParameters() { // No regularization parameters for Holzapfel
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
    //! Data
    Real _a1;
    Real _a2;
    Real _a3;
    Real _a4;
    Real _b1;
    Real _b2;
    Real _b3;
    Real _b4;

    const vector<Vector3d > _fibers;
	
  }; // class Holzapfel
  
}
#endif // _HOLZAPFEL_H_
