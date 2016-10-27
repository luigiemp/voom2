/*! 
  \file Guccione.h
  \brief Interface for Guccione passive myocardium finite deformation
  hyperelasticity model.
*/

#ifndef _GUCCIONE_H_
#define _GUCCIONE_H_

#include "MechanicsMaterial.h"

namespace voom {
  
  class Guccione : public MechanicsMaterial
  {
  public: 
    // Constructors/destructors:
  Guccione(int ID):  MechanicsMaterial(ID), _a1(1.0), _b1(1.0),  _b2(1.0), _b3(1.0) {};
  Guccione(int ID, Real a1, Real b1, Real b2, Real b3, vector<Vector3d > Fibers): 
    MechanicsMaterial(ID), _a1(a1), _b1(b1), _b2(b2), _b3(b3), _fibers(Fibers) {};
  Guccione(Guccione* BaseMaterial): 
    MechanicsMaterial(BaseMaterial->_matID), _a1(BaseMaterial->_a1), _b1(BaseMaterial->_b1), _b2(BaseMaterial->_b2), _b3(BaseMaterial->_b3),
      _fibers(BaseMaterial->_fibers) {};

    // Clone
    virtual Guccione* clone() const {
      return new Guccione(*this);
    }
      
    // Default copy constructor (compiler should already provide exactly this)
    Guccione(const Guccione & Old): 
    MechanicsMaterial(Old._matID), _a1(Old._a1), _b1(Old._b1), _b2(Old._b2), _b3(Old._b3), _fibers(Old._fibers) {};

    void setMaterialParameters(const vector<Real > & MatPar) {
      _a1 = MatPar[0];
    }
    void setInternalParameters(const vector<Real > & IntPar) {
      _b1 = IntPar[0];
      _b2 = IntPar[1];
      _b3 = IntPar[2];
    }
    void setRegularizationParameters(const vector<Real > &) {}; // No regularization parameters for Guccione


    vector<Real > getMaterialParameters() {
      vector<Real > MatProp(1, 0.0);
      MatProp[0] = _a1;
      return MatProp;
    }
    vector<Real > getInternalParameters() {
      vector<Real > IntPar(3, 0.0);
      IntPar[0] = _b1;
      IntPar[1] = _b2;
      IntPar[2] = _b3;
      return IntPar;
    }
    vector<Real > getRegularizationParameters() { // No regularization parameters for Guccione
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
    Real _b1;
    Real _b2;
    Real _b3;

    const vector<Vector3d > _fibers;
	
  }; // class Guccione
  
}
#endif // _GUCCIONE_H_
