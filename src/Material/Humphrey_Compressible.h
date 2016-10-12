/*!
  \file Humphrey.h
  \brief Interface for the Humphrey Model for the passive myocardium with a Volumetric Response
*/

#ifndef _HUMPHREYCOMPRESSIBLE_H_
#define _HUMPHREYCOMPRESSIBLE_H_

#include "MechanicsMaterial.h"

namespace voom {

  class Humphrey_Compressible : public MechanicsMaterial
  {
  public:
    // Constructors/destructors:
  Humphrey_Compressible(int ID): MechanicsMaterial(ID), _C1(1.0), _C2(1.0), _C3(1.0), _C4(1.0), _C5(1.0), _C6(1.0), _C7(1.0) {};
  Humphrey_Compressible(int ID, Real C1, Real C2, Real C3, Real C4, Real C5, Real C6, Real C7, vector<Vector3d> Fibers): MechanicsMaterial(ID), _C1(C1), _C2(C2), _C3(C3), _C4(C4), _C5(C5), _C6(C6), _C7(C7), _fibers(Fibers) {};
    Humphrey_Compressible(Humphrey_Compressible* BaseMaterial):
    MechanicsMaterial(BaseMaterial->_matID), _C1(BaseMaterial->_C1), _C2(BaseMaterial->_C2), _C3(BaseMaterial->_C3), _C4(BaseMaterial->_C4), _C5(BaseMaterial->_C5), _C6(BaseMaterial->_C6), _C7(BaseMaterial->_C7), _fibers(BaseMaterial->_fibers) {};

    // Clone
    virtual Humphrey_Compressible* clone() const {
      return new Humphrey_Compressible(*this);
    }

    // Default copy constructor (compiler should already provide exactly this)
    Humphrey_Compressible(const Humphrey_Compressible & Old):
    MechanicsMaterial(Old._matID), _C1(Old._C1), _C2(Old._C2), _C3(Old._C3), _C4(Old._C4), _C5(Old._C5), _C6(Old._C6), _C7(Old._C7), _fibers(Old._fibers) {};

    Real getC1() {return _C1;};
    Real getC2() {return _C2;};
    Real getC3() {return _C3;};
    Real getC4() {return _C4;};
    Real getC5() {return _C5;};
    Real getC6() {return _C6;};
    Real getC7() {return _C7;};

    void setC1(Real C1) {_C1 = C1;};
    void setC2(Real C2) {_C2 = C2;};
    void setC3(Real C3) {_C3 = C3;};
    void setC4(Real C4) {_C4 = C4;};
    void setC5(Real C5) {_C5 = C5;};
    void setC6(Real C6) {_C6 = C6;};
    void setC7(Real C7) {_C7 = C7;};

    void setMaterialParameters(const vector<Real > & MatParam) {
      _C1 = MatParam[0];
      _C2 = MatParam[1];
      _C3 = MatParam[2];
      _C4 = MatParam[3];
      _C5 = MatParam[4];
      _C6 = MatParam[5];
      _C7 = MatParam[6];
    }
    void setInternalParameters(const vector<Real > & IntParam) {}; // No internal parameters for Humphrey
    void setRegularizationParameters(const vector<Real > &)    {}; // Humphrey does not have regularization parameters

    vector<Real > getMaterialParameters() {
      vector<Real > MatProp(7, 0.0);
      MatProp[0] = _C1;
      MatProp[1] = _C2;
      MatProp[2] = _C3;
      MatProp[3] = _C4;
      MatProp[4] = _C5;
      MatProp[5] = _C6;
      MatProp[6] = _C7;
      return MatProp;
    }
    vector<Real > getInternalParameters() { // No internal parameters for Humphrey
      vector<Real > IntParam;
      return IntParam;
    }
    vector<Real > getRegularizationParameters() { // No regularization parameters for Humphrey
      vector<Real > RegParam;
      return RegParam;
    }


    // Operators
    //! Based on deformation gradient tensor F, calculates state of material
    void compute(FKresults & R, const Matrix3d & F);;

    //! Tells if material has history variables and needs to be duplicated at each quadrature point
    // It is used in the Model derived classes
    bool HasHistoryVariables() { return false; };

  private:
    //! Data
    Real _C1;
    Real _C2;
    Real _C3;
    Real _C4;
    Real _C5;
    Real _C6;
    Real _C7;

    const vector<Vector3d > _fibers;

  }; // class Humphrey_Compressible

}
#endif // _HUMPHREYCOMPRESSIBLE_H_
