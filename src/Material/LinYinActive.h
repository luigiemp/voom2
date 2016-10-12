/*!
  \file LinYinActive.h
  \brief Interface for the Humphrey Model for the passive myocardium
*/

#ifndef _LINYINACTIVE_H_
#define _LINYINACTIVE_H_

#include "MechanicsMaterial.h"

namespace voom {

  class LinYinActive : public MechanicsMaterial
  {
  public:
    // Constructors/destructors:
    LinYinActive(int ID): MechanicsMaterial(ID), _C0(1.0), _C1(1.0), _C2(1.0), _C3(1.0), _C4(1.0), _C5(1.0) {};
    LinYinActive(int ID, Real C0, Real C1, Real C2, Real C3, Real C4, Real C5, vector<Vector3d> Fibers): MechanicsMaterial(ID), _C0(C0), _C1(C1), _C2(C2), _C3(C3), _C4(C4), _C5(C5), _fibers(Fibers) {};
    LinYinActive(LinYinActive* BaseMaterial):
    MechanicsMaterial(BaseMaterial->_matID), _C0(BaseMaterial->_C0), _C1(BaseMaterial->_C1), _C2(BaseMaterial->_C2), _C3(BaseMaterial->_C3), _C4(BaseMaterial->_C4), _C5(BaseMaterial->_C5), _fibers(BaseMaterial->_fibers) {};

    // Clone
    virtual LinYinActive* clone() const {
      return new LinYinActive(*this);
    }

    // Default copy constructor (compiler should already provide exactly this)
    LinYinActive(const LinYinActive & Old):
    MechanicsMaterial(Old._matID), _C0(Old._C0), _C1(Old._C1), _C2(Old._C2), _C3(Old._C3), _C4(Old._C4), _C5(Old._C5), _fibers(Old._fibers) {};

    Real getC0() {return _C0;};
    Real getC1() {return _C1;};
    Real getC2() {return _C2;};
    Real getC3() {return _C3;};
    Real getC4() {return _C4;};
    Real getC5() {return _C5;};

    void setC0(Real C0) {_C0 = C0;};
    void setC1(Real C1) {_C1 = C1;};
    void setC2(Real C2) {_C2 = C2;};
    void setC3(Real C3) {_C3 = C3;};
    void setC4(Real C4) {_C4 = C4;};
    void setC5(Real C5) {_C5 = C5;};

    void setMaterialParameters(const vector<Real > & MatParam) {
      _C0 = MatParam[0];
      _C1 = MatParam[1];
      _C2 = MatParam[2];
      _C3 = MatParam[3];
      _C4 = MatParam[4];
      _C5 = MatParam[5];
    }
    void setInternalParameters(const vector<Real > & IntParam) {}; // No internal parameters for LinYinActive
    void setRegularizationParameters(const vector<Real > &)    {}; // LinYinActive does not have regularization parameters

    vector<Real > getMaterialParameters() {
      vector<Real > MatProp(6, 0.0);
      MatProp[0] = _C0;
      MatProp[1] = _C1;
      MatProp[2] = _C2;
      MatProp[3] = _C3;
      MatProp[4] = _C4;
      MatProp[5] = _C5;
      return MatProp;
    }
    vector<Real > getInternalParameters() { // No internal parameters for LinYinActive
      vector<Real > IntParam;
      return IntParam;
    }
    vector<Real > getRegularizationParameters() { // No regularization parameters for LinYinActive
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
    Real _C0;
    Real _C1;
    Real _C2;
    Real _C3;
    Real _C4;
    Real _C5;

    const vector<Vector3d > _fibers;

  }; // class LinYinActive

}
#endif // _LINYINACTIVE_H_
