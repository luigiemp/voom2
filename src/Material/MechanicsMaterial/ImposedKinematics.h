/*! 
  \file ImposedKinematics.h
  \brief Interface for a finite deformation hyperelastic material
         to impose specified kinematics along chosen directions
*/

#ifndef _IMPOSEDKINEMATICS_H_
#define _IMPOSEDKINEMATICS_H_

#include "MechanicsMaterial.h"

namespace voom {
  
  class ImposedKinematics : public MechanicsMaterial
  {
  public: 
    // Constructors/destructors:
  ImposedKinematics(int ID, vector<Real > Alphas, vector<Real > Stretches, vector<Vector3d > Directions, Real Beta):  
    MechanicsMaterial(ID), _alphas(Alphas), _stretches(Stretches), _directions(Directions), _beta(Beta) {};
  
  ImposedKinematics(ImposedKinematics* BaseMaterial): 
    MechanicsMaterial(BaseMaterial->_matID), _alphas(BaseMaterial->_alphas), _stretches(BaseMaterial->_stretches), _directions(BaseMaterial->_directions), _beta(BaseMaterial->_beta) {};

    // Clone
    virtual ImposedKinematics* clone() const {
      return new ImposedKinematics(*this);
    }
      
    // Default copy constructor (compiler should already provide exactly this)
    ImposedKinematics(const ImposedKinematics & Old): 
    MechanicsMaterial(Old._matID), _alphas(Old._alphas), _stretches(Old._stretches), _directions(Old._directions), _beta(Old._beta) {};

    void setMaterialParameters(const vector<Real > & MatPar) {
      _alphas = MatPar;
    }
    void setInternalParameters(const vector<Real > & IntPar) { 
      _stretches = IntPar;
    }
    void setRegularizationParameters(const vector<Real > & RegPar) { 
      _beta   = RegPar[0];
    };

    vector<Real > getMaterialParameters() {
      return _alphas;
    }
    vector<Real > getInternalParameters() {
      return _stretches;
    }
    vector<Real > getRegularizationParameters() {
      vector<Real > RegProp(1, 0.0);
      RegProp[0] = _beta;
      return RegProp;      
    };

    vector<Vector3d> getDirectionVectors() {
      return _directions;
    }

    // Operators
    //! Based on deformation gradient tensor F, calculates state of material
    void compute(FKresults & R, const Matrix3d & F);

    //! Tells if material has history variables and needs to be duplicated at each quadrature point
    // It is used in the Model derived classes
    bool HasHistoryVariables() { return false; };
    
  private:
    //! Data
    vector<Real > _alphas;
    vector<Real > _stretches;
    vector<Vector3d > _directions;
    Real _beta;
	
  }; // class ImposedKinematics
  
}
#endif // _IMPOSEDKINEMATICS_H_
