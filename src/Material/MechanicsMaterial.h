//-*-C++-*-
#ifndef __Mechanics_Material_h__
#define __Mechanics_Material_h__

#include "VoomMath.h"

namespace voom
{
  class MechanicsMaterial
  {
  public:
    struct FKresults
    {
      // Finite kinematics request and result type
      FKresults() : K(3,3,3,3), Dmat(), DDmat()
      {
	W = 0.0;
	P << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
	request = 0;
      };
    
      Real W;
      Matrix3d P;
      FourthOrderTensor K;      // Initialized automatically to zero
      ThirdOrderTensor  Dmat;   
      FourthOrderTensor DDmat;  
      int request;
    
    }; // struct FKresults

    //! Constructor
    MechanicsMaterial(int MatID): _matID(MatID) {;}

    //! Destructor
    virtual ~MechanicsMaterial(){;}

    //! Clone
    virtual MechanicsMaterial* clone() const = 0;

    //! Compute function
    virtual void compute(FKresults & R, const Matrix3d & F) = 0;

    //! Consistency Check for all Mecahnics Material Classes
    void checkConsistency(FKresults & R, const Matrix3d & F, 
			  const Real h = 1.0e-7, const Real tol = 1.0e-6);

    //! SetMaterialParameters function
    virtual void setMaterialParameters(const vector<Real > &) = 0;
    virtual void setInternalParameters(const vector<Real > &) = 0;
    virtual void setRegularizationParameters(const vector<Real > &) = 0;

     //! GetMaterialParameters function
    virtual vector<Real > getMaterialParameters() = 0;
    virtual vector<Real > getInternalParameters() = 0;
    virtual vector<Real > getRegularizationParameters() = 0;

    // Get Material ID
    int getMatID() {
      return _matID;
    }

    //! Tells if material has history variables and needs to be replicated at each quadrature point
    // It is used in the Model derived classes
    virtual bool HasHistoryVariables() = 0;

    // Set N - only implemented for anisotropic materials
    virtual void setN(Vector3d N) {}; 

  protected:
    int _matID;
    
  private:
    // No material can be initialized without ID
    MechanicsMaterial(){;};
    
  }; // class MechanicsMaterial

} // namespace voom

#endif
