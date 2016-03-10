//-*-C++-*-
#ifndef __Filament_Material_h__
#define __Filament_Material_h__

#include "VoomMath.h"

namespace voom
{
  class FilamentMaterial
  {
  public:
    struct Filresults
    {
      // Finite kinematics request and result type                                                                                                      
      Filresults() : k()
      {
	W = 0.0;
	f << 0.0, 0.0, 0.0;
	request = 0;
     };
      
      Real W;
      Vector3d f;
      Real k;      // Initialized automatically to zero                                                                                    
      int request;
      
    }; // 
    
    //! Constructor
    FilamentMaterial(int MatID): _matID(MatID) {;}
    
    //! Destructor
    virtual ~FilamentMaterial(){;}

    //! Clone
    virtual FilamentMaterial* clone() const = 0;

    //! Filament compute function
    virtual void compute(Filresults & R, const Vector3d & d) = 0;
    
    //! Consistency Check for all Mecahnics Material Classes
    void checkConsistency(Filresults & R, const Vector3d & d,
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

    // TODO: DELETE THESE THREE FUNCTIONS
    virtual void setTimestep(double deltaT){;}
    virtual void setActivationMultiplier(double activation){;}
    virtual void updateStateVariables(){;}

  protected:
    int _matID;
    
  private:
    // No material can be initialized without ID
    FilamentMaterial(){;};
    
  }; // class FilamentMaterial

} // namespace voom

#endif
