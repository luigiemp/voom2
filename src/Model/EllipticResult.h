//-*-C++-*-
#ifndef __EllipticResult_h__
#define __EllipticResult_h__
#include "Model.h"

namespace voom {
  //! Elliptic Result is a class which contain standard elliptic result and define basic interface to access/update them

  struct EllipticResult 
  {
    // Default constructor
    EllipticResult(): _request(0), _energy(0.0) {};

    // Virtual functions to set Elliptic results members (e.g. to be used in the EllipticModel compute function)
    virtual void setRequest(int request) { _request = request; };
    virtual void setEnergy(Real energy) {_energy = energy; };
    virtual void setResidual(int LocalIndex, Real value) = 0;
    // setStiffness is not needed right now, may be defined later
    // virtual void setStiffness(int LocalRow, int LocalCol, Real value) = 0;

    // Reset function
    virtual void resetResidualToZero() = 0;
    virtual void resetStiffnessToZero() = 0;
    
    // Virtual functions to add to Elliptic result members
    virtual void addEnergy(Real DeltaEnergy) {_energy += DeltaEnergy; };
    virtual void addResidual(int LocalIndex, Real value) = 0;
    virtual void addStiffness(int LocalRow, int LocalCol, Real value) = 0;
    virtual void FinalizeGlobalStiffnessAssembly() = 0;

    // Virtual functions to get Elliptic result members
    virtual int  getRequest() {return _request; };
    virtual Real getEnergy()  {return _energy; };
    virtual Real getResidual(int LocalIndex) = 0;
    virtual Real getStiffness(int LocalRow, int LocalCol) = 0;

  private:
    int  _request;
    Real _energy;

  }; // End of declaration of EllipticResult base class

} // namespace voom

#endif
