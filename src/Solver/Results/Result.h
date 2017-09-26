//-*-C++-*-
#ifndef __Result_h__
#define __Result_h__

#include "voom.h"

namespace voom {
  //! Result is a class which contains standard result and defines basic interface to access/update them

  class Result 
  {
  private: 
    Result() {}; // Default constructor is private because it should never be used by a derived class.

  public:

    // Constructor
    Result(int PbDoF): _request(0), _energy(0.0), _pbDoF(PbDoF) {}; // Notice: it does not initialize the field vector

    // Initialization function
    virtual void initializeResults(int NumEntries) = 0;
    void resetResults(int ResetOrder) {
      if (ResetOrder & ENERGY) {
	this->setEnergy(0.0); };
      if (ResetOrder & FORCE) {
	this->resetResidualToZero(); };
      if (ResetOrder & STIFFNESS) {
	this->resetStiffnessToZero(); };
      if (ResetOrder & DMATPROP) {
	this->resetGradgToZero(); 
	this->resetHgToZero();
	this->resetResidualToZero(); };
      }

    // NumMat and PbDoF
    virtual int getNumMatProp() = 0;
    int getPbDoF() { return _pbDoF; };

    // Request, Energy
    void setRequest(int request) { _request = request; };
    int  getRequest() {return _request; };
    void setEnergy(Real energy)  { _energy = energy;   };
    Real getEnergy()  {return _energy; };
    void addEnergy(Real DeltaEnergy) {_energy += DeltaEnergy; };
    
    // Residual
    virtual void resetResidualToZero() = 0;
    virtual void setResidual(int ind, Real value) = 0;
    virtual void addResidual(int ind, Real value) = 0;
    virtual Real getResidual(int ind) = 0;
    
    // Stiffness
    virtual void resetStiffnessToZero() = 0;
    virtual void addStiffness(int indRow, int indCol, Real value) = 0;
    virtual void FinalizeGlobalStiffnessAssembly() = 0;
    virtual Real getStiffness(int indRow, int indCol) = 0;

    // Gradg
    virtual void resetGradgToZero() = 0;
    virtual void addGradg(int ind, Real value) = 0;
    virtual Real getGradg(int ind) = 0;

    // Hg
    virtual void resetHgToZero() = 0;
    virtual void addHg(int indRow, int indCol, Real value) = 0;
    virtual void FinalizeHgAssembly() = 0;
    virtual Real getHg(int indRow, int indCol) = 0;

    
    

  protected:
    //! Solution value at all nodes, local and ghost
    //! Displacement are stored unrolled, [phi_x, phi_y, phi_z]
    int           _pbDoF;
    int           _request;
    Real          _energy;
 
  }; // End of declaration of Result base class

} // namespace voom

#endif
