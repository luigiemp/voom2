//-*-C++-*-
#ifndef __Result_h__
#define __Result_h__
#include "Model.h"

namespace voom {
  //! Result is a class which contains standard result and defines basic interface to access/update them

  class Result 
  {

  public:

    // Default constructor
    Result(): _request(0), _energy(0.0) {}; // Notice: it does not initialize the field vector

    // Return NumMat and PbDoF
    virtual int getNumMatProp() = 0;
    virtual int getPbDoF() = 0;

    // Virtual functions to set results members (e.g. to be used in the EllipticModel compute function)
    virtual void setRequest(int request) { _request = request; };
    virtual void setEnergy(Real energy) {_energy = energy; };
    virtual void setResidual(int ind, Real value) = 0;

    // Reset function
    virtual void resetResidualToZero() = 0;
    virtual void resetStiffnessToZero() = 0;
    virtual void resetGradgToZero() = 0;
    virtual void resetHgToZero() = 0;
    
    // Virtual functions to add to result members
    virtual void addEnergy(Real DeltaEnergy) {_energy += DeltaEnergy; };
    virtual void addResidual(int ind, Real value) = 0;
    virtual void addStiffness(int indRow, int indCol, Real value) = 0;
    virtual void FinalizeGlobalStiffnessAssembly() = 0;
    virtual void setStiffnessFromTriplets(vector<Triplet<Real > > &) = 0;

    virtual void addGradg(int ind, Real value) = 0;
    virtual void addHg(int indRow, int indCol, Real value) = 0;
    virtual void setHgFromTriplets(vector<Triplet<Real > > &) = 0;

    // Virtual functions to get result members
    virtual int  getRequest() {return _request; };
    virtual Real getEnergy()  {return _energy; };
    virtual Real getResidual(int ind) = 0;
    virtual Real getStiffness(int indRow, int indCol) = 0;
    virtual Real getGradg(int ind) = 0;
    virtual Real getHg(int indRow, int indCol) = 0;

    // Functions relating to field
    virtual void fieldResize(int size) {_field.resize(size);};
    virtual void initializeField(const vector<Real> field){_field = field;};
    virtual int getFieldSize() {return _field.size();};


    int  _request;
    Real _energy;
 
    //! Solution value at all nodes, local and ghost
    //! Displacement are stored unrolled, [phi_x, phi_y, phi_z]
    vector<Real > _field;

  }; // End of declaration of Result base class

} // namespace voom

#endif
