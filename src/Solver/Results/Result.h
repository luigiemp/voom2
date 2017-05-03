//-*-C++-*-
#ifndef __Result_h__
#define __Result_h__

namespace voom {
  //! Result is a class which contains standard result and defines basic interface to access/update them

  class Result 
  {
  private : 
    Result() {}; // Default constructor is private because it should never be used by a derived class.

  public:

    // Default constructor
    Result(int PbDoF): _request(0), _energy(0.0), _pbDoF(PbDoF) {}; // Notice: it does not initialize the field vector

    // Return NumMat and PbDoF
    virtual int getNumMatProp() = 0;
    int getPbDoF() { return _pbDoF; };

    // Virtual functions to set results members (e.g. to be used in the EllipticModel compute function)
    void setRequest(int request) { _request = request; };
    void setEnergy(Real energy) {_energy = energy; };
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

    virtual void addGradg(int ind, Real value) = 0;
    virtual void addHg(int indRow, int indCol, Real value) = 0;
    virtual void FinalizeHgAssembly() = 0;

    // Virtual functions to get result members
    int  getRequest() {return _request; };
    Real getEnergy()  {return _energy; };
    virtual Real getResidual(int ind) = 0;
    virtual Real getStiffness(int indRow, int indCol) = 0;
    virtual Real getGradg(int ind) = 0;
    virtual Real getHg(int indRow, int indCol) = 0;

    // Functions relating to field
    virtual void fieldResize(int size) = 0;
    virtual void initializeField(const vector<Real > field) = 0;
    virtual void initializeField(int ind, Real value) = 0;
    virtual int  getFieldSize() = 0;
    virtual void linearizedUpdate(const Real* locaValues, Real fact) = 0; 
    virtual void linearizedUpdate(const int id, const int LocalDoF, const Real value, const int nodeDoF = 3) = 0;
    virtual void setField(uint GlobalDoF, Real value) = 0; 
    virtual void printField(const int nodeDoF = 3) = 0;
    virtual void writeField(string OutputFile, int step) = 0;

    //! Solution value at all nodes, local and ghost
    //! Displacement are stored unrolled, [phi_x, phi_y, phi_z]
    int           _pbDoF;
    int           _request;
    Real          _energy;
 
  }; // End of declaration of Result base class

} // namespace voom

#endif
