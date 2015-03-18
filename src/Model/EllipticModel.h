//-*-C++-*-
#ifndef __EllipticModel_h__
#define __EllipticModel_h__

#include "Model.h"
#include "EllipticResult.h"
// #include "EpetraEllipticResult.h"

namespace voom {
  //! Elliptic model is the parent class for FiniteMechanics, LinearizedMechanics, Poisson  
  class EllipticModel: public Model 
  {
  public :
    // Default constructor
    EllipticModel(Mesh* aMesh, const uint NodeDoF) 
      : Model(aMesh, NodeDoF) {};

    //! Model comput function
    virtual void compute(EllipticResult & R) = 0;

    //! Check consistency - common to all elliptic models
    void checkConsistency(EllipticResult & R, Real perturbationFactor, int request = 6,
			  Real h = 1e-6, Real tol = 1e-6);

    //! Return number of unique material pointers
    virtual uint getNumMat() = 0;

    //! Return component of residual vector given nodal ID and number of DoF index
    // This function would be needed if we decide to NOT store DoF sequentially in the residual vector
    // Right now we assume that DoF are stored sequentially: all DoF for the first node followed by all DoF for the second node and so on.
    // virtual Real getResidualComponent(EllipticResult & R, int ID, int dim) = 0;
    // virtual Real getStiffnessComponent(EllipticResult & R, int RowID, int RowIndex, int ColID, int ColIndex) = 0;
    
  }; // End of declaration of EllipticModel base class

} // namespace voom

#endif
