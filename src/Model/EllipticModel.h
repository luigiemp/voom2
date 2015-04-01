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

    // SetPrev field
    virtual void setPrevField() = 0;
    
  }; // End of declaration of EllipticModel base class

} // namespace voom

#endif
