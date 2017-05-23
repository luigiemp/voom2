//-*-C++-*-
/*!
  \file Model.h
  \brief A base model class creates basic entities needs for solving in parallel
  and stores generic boundary conditions. 
 */
#ifndef __Model_h__
#define __Model_h__

#include "voom.h"
#include "Body.h"
#include "Result.h"

namespace voom {
  class Model {
  private :
    Model() {}; // Default constructor is private because it should never be used by a derived class.
  
  public:
    //! Constructor using an input file ... to be changed after implementation of input class - no longer current
    // Model(Mesh* myMesh, const string inputFile, const uint NodeDoF);

    //! Constructor to be used
    Model(vector<Body *> Bodies, Result* ModelResult): _bodies(Bodies) _modelResult(ModelResult) {} ;

    //! Destructor
    ~Model() {};

    //! Model compute function
    void compute(Result* R);

    //! Check consistency - common to all models
    void checkConsistency(Result* R,  int nodeNum, int nodeDoF, 
			  Real perturbationFactor, 
			  int request = 6, Real h = 1e-6, Real tol = 1e-6);

  protected:
    vector<Body *> _bodies;
    Result*        _modelResult;
  
  }; // Model

} // namespace voom

#endif
