//-*-C++-*-
/*!
  \file Model.h
  \brief A base model class creates basic entities needs for solving in parallel
  and stores generic boundary conditions. 
 */
#ifndef __Model_h__
#define __Model_h__

#include "voom.h"
#include "State.h"
#include "Result.h"
#include "Body.h"

namespace voom {
  class Model {
  private :
    Model() {}; // Default constructor is private because it should never be used by a derived class.
  
  public:
    //! Constructor to be used
    Model(vector<Body *> Bodies, State* myState): _bodies(Bodies), _myState(myState) {} ;

    //! Destructor
    ~Model() {};

    //! Model compute function
    void compute(Result* R);

    //! Check consistency - common to all models
    void checkConsistency(Result* R, Real perturbationFactor, 
			  int request = 6, Real h = 1e-6, Real tol = 1e-6);

    //! Return vector with all bodies in material
    vector<Body *> getBodies() { return _bodies; } ;

  protected:
    vector<Body *> _bodies;
    State* _myState;
  
  }; // Model

} // namespace voom

#endif
