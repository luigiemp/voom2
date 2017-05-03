//-*-C++-*-

#ifndef __Body_h__
#define __Body_h__

#include "voom.h"
#include "Mesh.h"
#include "Result.h"

namespace voom {
  class Body {
  private :
    Body() {}; // Default constructor is private because it should never be used by a derived class.
  
  public:
    //! Constructor to be used
    Body(Mesh* myMesh, const int NodeDoF)
      : _myMesh(myMesh), _nodeDoF(NodeDoF) {
    }

    //! Destructor
    virtual ~Body() {};

    //! Body compute function
    virtual void compute(Result* R) = 0;

    //! Check consistency - common to all bodies
    void checkConsistency(Result* R, Real perturbationFactor, int request = 6,
			  Real h = 1e-6, Real tol = 1e-6);

    //! Return number of unique material pointers
    virtual int getNumMat() = 0;

    //! Get DoF per Node
    int getDoFperNode() {
      return _nodeDoF;
    }
    
    //! GetMesh
    Mesh * getMesh() {
      return _myMesh;
    }

    //! Initialize field
    // From constant value
    virtual void initializeField(Result* R, Real value = 1.0) = 0;
 
    //! From array
    virtual void initializeField(Result* R, const vector<Real > & Field) = 0;

    //! Finalize compute
    virtual void FinalizeCompute() = 0;



  protected:
    Mesh*        _myMesh;

    //! DoF per node
    int         _nodeDoF;
  
  }; // Body

} // namespace voom

#endif
