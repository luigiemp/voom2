//-*-C++-*-
#ifndef __MechanicsLoad_h__
#define __MechanicsLoad_h__

#include "VoomMath.h"

namespace voom
{
  class MechanicsLoad
  {
  public:
    struct SourceLoadResults // No follower forces for now
    {
      SourceLoadResults()
      {
	W = 0.0;
	request = 0;
      };
    
      Real W;
      vector<Real > F;
      int request;
    
    }; // struct SourceLoadResults

    //! Constructor
    MechanicsLoad(){;}

    //! Destructor
    virtual ~MechanicsLoad(){;}

    //! Compute function
    virtual void compute(SourceLoadResults & R, vector<Real > u) = 0;

    //! Consistency Check for all Mecahnics Material Classes
    void checkConsistency(SourceLoadResults & R, vector<Real > u, 
			  const Real h = 1.0e-8, const Real tol = 1.0e-7);

    //! Return GLOBAL DoF to which load is applied
    vector<uint > getNeumannDoF() {
      return _newmannDoF;
    }

  private:
    vector<uint > _newmannDoF;
    
  }; // class MechanicsLoad

} // namespace voom

#endif
