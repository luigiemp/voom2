//-*-C++-*-
#ifndef __Pair_Material_h__
#define __Pair_Material_h__

#include "VoomMath.h"

namespace voom
{
  class PairMaterial
  {
  public:
    struct PairMresults
    {
      // Finite kinematics request and result type
      PairMresults()
      {
	W = 0.0;
	F = 0.0;
	K = 0.0;
	request = 0;
      };
    
      Real W;
      Real F;  // dW/dr
      Real K;  // ddW/ddr
      int request;
    
    }; // struct PairMresults

    //! Destructor
    virtual ~PairMaterial(){;}

    //! Clone
    virtual PairMaterial* clone() const = 0;

    //! Compute function
    virtual void compute(PairMresults & R, Real r) = 0;

    //! Consistency Check for all Mecahnics Material Classes
    void checkConsistency(PairMresults & R, Real r,
			  Real h = 1.0e-7, Real tol = 1.0e-6);

    //! SetMaterialParameters function
    virtual void setMaterialParameters(const vector<Real > &) = 0;

     //! GetMaterialParameters function
    virtual vector<Real > getMaterialParameters() = 0;
     
  }; // class PairMaterial

} // namespace voom

#endif
