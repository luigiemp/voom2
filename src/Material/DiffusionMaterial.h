//-*-C++-*-
#ifndef __Diffusion_Material_h__
#define __Diffusion_Material_h__

#include "VoomMath.h"

namespace voom
{
  class DiffusionMaterial
  {
    
  public:
    struct DiffusionResults
    {
      // DiffusionResults
      DiffusionResults()
      {
	A = Matrix3d::Zero();
      };
      
      Matrix3d A; // conductivity matrix
    }; // struct DiffusionResults

    //! Constructor
    DiffusionMaterial(){;}

    //! Destructor
    virtual ~DiffusionMaterial(){;}

    //! Clone
    virtual DiffusionMaterial* clone() const = 0;

    //! Tells if material has history variables and needs to be replicated at each quadrature point
    // It is used in the Model derived classes
    virtual bool HasHistoryVariables() = 0;

    //! Compute function - useful only if A depends on scalar field (e.g. T)
    virtual void compute(DiffusionResults & R, Real scalar = 0.0) = 0;
    
  }; // class DiffusionMaterial

} // namespace voom

#endif
