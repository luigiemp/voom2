/*! 
  \file IsotropicDiffusion.h
  \brief Interface for an isotropic diffusion material
*/

#ifndef _Isotropic_Diffusion_
#define _Isotropic_Diffusion_

#include "DiffusionMaterial.h"

namespace voom {
  
  class IsotropicDiffusion : public DiffusionMaterial
  {
    public: 
      // Constructors/destructors:
      IsotropicDiffusion(): _k(1.0) {};
      IsotropicDiffusion(Real k): _k(k) {};
      IsotropicDiffusion(IsotropicDiffusion* BaseMaterial): _k(BaseMaterial->_k) {};
    
      // Compute function only fills in A matrix in R
      void compute(DiffusionResults & R, Real T = 0.0) {
	R.A << _k, 0.0, 0.0, 0.0, _k, 0.0, 0.0, 0.0, _k;
      }
      
      // Clone
      virtual IsotropicDiffusion* clone() const {
	return new IsotropicDiffusion(*this);
      }
      
      // Default copy constructor (compiler should already provide exactly this)
      IsotropicDiffusion(const IsotropicDiffusion & Old): _k(Old._k) {};

      //! Tells if material has history variables and needs to be duplicated at each quadrature point
      // It is used in the Model derived classes
      bool HasHistoryVariables() { return false; };

   private:
      //! Data
      Real _k; // conductivity

  }; // class IsotropicDiffusion 

}
#endif // _Isotropic_Diffusion_
