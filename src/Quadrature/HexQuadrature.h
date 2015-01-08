//-*-C++-*-
#if !defined(__HexQuadrature_h__)
#define __HexQuadrature_h__

#include "Quadrature.h"

namespace voom
{
  class HexQuadrature : public Quadrature
  {
    public:
      //! Constructors
      HexQuadrature(unsigned int order);
      //! Test
      bool check(unsigned int degree) const;

    private:
      //! No quadrature rule can be instantiated without quadrature order
      HexQuadrature() {};
    
  }; // class HexQuadrature
	
} // namespace voom

#endif // __HexQuadrature_h__
