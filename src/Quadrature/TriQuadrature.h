//-*-C++-*-
#if !defined(__TriQuadrature_h__)
#define __TriQuadrature_h__

#include "Quadrature.h"

namespace voom
{
  class TriQuadrature : public Quadrature
  {
    public:
      //! Constructors
      TriQuadrature(unsigned int order);
      //! Test
      bool check(unsigned int degree) const;

    private:
      //! No quadrature rule can be instantiated without quadrature order
      TriQuadrature() {};
    
  }; // class TriQuadrature
	
} // namespace voom

#endif // __TriQuadrature_h__
