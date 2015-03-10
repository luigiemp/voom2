//-*-C++-*-
#if !defined(__TetQuadrature_h__)
#define __TetQuadrature_h__

#include "Quadrature.h"

namespace voom
{
  class TetQuadrature : public Quadrature
  {
    public:
      //! Constructors
      TetQuadrature(unsigned int order);
      //! Test
      bool check(unsigned int degree) const;

    private:
      //! No quadrature rule can be instantiated without quadrature order
      TetQuadrature() {};
    
  }; // class TetQuadrature
	
} // namespace voom

#endif // __TetQuadrature_h__
