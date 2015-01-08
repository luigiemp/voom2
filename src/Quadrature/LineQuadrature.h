// -*- C++ -*-
#if !defined(__LineQuadrature_h__)
#define __LineQuadrature_h__

#include "Quadrature.h"

namespace voom
{
  class LineQuadrature: public Quadrature
  {
  public:
    //! default constructor
    LineQuadrature(unsigned int order);
    //! Test
    bool check(unsigned int degree) const;

  private:
    //! No quadrature rule can be instantiated without quadrature order
    LineQuadrature() {};
  };
	
} // namespace voom

#endif // __LineQuadrature_h__
