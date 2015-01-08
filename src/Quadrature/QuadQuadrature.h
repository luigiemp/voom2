// -*- C++ -*-
/*! 
  \file QuadQuadrature.h
  
  \brief Quadrature rule for a quadrilateral element. The points and weights
  are computed using values from Line Quadrature.
  
*/

#if !defined(__QuadQuadrature_h__)
#define __QuadQuadrature_h__

#include "Quadrature.h"

namespace voom
{
  class QuadQuadrature : public Quadrature
  {
  public:
    //! Constructor
    QuadQuadrature(unsigned int order);
    //! Test
    virtual ~QuadQuadrature() {}
    bool check(unsigned int degree) const;
    
  private:
    //! No quadrature rule can be initialized without quadrature order
    QuadQuadrature() {};
  };
  
} // namespace voom

#endif // __QuadQuadrature_h__
