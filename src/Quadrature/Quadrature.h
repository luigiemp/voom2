/*! 
  \file Quadrature.h

  \brief Virtual Base Class for the concept of a quadrature rule.

*/

#if !defined(__Quadrature_h__)
#define __Quadrature_h__

#include "VoomMath.h"

namespace voom
{
  class Quadrature
  {
    public:
      virtual ~Quadrature() {;};
      // Operators
      const vector<Real >     & getQuadWeights() const {return _quadWeights; };
      const vector<VectorXd > & getQuadPoints()  const {return _quadPoints;  };
      // Test
      virtual bool check(unsigned int degree) const = 0;

    protected:
      vector<Real >     _quadWeights;
      vector<VectorXd > _quadPoints;
      
  }; // Class quadrature
	
}; // namespace voom

#endif // __Quadrature_h__
