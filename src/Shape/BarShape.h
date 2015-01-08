//-*-C++-*-
/*!
  \file BarShape.h

  \brief Shape function implementation for bar element.
*/

#if !defined(__BarShape_h__)
#define __BarShape_h__

#include "Shape.h"

namespace voom{

  class BarShape: public Shape {
  public:
    //! BarShape constructor fills in N and DN
    BarShape(const VectorXd Point) {
      _N.resize(2, 0.0);
      _DN.resize(2, 0.0);
      update(Point);
    }

    //! Update recomputes N and DN at new Point
    void update(const VectorXd & Point);

    //! Get number of shape functions
    uint getShapeFunctionNum() {return 2; };

    //! getN returns shape function values at a given node a.
    Real getN(const uint a) {
      return _N[a];
    };
    
    //! GetDN returns shape function derivatives at node a, in direction i.
    Real getDN(const uint a, const uint i = 0) {
      assert(i==0);
      return _DN[a];
    };
    
  private:
    vector<Real > _N;
    vector<Real > _DN;
  
  };

} // namespace voom

#endif
