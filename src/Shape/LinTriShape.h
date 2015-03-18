//-*-C++-*_
/*!
  \file LinTriShape.h
  
  \brief Shape function implementation for a linear triangular element
*/

#if !defined(__LinTriShape_h__)
#define __LinTriShape_h__

#include "Shape.h"

namespace voom{

  class LinTriShape: public Shape {

  public:
    //! LinTriShape constructor fills in N and DN
    LinTriShape(const Vector2d & Point) {
      _N.resize(3, 0.0);
      _DN.resize(3, Vector2d::Zero());
      update(Point);
    }

    //! Update recomputes N and DN at new Point
    void update(const VectorXd & Point);

    //! Get number of shape functions
    uint getShapeFunctionNum() {return 3; };

    //! getN returns shape function values at a given node a.
    Real getN(const uint a) {
      return _N[a];
    };
    
    //! GetDN returns shape function derivatives at node a, in direction i.
    Real getDN(const uint a, const uint i) {
      return _DN[a](i);
    };
    
  private:
    vector<Real > _N;
    vector<Vector2d > _DN;
  
  };

} // namespace voom

#endif
