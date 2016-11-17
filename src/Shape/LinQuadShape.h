//-*-C++-*_
/*!
  \file LinQuadShape.h

  \brief Shape function implementation for a linear quadrilateral element
*/

#if !defined(__LinQuadShape_h__)
#define __LinQuadShape_h__

#include "Shape.h"

namespace voom{

  class LinQuadShape: public Shape {

  public:
    //! LinTriShape constructor fills in N and DN
    ShapeQ4(const Vector2d & Point) {
      _N.resize(4, 0.0);
      _DN.resize(4, Vector2d::Zero());
      update(Point);
    }

    //! Update recomputes N and DN at new Point
    void update(const VectorXd & Point);

    //! Get number of shape functions
    uint getShapeFunctionNum() {return 4; };

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
