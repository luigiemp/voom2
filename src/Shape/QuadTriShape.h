//-*-C++-*_
/*!
  \file QuadTriShape.h
  
  \brief Shape function implementation for a quadratic triangular element
*/

#if !defined(__QuadTriShape_h__)
#define __QuadTriShape_h__

#include "Shape.h"

namespace voom{

  class QuadTriShape: public Shape {

  public:
    //! QuadTriShape constructor fills in N and DN
    QuadTriShape(const VectorXd & Point) {
      _N.resize(6, 0.0);
      _DN.resize(6, Vector2d::Zero());
      update(Point);
    }

    //! Update recomputes N and DN at new Point
    void update(const VectorXd & Point);

    //! Get number of shape functions
    uint getShapeFunctionNum() {return 6; };

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
