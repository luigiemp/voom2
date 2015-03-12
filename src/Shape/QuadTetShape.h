//-*-C++-*_
/*!
  \file QuadTetShape.h
  
  \brief Shape function implementation for a quadratic tetrahedral element
*/

#if !defined(__QuadTetShape_h__)
#define __QuadTetShape_h__

#include "Shape.h"

namespace voom{

  class QuadTetShape: public Shape {

  public:
    //! QuadTetShape constructor fills in N and DN
    QuadTetShape(const VectorXd & Point) {
      _N.resize(10, 0.0);
      _DN.resize(10, Vector3d::Zero());
      update(Point);
    }

    //! Update recomputes N and DN at new Point
    void update(const VectorXd & Point);

    //! Get number of shape functions
    uint getShapeFunctionNum() {return 10; };

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
    vector<Vector3d > _DN;
  
  };

} // namespace voom

#endif
