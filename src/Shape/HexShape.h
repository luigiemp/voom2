//-*-C++-*_
/*!
  \file HexShape.h
  
  \brief Shape function implementation for a hexahedral element
*/

#if !defined(__HexShape_h__)
#define __HexShape_h__

#include "Shape.h"

namespace voom{

  class HexShape: public Shape {

  public:
    //! HexShape constructor fills in N and DN
    HexShape(const VectorXd & Point) {
      _N.resize(8, 0.0);
      _DN.resize(8, Vector3d::Zero());
      update(Point);
    }

    //! Update recomputes N and DN at new Point
    void update(const VectorXd & Point);

    //! Get number of shape functions
    uint getShapeFunctionNum() {return 8; };

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
