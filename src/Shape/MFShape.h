//-*-C++-*-
/*!
  \file MFShape.h

  \brief Parent Meshfree shape function class. Implements a check partition of
  of unity function which also includes reciprocating conditions
*/

#if !defined(__MFShape_h__)
#define __MFShape_h__

#include "Shape.h"

namespace voom{
  class MFShape: public Shape {
  protected:
    //! Nodes which lay withing support radius
    const vector<VectorXd>&  _nodes;

    //! Shape functions
    vector<Real>             _N;

    //! Shapefunction derivatives
    vector<VectorXd>         _DN;

    //! Number of shape functions
    uint                     _shapeNum;

  public:
    //! Default Constructor
    MFShape(const vector<VectorXd>& nodes);

    //! Default Destructor
    ~MFShape() {;}

    //! Get number of shape functions
    uint getShapeFunctionNum() {return _shapeNum; };

    //! getN returns shape function values at a given node a.
    Real getN(const uint a) {
      return _N[a];
    };
    
    //! GetDN returns shape function derivatives at node a, in direction i.
    Real getDN(const uint a, const uint i) {
      return _DN[a](i);
    };

    //! Update recomputes N and DN at new Point
    virtual void update(const VectorXd & Point) = 0;

    //! Overloaded chech partition of unity function
    bool checkReproducingCondition(VectorXd Point, const Real tol = 1e-7);
  };
}

#endif
