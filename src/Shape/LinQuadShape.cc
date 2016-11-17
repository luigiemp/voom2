// -*- C++ -*-
#include "LinQuadShape.h"

namespace voom {

//! Update recomputes N and DN at a new Point
  void LinTriShape::update(const VectorXd & Point)
  {
    // Shape Functions
    _N[0] = 0.25 * (1 - Point(0)) * (1 - Point(1));
    _N[1] = 0.25 * (1 + Point(0)) * (1 - Point(1));
    _N[2] = 0.25 * (1 + Point(0)) * (1 + Point(1));
    _N[3] = 0.25 * (1 - Point(0)) * (1 + Point(1));

    // Shape functions derivatives
    _DN[0](0) = -0.25 * (1 - Point(1));
    _DN[0](1) = -0.25 * (1 - Point(0));

    _DN[1](0) =  0.25 * (1 - Point(1));
    _DN[1](1) = -0.25 * (1 + Point(0));

    _DN[2](0) =  0.25 * (1 + Point(1));
    _DN[2](1) =  0.25 * (1 + Point(0));

    _DN[3](0) = -0.25 * (1 + Point(1));
    _DN[3](1) =  0.25 * (1 - Point(0));
  } // update

} // namespace voom
