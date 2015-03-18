#include "LinTriShape.h"

namespace voom {

//! Update recomputes N and DN at a new Point
  void LinTriShape::update(const VectorXd & Point)
  {
    Real csi = 1.0 - Point(0) - Point(1);
    // Shape functions
    _N[0] = Point(0);
    _N[1] = Point(1);
    _N[2] = csi;

    // Shape functions derivatives
    _DN[0](0) = 1.0;
    _DN[0](1) = 0.0;

    _DN[1](0) = 0.0;
    _DN[1](1) = 1.0;

    _DN[2](0) =-1.0;
    _DN[2](1) =-1.0;

  } // update
						       
} // namespace voom
