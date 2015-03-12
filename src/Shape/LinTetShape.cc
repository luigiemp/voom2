#include "LinTetShape.h"

namespace voom {

//! Update recomputes N and DN at a new Point
  void LinTetShape::update(const VectorXd & Point)
  {
    // Shape functions
    _N[0] = Point(0);
    _N[1] = Point(1);
    _N[2] = Point(2);
    _N[3] = 1.0 - Point(0) - Point(1) - Point(2);

    // Shape functions derivatives
    _DN[0](0) = 1.0;
    _DN[0](1) = 0.0;
    _DN[0](2) = 0.0;

    _DN[1](0) = 0.0;
    _DN[1](1) = 1.0;
    _DN[1](2) = 0.0;

    _DN[2](0) = 0.0;
    _DN[2](1) = 0.0;
    _DN[2](2) = 1.0;

    _DN[3](0) =-1.0;
    _DN[3](1) =-1.0;
    _DN[3](2) =-1.0;
    
  } // update
						       
} // namespace voom
