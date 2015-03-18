#include "QuadTriShape.h"

namespace voom {

//! Update recomputes N and DN at a new Point
  void QuadTriShape::update(const VectorXd & Point)
  {
    Real csi = 1.0 - Point(0) - Point(1);
    // Shape functions
    _N[0] = Point(0)*(2.0*Point(0) - 1.0);
    _N[1] = Point(1)*(2.0*Point(1) - 1.0);
    _N[2] = csi*(2.0*csi - 1.0);
    _N[3] = 4.0*Point(0)*Point(1);
    _N[4] = 4.0*Point(1)*csi;
    _N[5] = 4.0*Point(0)*csi;

    // Shape functions derivatives
    _DN[0](0) = 4.0*Point(0) - 1.0;
    _DN[0](1) = 0.0;

    _DN[1](0) = 0.0;
    _DN[1](1) = 4.0*Point(1) - 1.0;

    _DN[2](0) = 1.0 - 4.0*csi;
    _DN[2](1) = 1.0 - 4.0*csi;

    _DN[3](0) = 4.0*Point(1);
    _DN[3](1) = 4.0*Point(0);

    _DN[4](0) =-4.0*Point(1);
    _DN[4](1) = 4.0*(csi - Point(1));

    _DN[5](0) = 4.0*(csi - Point(0));
    _DN[5](1) =-4.0*Point(0);

  } // update
						       
} // namespace voom
