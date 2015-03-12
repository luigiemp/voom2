#include "QuadTetShape.h"

namespace voom {

//! Update recomputes N and DN at a new Point
  void QuadTetShape::update(const VectorXd & Point)
  {
    Real csi = 1.0 - Point(0) - Point(1) - Point(2);
    // Shape functions
    _N[0] = Point(0)*(2.0*Point(0) - 1.0);
    _N[1] = Point(1)*(2.0*Point(1) - 1.0);
    _N[2] = Point(2)*(2.0*Point(2) - 1.0);
    _N[3] = csi*(2.0*csi - 1.0);
    _N[4] = 4.0*Point(0)*Point(1);
    _N[5] = 4.0*Point(1)*Point(2);
    _N[6] = 4.0*Point(2)*Point(0);
    _N[7] = 4.0*Point(0)*csi;
    _N[8] = 4.0*Point(1)*csi;
    _N[9] = 4.0*Point(2)*csi;

    // Shape functions derivatives
    _DN[0](0) = 4.0*Point(0) - 1.0;
    _DN[0](1) = 0.0;
    _DN[0](2) = 0.0;

    _DN[1](0) = 0.0;
    _DN[1](1) = 4.0*Point(1) - 1.0;
    _DN[1](2) = 0.0;

    _DN[2](0) = 0.0;
    _DN[2](1) = 0.0;
    _DN[2](2) = 4.0*Point(2) - 1.0;

    _DN[3](0) = 1.0 - 4.0*csi;
    _DN[3](1) = 1.0 - 4.0*csi;
    _DN[3](2) = 1.0 - 4.0*csi;

    _DN[4](0) = 4.0*Point(1);
    _DN[4](1) = 4.0*Point(0);
    _DN[4](2) = 0.0;

    _DN[5](0) = 0.0;
    _DN[5](1) = 4.0*Point(2);
    _DN[5](2) = 4.0*Point(1);

    _DN[6](0) = 4.0*Point(2);
    _DN[6](1) = 0.0;
    _DN[6](2) = 4.0*Point(0);

    _DN[7](0) = 4.0*(csi - Point(0));
    _DN[7](1) =-4.0*Point(0);
    _DN[7](2) =-4.0*Point(0);

    _DN[8](0) =-4.0*Point(1);
    _DN[8](1) = 4.0*(csi - Point(1));
    _DN[8](2) =-4.0*Point(1);

    _DN[9](0) =-4.0*Point(2);
    _DN[9](1) =-4.0*Point(2);
    _DN[9](2) = 4.0*(csi - Point(2));

  } // update
						       
} // namespace voom
