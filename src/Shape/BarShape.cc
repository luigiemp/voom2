#include "BarShape.h"

namespace voom{

  void BarShape::update(const VectorXd & Point)
  {
    _N[0] = 0.5*(1.0 - Point(0));
    _N[1] = 0.5*(1.0 + Point(0));

    _DN[0] =-0.5;
    _DN[1] = 0.5;
  }

} // namespace voom
