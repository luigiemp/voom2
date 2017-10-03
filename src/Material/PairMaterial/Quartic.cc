#include "Quartic.h"

namespace voom {

  // Operators
  void Quartic::compute(PairMresults & R, Real r)
  {
    Real DeltaSquare =  square(r - _r0);
    R.W = 0.25  * _k * square(DeltaSquare);
    R.F = _k * DeltaSquare * (r - _r0) ;
    R.K = 3.0 * _k * DeltaSquare;
  } // Quartic::compute

} // namespace voom

