#include "Harmonic.h"

namespace voom {

  // Operators
  void Harmonic::compute(PairMresults & R, Real r)
  {
    R.W = 0.5  * _k * square(r - _r0);
    R.F = _k * (r - _r0);
    R.K = _k;
  } // Harmonic::compute

} // namespace voom

