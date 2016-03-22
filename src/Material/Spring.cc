#include "Spring.h"

namespace voom {

  // Operators
  void Spring::compute(Filresults & R, const Vector3d & d)
  {
    
        
    if( R.request & ENERGY )
    {
      //R.W = 0.5*_k*sqr(norm2(d-_d0));
      R.W = 0.0;
    }
    
    if( R.request & FORCE ) 
    {
      R.f << 0.0,0.0,0.0;
    }
    if( R.request & STIFFNESS )
    {
      R.k = 0.0;
    } // i
    
  } // Spring::compute


} // namespace voom

