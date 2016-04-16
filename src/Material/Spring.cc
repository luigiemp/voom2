#include "Spring.h"
#include "VoomMath.h"

namespace voom {

  // Operators
  void Spring::compute(Filresults & R,  Vector3d & d, const  Vector3d & d0)
  {
    
        
    if( R.request & ENERGY )
    {
      R.W = 0.5*_k*((d-d0).squaredNorm());
    }
    
    if( R.request & FORCE ) 
    {
      Real dl = (d-d0).norm();
      R.f << _k*dl*d(0)/d.norm() , _k*dl*d(1)/d.norm() , _k*dl*d(2)/d.norm();
    }
    if( R.request & STIFFNESS )
    {
      R.k = _k;
    } 
    
  } // Spring::compute


} // namespace voom

