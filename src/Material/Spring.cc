#include "Spring.h"
#include "VoomMath.h"

namespace voom {

  // Operators
  void Spring::compute(Filresults & R,  vector<Vector3d> & x, const  vector<Vector3d> & X)
  {

    Vector3d nodeA = x[0];
    Vector3d nodeB = x[1];

    Vector3d nodeA0 = X[0];
    Vector3d nodeB0 = X[1];
    
    Vector3d disp = nodeB - nodeA;
    Vector3d disp0 = nodeB0 - nodeA0;

    if( R.request & ENERGY )
    {
      R.W = 0.5*_k*((disp-disp0).squaredNorm());
    }
    
    if( R.request & FORCE ) 
    {
      Real dl = (disp-disp0).norm();
      R.f1 << _k*dl*disp(0)/disp.norm() , _k*dl*disp(1)/disp.norm() , _k*dl*disp(2)/disp.norm();
    }
    if( R.request & STIFFNESS )
    {
      R.k = _k;
    } 
    
  } // Spring::compute


} // namespace voom

