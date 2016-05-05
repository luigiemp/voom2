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

    Real d = disp.norm();
    Real d0 = disp0.norm();

    if( R.request & ENERGY )
    {
      R.W = 0.5*_k*pow((d-d0),2);
    }
    
    if( R.request & FORCE ) 
    {
      Real dl = (disp-disp0).norm();
      R.f1 << _k*(d-d0)*(nodeA(0)-nodeB(0))/d , _k*(d-d0)*(nodeA(1)-nodeB(1))/d , _k*(d-d0)*(nodeA(2)-nodeB(2))/d;
    }
    if( R.request & STIFFNESS )
    {
      R.k = _k;
    } 
    
  } // Spring::compute


   // Operators
  void Spring::compute(Filresults & R,  vector<Vector3d> & x)
  {

    Vector3d nodeA = x[0];
    Vector3d nodeB = x[1];

    Vector3d disp = nodeB - nodeA;
    
    Real d = disp.norm();
    Real d0 = _L0;

    if( R.request & ENERGY )
    {
      R.W = 0.5*_k*pow((d-d0),2);
    }
    
    if( R.request & FORCE ) 
    {
      R.f1 << _k*(d-d0)*(nodeA(0)-nodeB(0))/d , _k*(d-d0)*(nodeA(1)-nodeB(1))/d , _k*(d-d0)*(nodeA(2)-nodeB(2))/d;
    }
    if( R.request & STIFFNESS )
    {
      R.k = _k;
    } 
    
  } // Spring::compute


} // namespace voom

