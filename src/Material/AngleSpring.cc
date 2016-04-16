#include "AngleSpring.h"
#include "VoomMath.h"

namespace voom {


  void AngleSpring::compute(Filresults & R,  Vector3d & d,const Vector3d & d0)
  {

  }
  /*
  // Operators
  void AngleSpring::compute(Filresults & R,  Vector3d nodeA, Vector3d nodeB, Vector3d nodeC)
  {
        
    Vector3D tBA;
    tBA = nodeA-nodeB;
    double LBA = norm2(tBA);
    tBA /= LBA;

    Vector3D tBC;
    tBC = nodeC-nodeB;
    double LBC = norm2(tBC);
    tBC /= LBC;

    double cosABC = dot(tBA,tBC);

    


    if( R.request & ENERGY )
    {
      R.W = _kappa*(1.0+cosABC);
    }
    
    if( R.request & FORCE ) 
    {
      Vector3D fA;
      fA = -cosABC * tBA;
      fA += tBC;
      fA *=  _kappa/LBA;

      Vector3D fC;
      fC = -cosABC * tBC;
      fC += tBA;
      fC *=  _kappa/LBC;

      R.f = -fA -fC;
    }
    if( R.request & STIFFNESS )
    {
      R.k = _kappa;
    } 
    
  } // Spring::compute
  */

} // namespace voom

