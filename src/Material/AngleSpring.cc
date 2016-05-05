#include "AngleSpring.h"
#include "VoomMath.h"

namespace voom {


  void AngleSpring::compute(Filresults & R, vector< Vector3d> & x,const vector<Vector3d> & X)
  {
    Vector3d nodeA = x[2];
    Vector3d nodeB = x[0];
    Vector3d nodeC = x[1];
    
    // Use X is filament not straight in ref position

    Vector3d tBA;
    tBA = nodeA-nodeB;
    double LBA = tBA.norm();
    tBA /= LBA;

    Vector3d tBC;
    tBC = nodeC-nodeB;
    double LBC = tBC.norm();
    tBC /= LBC;

    double cosABC = tBA.dot(tBC);

    if( R.request & ENERGY )
      {
	R.W = _kappa*(1.0+cosABC);
      }

    if( R.request & FORCE )
      {
	Vector3d fA;
	fA = -cosABC * tBA;
	fA += tBC;
	fA *=  _kappa/LBA;

	Vector3d fC;
	fC = -cosABC * tBC;
	fC += tBA;
	fC *=  _kappa/LBC;

	R.f1 = fA;
	R.f2 = fC;
      }
    if( R.request & STIFFNESS )
      {
	R.k = _kappa;
      }

    
  }// end compute
  

  void AngleSpring::compute(Filresults & R,  vector<Vector3d> & x){
  };

} // namespace voom

