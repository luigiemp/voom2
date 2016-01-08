#include "Jacobian.h"

namespace voom {

  // Operators
  void Jacobian::compute(FKresults & R, const Matrix3d & F, Vector3d * Fiber)
  {
    // Needed for all requests
    Real DetF = F.determinant();
    Matrix3d invF;

    if( (R.request & FORCE) || (R.request & STIFFNESS) )
    {
      invF = F.inverse();
    }
    
    if( R.request & ENERGY )
    {
      R.W = DetF;
    }
    if( R.request & FORCE ) 
    {
      R.P = DetF*invF.transpose();
    }
    if( R.request & STIFFNESS )
    {
      R.K.resetIterator();
      for (unsigned int L = 0; L<3; L++) {
	for (unsigned int k = 0; k<3; k++) {
	  for (unsigned int J = 0; J<3; J++) {
	    for (unsigned int i = 0; i<3; i++) {
	      // R.K.set(i,J,k,L, -invF(J,k)*invF(L,i)*(_lambda*LogDetF - _mu) + _lambda*invF(J,i)*invF(L,k) + _mu*ID(i,k)*ID(J,L) ); 
	      R.K.sequentialSet( DetF*( invF(J,i)*invF(L,k) - invF(J,k)*invF(L,i) ) );
	      R.K.incrementIterator();
	    } // L
	  } // k
	} // J
      } // i
    
    } // STIFFNESS

    if( R.request & DMATPROP ) 
    {
      R.Dmat.resize(0, 0, 0);     // Already initialized to zero
      R.DDmat.resize(0, 0, 0, 0); // Already initialized to zero 
    } // DMATPROP

  } // Jacobian::compute


} // namespace voom

