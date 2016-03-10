#include "Spring.h"

namespace voom {

  // Operators
  void Spring::compute(Filresults & R, const Vector3D & d)
  {
    
        
    if( R.request & ENERGY )
    {
      R.W = 0.5*_k*sqr(norm2(d-_d0));
    }
    if( R.request & FORCE ) 
    {
      R.F = _k*(norm2(d-_d0))/(d());
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

