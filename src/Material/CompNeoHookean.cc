#include "CompNeoHookean.h"

namespace voom {

  // Copy constructor
  //CompNeoHookean::CompNeoHookean(const CompNeoHookean &Input) :  _data(Input._data)
  //{;}


  // Operators
  void CompNeoHookean::compute(FKresults & R, const Matrix3d & F)
  {
    // Needed for all requests
    Real LogDetF = log(F.determinant());
    Matrix3d invF, ID;
    ID << 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0;

    if( (R.request & FORCE) || (R.request & STIFFNESS) )
    {
      invF = F.inverse();
    }
    
    if( R.request & ENERGY )
    {
      R.W = LogDetF*(0.5*_lambda*LogDetF - _mu) + 0.5*_mu*( (F.transpose()*F).trace() - 3.0);
    }
    if( R.request & FORCE ) 
    {
      R.P = (_lambda*LogDetF - _mu)*(invF.transpose()) + _mu*F;
    }
    if( R.request & STIFFNESS )
    {
      R.K.resetIterator();
      for (unsigned int L = 0; L<3; L++) {
	for (unsigned int k = 0; k<3; k++) {
	  for (unsigned int J = 0; J<3; J++) {
	    for (unsigned int i = 0; i<3; i++) {
	      // R.K.set(i,J,k,L, -invF(J,k)*invF(L,i)*(_lambda*LogDetF - _mu) + _lambda*invF(J,i)*invF(L,k) + _mu*ID(i,k)*ID(J,L) ); 
	      R.K.sequentialSet(-invF(J,k)*invF(L,i)*(_lambda*LogDetF - _mu) + _lambda*invF(J,i)*invF(L,k) + _mu*ID(i,k)*ID(J,L) );
	      R.K.incrementIterator();
	    } // L
	  } // k
	} // J
      } // i
    
    } // STIFFNESS

    if( R.request & DMATPROP ) 
    {
      R.Dmat.resize(2, 3, 3);     // Already initialized to zero
      R.DDmat.resize(2, 2, 3, 3); // Already initialized to zero 
      Matrix3d Pa = LogDetF*(invF.transpose());
      Matrix3d Pb = F - (invF.transpose());

      for (unsigned int i = 0; i<3; i++) {
	for (unsigned int J = 0; J<3; J++) {
	  (R.Dmat).set( 0, i, J, Pa(i,J) );
	  (R.Dmat).set( 1, i, J, Pb(i,J) );
	} // i
      } // J
    } // DMATPROP

  } // CompNeoHookean::compute


} // namespace voom

