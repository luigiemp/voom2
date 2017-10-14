#include "CompNeoHookean.h"

namespace voom {

  // Copy constructor
  //CompNeoHookean::CompNeoHookean(const CompNeoHookean &Input) :  _data(Input._data)
  //{;}


  // Operators
  void CompNeoHookean::compute(FKresults & R, const Matrix3d & F)
  {
    
    Matrix3d C, invF, Delta;
    Delta << 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0;
    C = F.transpose()*F;
    invF = F.inverse();
    Real I1 = C.trace();
    Real I3 = C.determinant(); // cout << "I3 = " << I3 << endl;
    Real I3Third = pow(I3, -1.0/3.0);
    Real LogDetF = log(F.determinant());



    if( R.request & ENERGY )
    {
      R.W = LogDetF*(0.5*_lambda*LogDetF - _mu) + 0.5*_mu*(  I1*I3Third - 3.0 );
      // cout << LogDetF*(0.5*_lambda*LogDetF - _mu) << " " <<  0.5*_mu*( (F.transpose()*F).trace() - 3.0) << endl;
    }
    if( (R.request & FORCE) || (R.request & DMATPROP) ) 
    {
      R.P = (_lambda*LogDetF - _mu)*(invF.transpose()) + _mu*I3Third*(F - (I1/3.0)*invF.transpose());
    }
    if( R.request & STIFFNESS )
    {
      R.K.resetIterator();
      for (unsigned int L = 0; L<3; L++) {
	for (unsigned int k = 0; k<3; k++) {
	  for (unsigned int J = 0; J<3; J++) {
	    for (unsigned int i = 0; i<3; i++) {
	      // R.K.set(i,J,k,L, -invF(J,k)*invF(L,i)*(_lambda*LogDetF - _mu) + _lambda*invF(J,i)*invF(L,k) + _mu*ID(i,k)*ID(J,L) ); 
	      R.K.sequentialSet(-invF(J,k)*invF(L,i)*(_lambda*LogDetF - _mu) + _lambda*invF(J,i)*invF(L,k) + 
				_mu*I3Third*( (-2.0/3.0)*invF(L,k)*F(i,J) + Delta(i,k)*Delta(J,L) - (2.0/3.0)*F(k,L)*invF(J,i) + (2.0/9.0)*I1*invF(L,k)*invF(J,i) 
					      + (I1/3.0)*invF(J,k)*invF(L,i) ) );
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
      Matrix3d Pb = I3Third*F - (I3Third*I1/3.0 + 1.0)*invF.transpose();

      for (unsigned int i = 0; i<3; i++) {
	for (unsigned int J = 0; J<3; J++) {
	  (R.Dmat).set( 0, i, J, Pa(i,J) );
	  (R.Dmat).set( 1, i, J, Pb(i,J) );
	  // cout << LogDetF << " " << Pa(i,J) << " " << Pb(i,J) << endl;
	} // i
      } // J
    } // DMATPROP

  } // CompNeoHookean::compute


} // namespace voom

