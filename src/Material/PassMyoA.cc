#include "PassMyoA.h"

namespace voom {

  // Operators
  void PassMyoA::compute(FKresults & R, const Matrix3d & F)
  {
    // Needed for all requests
    Matrix3d C, invF, FM, Delta;
    Delta << 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0;;
    C = F.transpose()*F;
    invF = F.inverse();
    FM = F*(_N*_N.transpose());

    Real I1 = C.trace();
    Real I3 = C.determinant();
    Real I4 = _N.dot(C*_N);
    Real I3third = pow(I3, 1.0/3.0);
    Real I1bar = I1/I3third;
    Real MCI4 = I4 - 1.0;
    if (MCI4 < 0.0) {
      MCI4 = 0.0;
    }

    if( R.request & ENERGY )
    {
      R.W = _alpha*pow(I1bar - 3.0, _a1) + _beta*pow(MCI4, _a2);
    }
    if( R.request & FORCE ) 
    {
      R.P = 2.0*_a1*_alpha*pow(I1bar - 3.0, _a1-1.0)*( F - (I1/3.0)*invF.transpose() )/I3third + 
	    2.0*_beta*_a2*pow(MCI4, _a2-1.0)*FM;
    }
    if( R.request & STIFFNESS )
    {
      R.K.resetIterator();
      for (unsigned int L = 0; L<3; L++) {
	for (unsigned int k = 0; k<3; k++) {
	  for (unsigned int J = 0; J<3; J++) {
	    for (unsigned int i = 0; i<3; i++) {
	      R.K.sequentialSet(
				2.0*_a1*_alpha*( 
						(_a1-1.0)*pow(I1bar - 3.0, _a1-2.0)*2.0*(F(k,L) - (I1/3.0)*invF(L,k) )*( F(i,J) - (I1/3.0)*invF(J,i) )/pow(I3third, 2.0) + 
						pow(I1bar - 3.0, _a1-1.0)*(1.0/I3third)*( Delta(i,k)*Delta(J,L) - (1.0/3.0)*(2.0*F(k,L)*invF(J,i) - I1*invF(J,k)*invF(L,i) + 2.0*(F(i,J)-(I1/3.0)*invF(J,i))*invF(L,k)) )
						 ) +
				2.0*_beta*_a2*(pow(MCI4, _a2-2.0)*(_a2-1.0)*FM(i,J)*2.0*FM(k,L) + pow(MCI4, _a2-1.0)*Delta(i,k)*_N(L)*_N(J) ) 
				);
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
      Matrix3d Pa = 2.0*_a1*pow(I1bar - 3.0, _a1-1.0)*( F - (I1/3.0)*invF.transpose() )/I3third;
      Matrix3d Pb = 2.0*_a2*pow(MCI4, _a2-1.0)*FM;

      for (unsigned int i = 0; i<3; i++) {
	for (unsigned int J = 0; J<3; J++) {
	  (R.Dmat).set( 0, i, J, Pa(i,J) );
	  (R.Dmat).set( 1, i, J, Pb(i,J) );
	} // i
      } // J
    } // DMATPROP

  } // PassMyoA::compute

    
} // namespace voom

