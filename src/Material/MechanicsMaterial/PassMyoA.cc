#include "PassMyoA.h"

namespace voom {

  // Operators
  void PassMyoA::compute(FKresults & R, const Matrix3d & F)
  {
    // Needed for all requests
    Vector3d f;
    f = _fibers[0];
    Matrix3d C, invF, FM, Delta;
    Delta << 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0;;
    C = F.transpose()*F;
    invF = F.inverse();
    FM = F*(f*f.transpose());

    Real I1 = C.trace();
    Real I3 = C.determinant(); // cout << "I3 = " << I3 << endl;
    Real I4 = f.dot(C*f);
    Real I3third = pow(I3, 1.0/3.0);
   
    Real I1bar = I1/I3third;
    if (I1bar < 3.0) { I1bar = 3.0; }; // Theoretically I1bar can never be less than 3.0 but small numerical error may oocur. If this happen, pow(I1bar-3, a1) return nan - see definition of pow in C++
    Real MCI4 = I4 - 1.0;
    if (MCI4 < 0.0) {
      MCI4 = 0.0;
    }
    // cout << "MCI4 = " << MCI4 << endl;

    if( R.request & ENERGY )
    {
      R.W = _alpha1*pow(I1bar - 3.0, _a1) + _alpha2*pow(MCI4, _a2) + _beta*(pow(I3, 2.0) + pow(I3,-2.0) - 2.0) + _gamma*_alpha1*(I1bar - 3.0);
      // cout <<  _alpha1*pow(I1bar - 3.0, _a1) << " " <<  _alpha2*(pow(I3, 2.0) + pow(I3,-2.0) - 2.0) << " " <<  _beta*pow(MCI4, _a2) << endl;
      // cout <<  _alpha1*pow(I1bar - 3.0, _a1) << " " << pow(I1bar - 3.0, _a1) << " " << I1bar - 3.0  << endl;
      // cout << _alpha1 << " " << _a1 << endl;
    }
    if( (R.request & FORCE) || (R.request & DMATPROP) ) 
    {
      R.P = 2.0*_a1*_alpha1*pow(I1bar - 3.0, _a1-1.0)*( F - (I1/3.0)*invF.transpose() )/I3third + 
      	2.0*_alpha2*_a2*pow(MCI4, _a2-1.0)*FM +
      	4.0*_beta*(I3 - pow(I3,-3.0) )*I3*invF.transpose() +
	2.0*_gamma*_alpha1*( F - (I1/3.0)*invF.transpose() )/I3third;
      
      // a1 = 1
      // R.P = 2.0*_a1*_alpha1*( F - (I1/3.0)*invF.transpose() )/I3third + 
      // 	2.0*_alpha2*_a2*pow(MCI4, _a2-1.0)*FM +
      // 	4.0*_beta*(I3 - pow(I3,-3.0) )*I3*invF.transpose();

      // a1 = 2
      // R.P = 2.0*_a1*_alpha1*(I1bar - 3.0)*( F - (I1/3.0)*invF.transpose() )/I3third + 
      // 	2.0*_alpha2*_a2*pow(MCI4, _a2-1.0)*FM +
      // 	4.0*_beta*(I3 - pow(I3,-3.0) )*I3*invF.transpose();
    }
    if( R.request & STIFFNESS )
    {
      R.K.resetIterator();
      for (unsigned int L = 0; L<3; L++) {
	for (unsigned int k = 0; k<3; k++) {
	  for (unsigned int J = 0; J<3; J++) {
	    for (unsigned int i = 0; i<3; i++) {
	      R.K.sequentialSet(
	      			2.0*_a1*_alpha1*( 
	      					 (_a1-1.0)*pow(I1bar - 3.0, _a1-2.0)*2.0*(F(k,L) - (I1/3.0)*invF(L,k) )*( F(i,J) - (I1/3.0)*invF(J,i) )/pow(I3third, 2.0) + 
	      					pow(I1bar - 3.0, _a1-1.0)*(1.0/I3third)*( Delta(i,k)*Delta(J,L) - (1.0/3.0)*(2.0*F(k,L)*invF(J,i) - I1*invF(J,k)*invF(L,i) + 2.0*(F(i,J)-(I1/3.0)*invF(J,i))*invF(L,k)) )
	      					  ) +
	      			2.0*_alpha2*_a2*(pow(MCI4, _a2-2.0)*(_a2-1.0)*FM(i,J)*2.0*FM(k,L) + pow(MCI4, _a2-1.0)*Delta(i,k)*(f(L))*(f(J)) ) +
	      			2.0*_beta*( (I3 - pow(I3,-3.0) )*2.0*I3*( 2.0*invF(J,i)*invF(L,k) - invF(J,k)*invF(L,i) ) + 4.0*(1.0 + 3.0*pow(I3,-4.0) )*pow(I3,2.0)*invF(J,i)*invF(L,k) ) +
				2.0*_gamma*_alpha1*( (1.0/I3third)*( Delta(i,k)*Delta(J,L) - (1.0/3.0)*(2.0*F(k,L)*invF(J,i) - I1*invF(J,k)*invF(L,i) + 2.0*(F(i,J)-(I1/3.0)*invF(J,i))*invF(L,k)) ) )
	      			);
	      // a1 = 1
	      // R.K.sequentialSet(
	      // 			2.0*_a1*_alpha1*( 
	      // 					 (1.0/I3third)*( Delta(i,k)*Delta(J,L) - (1.0/3.0)*(2.0*F(k,L)*invF(J,i) - I1*invF(J,k)*invF(L,i) + 2.0*(F(i,J)-(I1/3.0)*invF(J,i))*invF(L,k)) )
	      // 					  ) +
	      // 			2.0*_alpha2*_a2*(pow(MCI4, _a2-2.0)*(_a2-1.0)*FM(i,J)*2.0*FM(k,L) + pow(MCI4, _a2-1.0)*Delta(i,k)*_N(L)*_N(J) ) +
	      // 			2.0*_beta*( (I3 - pow(I3,-3.0) )*2.0*I3*( 2.0*invF(J,i)*invF(L,k) - invF(J,k)*invF(L,i) ) + 4.0*(1.0 + 3.0*pow(I3,-4.0) )*pow(I3,2.0)*invF(J,i)*invF(L,k) )
	      // 			);
	      // a1 = 2
	      // R.K.sequentialSet(
	      // 			2.0*_a1*_alpha1*( 
	      // 					 2.0*(F(k,L) - (I1/3.0)*invF(L,k) )*( F(i,J) - (I1/3.0)*invF(J,i) )/pow(I3third, 2.0) + 
	      // 					(I1bar - 3.0)*(1.0/I3third)*( Delta(i,k)*Delta(J,L) - (1.0/3.0)*(2.0*F(k,L)*invF(J,i) - I1*invF(J,k)*invF(L,i) + 2.0*(F(i,J)-(I1/3.0)*invF(J,i))*invF(L,k)) )
	      // 					  ) +
	      // 			2.0*_alpha2*_a2*(pow(MCI4, _a2-2.0)*(_a2-1.0)*FM(i,J)*2.0*FM(k,L) + pow(MCI4, _a2-1.0)*Delta(i,k)*_N(L)*_N(J) ) +
	      // 			2.0*_beta*( (I3 - pow(I3,-3.0) )*2.0*I3*( 2.0*invF(J,i)*invF(L,k) - invF(J,k)*invF(L,i) ) + 4.0*(1.0 + 3.0*pow(I3,-4.0) )*pow(I3,2.0)*invF(J,i)*invF(L,k) )
	      // 			);
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
      Matrix3d Pa = 2.0*_a1*pow(I1bar - 3.0, _a1-1.0)*( F - (I1/3.0)*invF.transpose() )/I3third + 2.0*_gamma*( F - (I1/3.0)*invF.transpose() )/I3third;
      Matrix3d Pb = 2.0*_a2*pow(MCI4, _a2-1.0)*FM;
      // Matrix3d Pc = 4.0*(I3 - pow(I3,-3.0) )*I3*invF.transpose();
      // Matrix3d Pd = 2.0*( F - (I1/3.0)*invF.transpose() )/I3third;
      
      for (unsigned int i = 0; i<3; i++) {
	for (unsigned int J = 0; J<3; J++) {
	  (R.Dmat).set( 0, i, J, Pa(i,J) );
	  (R.Dmat).set( 1, i, J, Pb(i,J) );
	  // (R.Dmat).set( 2, i, J, Pc(i,J) );
	  // (R.Dmat).set( 3, i, J, Pd(i,J) );
	} // i
      } // J
    } // DMATPROP

  } // PassMyoA::compute

    
} // namespace voom

