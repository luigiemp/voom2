#include "ImposedKinematics.h"

namespace voom {

  // Operators
  void ImposedKinematics::compute(FKresults & R, const Matrix3d & F)
  {
    // Needed for all requests
    Matrix3d C, invF, Delta;
    Delta << 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0;
    C = F.transpose()*F;
    invF = F.inverse();
    Real I3 = C.determinant(); // cout << "I3 = " << I3 << endl;


    // Start from incompressibility term
    if( R.request & ENERGY ) {
      R.W = _beta*(pow(I3, 2.0) + pow(I3,-2.0) - 2.0) +
            _gamma*( C.trace() - 3.0 ); 
    }
    if( (R.request & FORCE) || (R.request & DMATPROP) ) {
      R.P = 4.0*_beta*(I3 - pow(I3,-3.0) )*I3*invF.transpose() +
	2.0*_gamma*F;
    }
    if( R.request & STIFFNESS )
    {
      R.K.resetIterator();
      for (unsigned int L = 0; L<3; L++) {
	for (unsigned int k = 0; k<3; k++) {
	  for (unsigned int J = 0; J<3; J++) {
	    for (unsigned int i = 0; i<3; i++) {
	      R.K.sequentialSet( 2.0*_beta*( (I3 - pow(I3,-3.0) )*2.0*I3*( 2.0*invF(J,i)*invF(L,k) - invF(J,k)*invF(L,i) ) + 4.0*(1.0 + 3.0*pow(I3,-4.0) )*pow(I3,2.0)*invF(J,i)*invF(L,k) ) + 2.0*_gamma*Delta(i,k)*Delta(J,L) );
	      R.K.incrementIterator();
	    } // L
	  } // k
	} // J
      } // i
    } // STIFFNESS



    // Initialize derivatives with respect to materil properties
    if( R.request & DMATPROP ) 
    {
      int NumMat = _alphas.size();
      R.Dmat.resize(NumMat, 3, 3);          // Already initialized to zero
      R.DDmat.resize(NumMat, NumMat, 3, 3); // Already initialized to zero 
    } // DMATPROP



    // Loop through all the energy terms of the form alpha*(v^T C v - \lambda)^2
    for (int a = 0; a < _alphas.size(); a++) {
      Vector3d v = _directions[a];
      Real Is = v.dot(C*v), alpha = _alphas[a], lambda = _stretches[a];
      Matrix3d FM = F*(v*v.transpose());
      Matrix3d Ptemp; 
      
       if( R.request & ENERGY ) {
	 R.W += alpha*pow(Is - lambda, 2.0); 
       }
       
       if( (R.request & FORCE) || (R.request & DMATPROP) ) {
	 Ptemp = 4.0*alpha*(Is - lambda)*FM;
	 R.P  += Ptemp;
       }
      
      if( R.request & STIFFNESS ) {
	for (unsigned int L = 0; L<3; L++) {
	  for (unsigned int k = 0; k<3; k++) {
	    for (unsigned int J = 0; J<3; J++) {
	      for (unsigned int i = 0; i<3; i++) {
		R.K.add(i, J, k, L, 4.0*alpha*(2.0*FM(i,J)*FM(k,L) + (Is - lambda)*Delta(i,k)*(v(L))*(v(J)) ) );
	      } // L
	    } // k
	  } // J
	} // i
      } // STIFFNESS
      
      if( R.request & DMATPROP ) {
	for (unsigned int i = 0; i<3; i++) {
	  for (unsigned int J = 0; J<3; J++) {
	    (R.Dmat).set(a, i, J, Ptemp(i,J)/alpha );
	  }
	}
      } // DMATPROP

    } // Loop through all the energy terms

    

  } // ImposedKinematics::compute

    
} // namespace voom

