#include "Humphrey_Compressible.h"

namespace voom {

  // Operators
  void Humphrey_Compressible::compute(FKresults & R, const Matrix3d & F)
  {
    // Needed for all requests
    Matrix3d C, invF, Delta;
    Delta << 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0;
    C = F.transpose()*F;
    invF = F.inverse();
    Real I3 = C.determinant();
    Real I3Third = pow(I3, -1.0/3.0);
    Real I1 = C.trace();
    Real I1bar = I1*I3Third;
    Real I4 = _fibers[0].transpose() * C * (_fibers[0]);
    Real sI4 = sqrt(I4);
    Matrix3d FNN = F * ((_fibers[0]) * _fibers[0].transpose());
    Real detF = F.determinant();
    Matrix3d dI1bardF = 2*I3Third*(F - (I1/3.0)*invF.transpose());
    Matrix3d dI4dF = 2*FNN;



    if( R.request & ENERGY )
    {
      R.W = _C1 * pow((sI4 - 1.0),2.0) + _C2 * pow((sI4 - 1.0),3.0) + _C3 * (I1bar - 3.0) + _C4 * (I1bar - 3.0) * (sI4 - 1.0) + _C5 * pow((I1bar-3.0),2.0) + _C6*(I1bar - 3.0 - 2.0 * log(detF)) + _C7 * (0.5 * pow(log(detF),2));
    }

    if(R.request & FORCE)
    {
      R.P = (_C3 + _C4*(sI4-1) + 2.0*_C5*(I1bar-3.0)) * dI1bardF + (1.0/sI4) * ((sI4-1.0) * (2.0*_C1 + 3.0 * _C2 * (sI4-1.0)) + _C4 * (I1bar-3.0))*FNN + _C6*(dI1bardF - 2.0 * invF.transpose()) + _C7 * (log(detF) * invF.transpose());

      // R.P = 2.0 * _C1 *(1.0/sI4)*(sI4-1.0)*FNN + 
      // 	3.0 * _C2 *(1.0/sI4)*pow(sI4-1.0, 2.0)*FNN + 
      //         _C3 *dI1bardF     +
      // 	      _C4 *( (sI4-1)* dI1bardF + (1.0/sI4)*(I1bar-3.0)*FNN) + 
      // 	2.0 * _C5 *(I1bar-3.0) * dI1bardF +
      // 	      _C6 *(dI1bardF - 2.0 * invF.transpose()) +
      // 	      _C7 *(log(detF) * invF.transpose());
    }

    if( R.request & STIFFNESS )
    {
      R.K.resetIterator();
      for (unsigned int L = 0; L<3; L++) {
        for (unsigned int k = 0; k<3; k++) {
          for (unsigned int J = 0; J<3; J++) {
            for (unsigned int i = 0; i<3; i++) {
	      R.K.sequentialSet( _fibers[0][L]*_fibers[0][J]*Delta(i,k)*( (1.0-1.0/sI4)*(2.0*_C1 + 3.0*_C2*(sI4-1.0)) + (1.0/sI4)*(I1bar-3.0)*_C4) + 
				 2.0 * _C1 * ( pow(sI4, -3.0)*FNN(i,J)*FNN(k,L) ) +
				 3.0 * _C2 * ( (sI4-1.0)*pow(sI4, -3.0) + (1.0-1.0/sI4)*(1.0/sI4) )*FNN(i,J)*FNN(k,L) +
				 2.0 * (_C3 + _C4*(sI4-1.0) + _C6 + 2.0*_C5*(I1bar-3.0))*I3Third*( (-2.0/3.0)*invF(L,k)*F(i,J) + Delta(i,k)*Delta(J,L) - (2.0/3.0)*F(k,L)*invF(J,i) + (2.0/9.0)*I1*invF(L,k)*invF(J,i) + (I1/3.0)*invF(J,k)*invF(L,i) ) +
				 _C4 * ( - (I1bar-3.0)*pow(sI4, -3.0)*FNN(i,J)*FNN(k,L) + (1.0/sI4)*FNN(i,J)*dI1bardF(k,L) + dI1bardF(i,J)*FNN(k,L)/sI4 ) +
				 2.0 * _C5 *( dI1bardF(i,J) * dI1bardF(k,L) ) + 
				 (invF(J,k) * invF(L,i))*(2.0*_C6 - log(detF)*_C7)  +
				 _C7 * invF(L,k) * invF(J,i) );

              // R.K.sequentialSet( 2.0 * _C1 * ((1.0-1.0/sI4)*(_fibers[0][L])*(_fibers[0][J])*Delta(i,k) + pow(sI4, -3.0)*FNN(i,J)*FNN(k,L) ) +
	      // 			 3.0 * _C2 * ((1.0-1.0/sI4)*(sI4-1.0)*(_fibers[0][L])*(_fibers[0][J])*Delta(i,k) + ( (sI4-1.0)*pow(sI4, -3.0) + (1.0-1.0/sI4)*(1.0/sI4) )*FNN(i,J)*FNN(k,L) ) +
	      // 			 2.0 * (_C3 + _C4*(sI4-1.0) + _C6 + 2.0*_C5*(I1bar-3.0))*I3Third*( (-2.0/3.0)*invF(L,k)*F(i,J) + Delta(i,k)*Delta(J,L) - (2.0/3.0)*F(k,L)*invF(J,i) + (2.0/9.0)*I1*invF(L,k)*invF(J,i) + (I1/3.0)*invF(J,k)*invF(L,i) ) +
	      // 			 _C4 * ( (1.0/sI4)*(I1bar-3.0)*_fibers[0][L]*_fibers[0][J]*Delta(i,k) - (I1bar-3.0)*pow(sI4, -3.0)*FNN(i,J)*FNN(k,L) + 
	      // 				 (1.0/sI4)*FNN(i,J)*dI1bardF(k,L) + 
	      // 				 dI1bardF(i,J)*FNN(k,L)/sI4 ) +
	      // 			 2.0 * _C5 *( dI1bardF(i,J) * dI1bardF(k,L) ) + 
	      // 			 (invF(J,k) * invF(L,i))*(2.0*_C6 - log(detF)*_C7)  +
	      // 			 _C7 * invF(L,k) * invF(J,i) );

              R.K.incrementIterator();
            } // L
          } // k
        } // J
      } // i

    } // STIFFNESS

    if( R.request & DMATPROP ) 
    {
      R.Dmat.resize(7, 3, 3);     // Already initialized to zero
      R.DDmat.resize(7, 7, 3, 3); // Already initialized to zero 
      
      Matrix3d Pa = (1.0/sI4) * (sI4-1.0) * 2.0 * FNN; 
      Matrix3d Pb = (1.0/sI4) * (sI4-1.0) * 3.0 * (sI4-1.0) * FNN;
      Matrix3d Pc = dI1bardF; 
      Matrix3d Pd = (sI4-1.0) * dI1bardF + (1.0/sI4) * (I1bar-3.0) * FNN;
      Matrix3d Pe = 2.0 * (I1bar-3.0) * dI1bardF;
      Matrix3d Pf = dI1bardF - 2.0 * invF.transpose();
      Matrix3d Pg = log(detF) * invF.transpose();

      for (unsigned int i = 0; i<3; i++) {
	for (unsigned int J = 0; J<3; J++) {
	  (R.Dmat).set( 0, i, J, Pa(i,J) );
	  (R.Dmat).set( 1, i, J, Pb(i,J) );
	  (R.Dmat).set( 2, i, J, Pc(i,J) );
	  (R.Dmat).set( 3, i, J, Pd(i,J) );
	  (R.Dmat).set( 4, i, J, Pe(i,J) );
	  (R.Dmat).set( 5, i, J, Pf(i,J) );
	  (R.Dmat).set( 6, i, J, Pg(i,J) );
	} // i
      } // J
     
    } // DMATPROP

  } // Humphrey_Compressible::compute


} // namespace voom
