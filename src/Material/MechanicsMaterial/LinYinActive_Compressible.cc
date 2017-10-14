#include "LinYinActive_Compressible.h"

namespace voom {

  // Operators
  void LinYinActive_Compressible::compute(FKresults & R, const Matrix3d & F)
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
      R.W =  _C1*(I1bar-3.0)*(I4-1) + _C2*pow((I1bar-3),2.0) + _C3*pow((I4-1),2.0) + _C4*(I1bar - 3.0 - 2.0 * log(detF)) + _C5 * (0.5 * pow(log(detF),2));
    }

    if(R.request & FORCE)
    {
      R.P = _C1*((I1bar-3)*dI4dF + dI1bardF*(I4-1)) + _C2*(2*(I1bar-3)*dI1bardF)+ _C3*(2*(I4-1)*dI4dF) + _C4*(dI1bardF - 2 * invF.transpose()) + _C5 * (log(detF) * invF.transpose());
    }

    if( R.request & STIFFNESS )
    {
      R.K.resetIterator();
      for (unsigned int L = 0; L<3; L++) {
        for (unsigned int k = 0; k<3; k++) {
          for (unsigned int J = 0; J<3; J++) {
            for (unsigned int i = 0; i<3; i++) {
              R.K.sequentialSet( 2.0 * ( _C1*(I4-1) + 2.0*_C2*(I1bar-3) + _C4 )*I3Third*( (-2.0/3.0)*invF(L,k)*F(i,J) + Delta(i,k)*Delta(J,L) - (2.0/3.0)*F(k,L)*invF(J,i) + (2.0/9.0)*I1*invF(L,k)*invF(J,i) + (I1/3.0)*invF(J,k)*invF(L,i) ) +
				 2.0*( _C1*(I1bar-3) +  2.0*_C3*(I4-1) )*(_fibers[0][L]*_fibers[0][J]*Delta(i,k)) +
				 _C1*( dI1bardF(i,J)*dI4dF(k,L) + dI1bardF(k,L)*dI4dF(i,J) ) + 
				 2.0*_C2*dI1bardF(i,J)*dI1bardF(k,L) +
				 2.0*_C3*(dI4dF(k,L)*dI4dF(i,J)) +
				 2.0*_C4*invF(J,k) * invF(L,i) +
				 _C5 * (invF(L,k) * invF(J,i) - log(detF) * invF(J,k) * invF(L,i)) );
              R.K.incrementIterator();
            } // L
          } // k
        } // J
      } // i

    } // STIFFNESS
  } // LinYinActive::compute
} // namespace voom
