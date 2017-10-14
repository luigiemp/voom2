#include "LinYinActive.h"

namespace voom {

  // Operators
  void LinYinActive::compute(FKresults & R, const Matrix3d & F)
  {
    // Needed for all requests
    Matrix3d C, invF, Delta;
    Delta << 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0;
    C = F.transpose()*F;
    invF = F.inverse();
    Real I1 = C.trace();
    Real I4 = _fibers[0].transpose() * C * (_fibers[0]);
    Real sI4 = sqrt(I4);
    Real detF = F.determinant();
 
    Matrix3d FNN = F * ((_fibers[0]) * _fibers[0].transpose());

    Matrix3d dI1dF = 2*F;
    Matrix3d dI4dF = 2*FNN;


    if( R.request & ENERGY )
    {
      R.W = _C0 + _C1*(I1-3.0)*(I4-1) + _C2*pow((I1-3),2.0) + _C3*pow((I4-1),2.0) + _C4*(I1-3) + _C5*(I4-1);
    }

    if(R.request & FORCE)
    {
      R.P = _C1*((I1-3)*dI4dF + dI1dF*(I4-1)) + _C2*(2*(I1-3)*dI1dF)+ _C3*(2*(I4-1)*dI4dF) + _C4*(dI1dF) +  _C5*(dI4dF);
    }

    if( R.request & STIFFNESS )
    {
      R.K.resetIterator();
      for (unsigned int L = 0; L<3; L++) {
        for (unsigned int k = 0; k<3; k++) {
          for (unsigned int J = 0; J<3; J++) {
            for (unsigned int i = 0; i<3; i++) {
              R.K.sequentialSet( 2.0 * ( _C1*(I4-1) + 2.0*_C2*(I1-3) + _C4 )* ( Delta(i,k)*Delta(J,L) ) +
				 2.0*( _C1*(I1-3) +  2.0*_C3*(I4-1) + _C5)*(_fibers[0][L]*_fibers[0][J]*Delta(i,k)) +
				 _C1*( dI1dF(i,J)*dI4dF(k,L) + dI1dF(k,L)*dI4dF(i,J) ) + 
				 2.0*_C2*dI1dF(i,J)*dI1dF(k,L) +
				 2.0*_C3*(dI4dF(k,L)*dI4dF(i,J)) );
              R.K.incrementIterator();
            } // L
          } // k
        } // J
      } // i

    } // STIFFNESS
  } // LinYinActive::compute
} // namespace voom
