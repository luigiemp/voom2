#include "LinYinActive_Compressible.h"

namespace voom {

  // Operators
  void LinYinActive_Compressible::compute(FKresults & R, const Matrix3d & F, Vector3d * Fiber)
  {
    // Needed for all requests
    Matrix3d C = F.transpose() * F;
    Real I1 = C.trace();
    Real I4 = Fiber->transpose() * C * (*Fiber);
    Real sI4 = sqrt(I4);
    Matrix3d FNN = F * ((*Fiber) * Fiber->transpose());
    Real detF = F.determinant();
    Matrix3d Finv = F.inverse();

    Matrix3d dI1dF = 2*F;
    Matrix3d dI4dF = 2*FNN;

    Matrix3d invF, ID;
    ID << 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0;

    if( R.request & ENERGY )
    {
      R.W =  _C1*(I1-3.0)*(I4-1) + _C2*pow((I1-3),2.0) + _C3*pow((I4-1),2.0) + _C4*(I1 - 3.0 - 2.0 * log(detF)) + _C5 * (0.5 * pow(log(detF),2));
    }

    if(R.request & FORCE)
    {
      R.P = _C1*((I1-3)*dI4dF + dI1dF*(I4-1)) + _C2*(2*(I1-3)*dI1dF)+ _C3*(2*(I4-1)*dI4dF) + _C4*(dI1dF - 2 * Finv.transpose()) + _C5 * (log(detF) * Finv.transpose());
    }

    if( R.request & STIFFNESS )
    {
      R.K.resetIterator();
      for (unsigned int L = 0; L<3; L++) {
        for (unsigned int k = 0; k<3; k++) {
          for (unsigned int J = 0; J<3; J++) {
            for (unsigned int i = 0; i<3; i++) {
              R.K.sequentialSet(_C1*(2*ID(i,k)*ID(J,L)*(I4-1)+ dI1dF(i,J)*dI4dF(k,L)+dI1dF(k,L)*dI4dF(i,J) +
                  (I1-3)*2*(*Fiber)[L]*(*Fiber)[J]*ID(i,k)) + 2*_C2*(dI1dF(k,L)*dI1dF(i,J)+(I1-3)*2*ID(i,k)*ID(J,L)) +
                  2*_C3*(dI4dF(k,L)*dI4dF(i,J)+(I4-1)*2*(*Fiber)[L]*(*Fiber)[J]*ID(i,k)) +
                  _C4 * 2 * (ID(i,k) * ID(J,L) + Finv(J,k) * Finv(L,i)) + 
		  _C5 * (Finv(L,k) * Finv(J,i) - log(detF) * Finv(J,k) * Finv(L,i)));
              R.K.incrementIterator();
            } // L
          } // k
        } // J
      } // i

    } // STIFFNESS
  } // LinYinActive::compute
} // namespace voom
