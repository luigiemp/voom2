#include "Humphrey.h"

namespace voom {

  // Operators
  void Humphrey::compute(FKresults & R, const Matrix3d & F)
  {
    // Needed for all requests
    Matrix3d C = F.transpose() * F;
    Real I1 = C.trace();
    Real I4 = _fibers[0].transpose() * C * (_fibers[0]);
    Real sI4 = sqrt(I4);
    Matrix3d FNN = F * ((_fibers[0]) * _fibers[0].transpose());

    Matrix3d invF, ID;
    ID << 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0;

    if( R.request & ENERGY )
    {
      R.W = _C1 * pow((sI4 - 1.0),2.0) + _C2 * pow((sI4 - 1.0),3.0) + _C3 * (I1 - 3.0) + _C4 * (I1 - 3.0) * (sI4 - 1.0) + _C5 * pow((I1-3.0),2.0);
    }

    if(R.request & FORCE)
    {
      R.P =  2.0*(_C3 + _C4*(sI4-1) + 2.0*_C5*(I1-3)) * F + (1.0/sI4) * ((sI4-1.0) * (2.0*_C1 + 3.0 * _C2 * (sI4-1.0)) + _C4 * (I1-3.0))*FNN;
    }

    if( R.request & STIFFNESS )
    {
      R.K.resetIterator();
      for (unsigned int L = 0; L<3; L++) {
        for (unsigned int k = 0; k<3; k++) {
          for (unsigned int J = 0; J<3; J++) {
            for (unsigned int i = 0; i<3; i++) {
              R.K.sequentialSet(2.0*(_C3 + _C4*(sI4-1.0) + 2.0*_C5*(I1-3.0))*ID(i,k)*ID(J,L) +
                             ((1.0-1.0/sI4)*(2.0*_C1 + 3.0*_C2*(sI4-1.0)) + _C4*(I1-3.0)/sI4)* ((_fibers[0])[L])*((_fibers[0])[J])*ID(i,k) +
                             ((2.0*_C4/sI4)*FNN(k,L) + 8.0*_C5*F(k,L))*F(i,J) +
                             FNN(i,J)*(FNN(k,L)*((2.0*_C1 + 3.0*_C2*(sI4-1.0))*pow(sI4,(-3.0)) + (1.0-1.0/sI4)*(3.0*_C2/sI4) -
                                       _C4*(I1-3.0)*pow(sI4,(-3.0))) + 2.0*_C4*F(k,L)/sI4));
              // R.K.set(i,J,k,L, -invF(J,k)*invF(L,i)*(_lambda*LogDetF - _mu) + _lambda*invF(J,i)*invF(L,k) + _mu*ID(i,k)*ID(J,L) );
              // R.K.sequentialSet(-invF(J,k)*invF(L,i)*(_lambda*LogDetF - _mu) + _lambda*invF(J,i)*invF(L,k) + _mu*ID(i,k)*ID(J,L) );
              R.K.incrementIterator();
            } // L
          } // k
        } // J
      } // i

    } // STIFFNESS
  } // CompNeoHookean::compute


} // namespace voom
