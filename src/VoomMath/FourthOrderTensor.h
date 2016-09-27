//-*-C++-*-
#ifndef __Fourth_Order_Tensor_h__
#define __Fourth_Order_Tensor_h__

#include "VoomMath.h"

namespace voom
{
  class FourthOrderTensor
  {
    
  public:
    FourthOrderTensor():
      _sizeI(0), _sizeJ(0), _sizeK(0), _sizeL(0)
    { 
      CurrentPos = A.begin();
    };
      
    FourthOrderTensor(uint sizeI, uint sizeJ, uint sizeK, uint sizeL): 
      _sizeI(sizeI), _sizeJ(sizeJ), _sizeK(sizeK), _sizeL(sizeL),
      A(sizeI*sizeJ*sizeK*sizeL, 0.0)  
    { 
      CurrentPos = A.begin();
    };

    double& operator()(uint i, uint j, uint k, uint l)
    {
      return A[i + _sizeI*(j + _sizeJ*(k + l*_sizeK))];
    }

    double operator()(uint i, uint j, uint k, uint l) const
    {
      return A[i + _sizeI*(j + _sizeJ*(k + l*_sizeK))];
    }
    
    void resize(uint sizeI, uint sizeJ, uint sizeK, uint sizeL)
    {
      _sizeI = sizeI; 
      _sizeJ = sizeJ;
      _sizeK = sizeK;
      _sizeL = sizeL;
      A.assign(sizeI*sizeJ*sizeK*sizeL, 0.0);   
      CurrentPos = A.begin();
    };

    void set(uint i, uint j, uint k, uint l, Real a) {
      A[i + _sizeI*(j + _sizeJ*(k + l*_sizeK))] = a;
    };

    Real get(uint i, uint j, uint k, uint l)  {
      return A[i + _sizeI*(j + _sizeJ*(k + l*_sizeK))];
    };

    // Can be used only after setting CurrentPos
    void sequentialSet(Real a) { *CurrentPos = a; };
    Real sequentialGet() { return *CurrentPos; };

    void incrementIterator() { CurrentPos++; };
    void resetIterator() { CurrentPos = A.begin(); };

  private:
    vector<Real >::iterator CurrentPos;
    uint _sizeI;
    uint _sizeJ;
    uint _sizeK;
    uint _sizeL;
    vector<Real > A;
    
  }; // class FourthOrderTensor

} // namespace voom

#endif
