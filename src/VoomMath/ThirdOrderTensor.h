//-*-C++-*-
#ifndef __Third_Order_Tensor_h__
#define __Third_Order_Tensor_h__

#include "VoomMath.h"

namespace voom
{
  class ThirdOrderTensor
  {
    
  public:
    ThirdOrderTensor():
      _sizeI(0), _sizeJ(0), _sizeK(0)
    { 
      CurrentPos = A.begin();
    };
      
    ThirdOrderTensor(uint sizeI, uint sizeJ, uint sizeK): 
      _sizeI(sizeI), _sizeJ(sizeJ), _sizeK(sizeK),
      A(sizeI*sizeJ*sizeK, 0.0)  
    { 
      CurrentPos = A.begin();
    };
    
    void resize(uint sizeI, uint sizeJ, uint sizeK)
    {
      _sizeI = sizeI; 
      _sizeJ = sizeJ;
      _sizeK = sizeK;
      A.assign(sizeI*sizeJ*sizeK, 0.0);   
      CurrentPos = A.begin();
    };

    void set(uint i, uint j, uint k, Real a) {
      A[i + _sizeI*(j + _sizeJ*k)] = a;
    };

    Real get(uint i, uint j, uint k)  {
      return A[i + _sizeI*(j + _sizeJ*k)];
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
    vector<Real > A;
    
  }; // class ThirdOrderTensor

} // namespace voom

#endif
