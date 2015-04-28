//-*-C++-*-
/*!
  \file LoopShellElement.h

  \brief  
*/

#if !defined(__LoopShellElement_h__)
#define __LoopShellElement_h__

#include "FEgeomElement.h"

namespace voom {

  class LoopShellElement: public FEgeomElement{
  public:

    LoopShellElement(const int elemID, const vector<int > & nodesID, 
		     const vector<VectorXd > & nodesX,
		  vector<Shape* > shape, Quadrature* quadrature);

    //! Get second derivatives of shape functions at quadrature point
    //! q, node a, direction (i,j)
    Real getDDN(uint q, uint a, uint i, uint j) {
      return _DDN[q*_nodesID.size() + a](i,j); 
    }

  protected:
    
    vector<Matrix2d> _DDN;


  }; // LoopShellElement

} // namespace voom

#endif
