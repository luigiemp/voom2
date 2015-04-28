//-*-C++-*-

#include "LoopShellElement.h"
#include "LoopShellShape.h"

namespace voom{

  LoopShellElement::LoopShellElement(const int elemID, const vector<int > & nodesID, 
				     const vector<VectorXd > & nodesX,
				     vector<Shape* > shape, Quadrature* quadrature): 
    FEgeomElement(elemID, nodesID)
  {
    
    // Compute shape functions, shape functions derivatives and QP weights
    
    // Quadrature Points
    const vector<Real > & quadWeight = quadrature->getQuadWeights();
    assert(shape.size() == quadWeight.size());

    // Problem size and dimension
    const uint numQP = shape.size(); //Shape function object for each Quadrature point    
    const uint nodePerElem = _nodesID.size(); //Includes 2ring neighbors
    assert(nodePerElem == nodesX.size());
    uint dim = 2;

    // Clear/resize QPweights, N, and DN
    _QPweights.resize(numQP, 0.0);
    _N.resize(nodePerElem * numQP, 0.0);
    _DN.resize(nodePerElem * numQP, VectorXd::Zero(dim));
    _DDN.resize(nodePerElem * numQP, Matrix2d::Zero());

    // Loop over quad points and fill N and DN
    for(uint q = 0; q < numQP; q++)
      {      
	// Store QPweight
	_QPweights[q] = quadWeight[q];
      
	// Computing derivatives and shape function values
	for(uint a = 0; a < nodePerElem; a++) {
	  _N[q*nodePerElem + a] = shape[q]->getN(a); // it can be in a different loop, here just for convenience
	  for(uint i = 0; i < dim; i++) {
	    _DN[q*nodePerElem + a](i) += shape[q]->getDN(a,i);
	    for(uint j=0; j < dim; j++) {
	      //shape[q] is a pointer to Shape class which does not have getDDN. The pointer
	      // has to be downcast to LoopShellShape* type to be able to access the method
	      _DDN[q*nodePerElem + a](i,j) += dynamic_cast<LoopShellShape*>(shape[q])->getDDN(a,i,j);
	    }//j loop
	  } // i loop
	} // a loop

      } // loop over quad points

  } // end LoopShellElement constructor

} // namespace voom
