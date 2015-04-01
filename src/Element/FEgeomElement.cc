//-*-C++-*-

#include "FEgeomElement.h"

namespace voom{

  FEgeomElement::FEgeomElement(const int elemID, const vector<int > & nodesID, 
			       const vector<VectorXd > & nodesX,
			       vector<Shape* > shape, Quadrature* quadrature): 
			       GeomElement(elemID, nodesID)
  {
    // Compute shape functions, shape functions derivatives and QP weights
    
    // Quadrature Points
    const vector<Real > & quadWeight = quadrature->getQuadWeights();
    assert(shape.size() == quadWeight.size());

    // Problem size and dimension
    const uint numQP = shape.size();    
    const uint nodePerElem = _nodesID.size();
    assert(nodePerElem == nodesX.size());
    uint dim = nodesX[0].size();

    // Clear/resize QPweights, N, and DN
    _QPweights.resize(numQP, 0.0);
    _N.resize(nodePerElem * numQP, 0.0);
    _DN.resize(nodePerElem * numQP, VectorXd::Zero(dim));

    // To be done differently !! This is just a quick fix to be used with the pressure elements
    if ( ((quadrature->getQuadPoints())[0]).size() == dim) { 
      // Loop over quad points and fill N and DN
      for(uint q = 0; q < numQP; q++)
	{      
	  // Compute trasformation Jacobian
	  MatrixXd J = MatrixXd::Zero(dim, dim), Jinv = MatrixXd::Zero(dim ,dim);
	  for(uint i = 0; i < dim; i++)
	    for(uint j = 0; j < dim; j++)
	      for(uint a = 0; a < nodePerElem; a++)
		J(i,j) += shape[q]->getDN(a,j)*nodesX[a](i);
      
	  Real detJ = fabs(J.determinant());
	  Jinv = J.inverse();
	  // Store QPweight
	  _QPweights[q] = quadWeight[q]*detJ;
	  // cout << "det J " << detJ << endl
      
	  // Computing spatial derivatives and shape function values
	  for(uint a = 0; a < nodePerElem; a++) {
	    _N[q*nodePerElem + a] = shape[q]->getN(a); // it can be in a different loop, here just for convenience
	    for(uint i = 0; i < dim; i++) {
	      for(uint j = 0; j < dim; j++) { 
		_DN[q*nodePerElem + a](i) += shape[q]->getDN(a,j) *Jinv(j,i);
	      } // j loop
	    } // i loop
	  } // a loop

	} // loop over quad points
    } // WARNING
    else {
      // Computing shape function values
      for(uint q = 0; q < numQP; q++)
      {  
	for(uint a = 0; a < nodePerElem; a++) {
	  _N[q*nodePerElem + a] = shape[q]->getN(a);
	  _DN[q*nodePerElem + a](0) = shape[q]->getDN(a,0);
	  _DN[q*nodePerElem + a](1) = shape[q]->getDN(a,1);
	}
	_QPweights[q] = quadWeight[q];
      }
    }

  } // end FEgeomElement constructor

} // namespace voom
