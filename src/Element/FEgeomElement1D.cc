#include "FEgeomElement1D.h"

namespace voom {

  // Constructor
  FEgeomElement1D::FEgeomElement1D(const int elemID, const vector<int > & nodesID, 
				   const vector<VectorXd > & nodesX, 
				   vector<Shape* > shape, Quadrature* quadrature,
				   const Real radius):
    GeomElement(elemID, nodesID)
  {
    
    Vector3d ang = nodesX[1];
    ang -= nodesX[0];
    Real mag = ang.norm();
    Real fac = 2.0/mag;

    // Quadrature Points
    const vector<Real > & quadWeight = quadrature->getQuadWeights();
    assert(shape.size() == quadWeight.size());
    // FEgeomElement1D always has 2 nodes
    assert(nodesID.size() == 2);
    assert(nodesX.size() == 2);
    const uint numQP = shape.size();    

    // Clear/resize QPweights, N, and DN
    _QPweights.resize(numQP, 0.0);
    _N.resize(2 * numQP, 0.0);
    _DN.resize(2 * numQP, Vector3d::Zero());
    
    // Loop over quad points
    for(uint q = 0; q < numQP; q++)
    {
      // Compute shape functions
      _N[q*2]     = shape[q]->getN(0);
      _N[q*2 + 1] = shape[q]->getN(1);

      // Compute shape function derivatives
      _DN[q*2](0)     = shape[q]->getDN(0,0)*ang(0)*fac;
      _DN[q*2](1)     = shape[q]->getDN(0,0)*ang(1)*fac; 
      _DN[q*2](2)     = shape[q]->getDN(0,0)*ang(2)*fac; 
      _DN[q*2 + 1](0) = shape[q]->getDN(1,0)*ang(0)*fac;
      _DN[q*2 + 1](1) = shape[q]->getDN(1,0)*ang(1)*fac; 
      _DN[q*2 + 1](2) = shape[q]->getDN(1,0)*ang(2)*fac; 
      
      // Compute quadrature weights
      _QPweights[q] =  M_PI * square(radius) * mag * 2.0 * quadWeight[q];
     
    } // q loop

  } // FEgeomElement1D

} //  namespace voom
