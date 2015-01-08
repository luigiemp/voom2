//-*-C++-*-

#include "LMEgeomElement.h"

namespace voom{

  LMEgeomElement::LMEgeomElement(int elemID, const vector<int > & nodesID, 
		                 const vector<VectorXd > & NodesX, 
		                 const vector<VectorXd > & MaterialPoints, 
				 const vector<Real > & Weights,
		                 Real beta, Real tol, uint MaxIter):
    GeomElement(elemID, nodesID), _nodesX(NodesX), _QPweights(Weights)
  {
    uint NumMP = Weights.size();
    assert(NumMP == MaterialPoints.size());
    _LMEshapes.reserve(NumMP);
    // Compute LME shape functions at given MaterialPoints
    for(uint q=0; q<NumMP; q++)
    {
      LMEShape* LMEshape = new LMEShape(_nodesX, MaterialPoints[q], 
					beta, tol, MaxIter);
      _LMEshapes.push_back(LMEshape);
    } // loop over material points

  } // end LMEgeomElement constructor

} // namespace voom
