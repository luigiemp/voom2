//-*-C++-*-

#include "RKPMgeomElement.h"

namespace voom{

  RKPMgeomElement::RKPMgeomElement(int elemID, const vector<int > & nodesID, 
				   const vector<VectorXd > & NodesX, 
				   const vector<VectorXd > & MaterialPoints, 
				   const vector<Real > & Weights,
				   Real support, Real radius,
				   Real supportHat):
    GeomElement(elemID, nodesID), _nodesX(NodesX), _QPweights(Weights)
  {
    uint NumMP = Weights.size();
    assert(NumMP == MaterialPoints.size());
    _RKPMshapes.reserve(NumMP);
    // Compute RKPM shape functions at given MaterialPoints
    for(uint q=0; q < NumMP; q++)
    {
      MRKPMShape* RKPMshape = new MRKPMShape(_nodesX, MaterialPoints[q], 
					     support, radius, supportHat);
      _RKPMshapes.push_back(RKPMshape);
    } // loop over material points

  } // end RKPMgeomElement constructor

} // namespace voom
