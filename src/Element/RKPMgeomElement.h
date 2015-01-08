//-*-C++-*-
/*!
  \file RKPMgeomElement.h

  \brief Base class implementation for RKPM geometry element class.
  This class stores the RKPM shape functions evaluated at selected QP or
  nodal positions at construction time.
  There is no physics implemented in this class. 
*/

#if !defined(__RKPMgeomElement_h__)
#define __RKPMgeomElement_h__

#include "GeomElement.h"
#include "MRKPMShape.h"

namespace voom {

  class RKPMgeomElement: public GeomElement {    
  public:
    RKPMgeomElement(int elemID, const vector<int > & nodesID, 
		   const vector<VectorXd > & NodesX, 
		   const vector<VectorXd > & MaterialPoints, 
		   const vector<Real > & Weights,
		   Real support, Real radius, Real supportHat);

    // RKPMgeomElement distructor to delete RKPMshape functions
    ~RKPMgeomElement() { 
      for(uint i = 0; i < _RKPMshapes.size(); i++)
	delete _RKPMshapes[i];
    }
    
    //! Get number of quadrature points
    uint getNumberOfQuadPoints() {return _QPweights.size(); }

    //! Get weight for quadrature point q
    Real getQPweights(uint q) {return _QPweights[q]; }

    //! Get shape functions values at quadrature point q, node a
    Real getN(uint q, uint a) {
      return _RKPMshapes[q]->getN(a); 
    }

    //! Get shape functions derivatives at quadrature point q, node a, direction i
    Real getDN(uint q, uint a, uint i) {
      return _RKPMshapes[q]->getDN(a,i); 
    }

  protected:
    //! Mesh Free nodes cluster considered in this element
    vector<VectorXd >      _nodesX;
    
    //! Quadrature weights
    vector<Real >          _QPweights;

    /*! 
      Shape functions objects evaluated at quadrature points
    */
    vector<MRKPMShape* >   _RKPMshapes;

  };

} // namespace voom

#endif
