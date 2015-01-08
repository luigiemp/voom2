//-*-C++-*-
/*!
  \file LMEgeomElement.h

  \brief Base class implementation for Max Ent geometry element class.
  This class stores the LME shape functions evaluated at selected QP 
  at construction time.
  There is no physics implemented in this class. 
*/

#if !defined(__LMEgeomElement_h__)
#define __LMEgeomElement_h__

#include "GeomElement.h"
#include "LMEShape.h"

namespace voom {

  class LMEgeomElement: public GeomElement {    
  public:
    LMEgeomElement(int elemID, const vector<int > & nodesID, 
		   const vector<VectorXd > & NodesX, 
		   const vector<VectorXd > & MaterialPoints, 
		   const vector<Real > & Weights,
		   Real beta, Real tol, uint MaxIter);

    // LMEgeomElement distructor to delete LMEshape functions
    ~LMEgeomElement() { 
      for(uint i = 0; i < _LMEshapes.size(); i++)
	delete _LMEshapes[i];
    }
    
    //! Get number of quadrature points
    uint getNumberOfQuadPoints() {return _QPweights.size(); }

    //! Get weight for quadrature point q
    Real getQPweights(uint q) {return _QPweights[q]; }

    //! Get shape functions values at quadrature point q, node a
    Real getN(uint q, uint a) {
      return _LMEshapes[q]->getN(a); 
    }

    //! Get shape functions derivatives at quadrature point q, node a, direction i
    Real getDN(uint q, uint a, uint i) {
      return _LMEshapes[q]->getDN(a,i); 
    }

  protected:
    //! Mesh Free nodes cluster considered in this element
    vector<VectorXd >      _nodesX;
    
    //! Quadrature weights
    vector<Real >          _QPweights;

    /*! 
      Shape functions objects evaluated at quadrature points
    */
    vector<LMEShape* >     _LMEshapes;

  };

} // namespace voom

#endif
