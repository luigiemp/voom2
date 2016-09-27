//-*-C++-*-
/*!
  \file GeomElement.h

  \brief Base class implementation for Geometry element.
  The derived classes will compute shape functions and shape functions 
  derivatives with respect to the global coordinate system.
  There is no physics implemented in this and the derived classes. 
*/

#if !defined(__GeomElement_h__)
#define __GeomElement_h__

#include "Quadrature.h"
#include "Shape.h"

namespace voom {

  class GeomElement {    

  public:
    /*!
      Constructor for parent geometry element.
      \param elemID is element ID
      \param nodesID is a vector of Node ID
    */
    GeomElement(const int elemID, const vector<int > & nodesID):
      _elemID(elemID), _nodesID(nodesID) {}
    
    //! Get element ID
    int getGeomElementID() {return _elemID; }

    //! Get nodes per element
    uint getNodesPerElement() {return _nodesID.size(); }

    //! Get node list
    const vector<int > & getNodesID() {return _nodesID; }

    // Virtual functions //
    //! Get number of quadrature points
    virtual uint getNumberOfQuadPoints() = 0;
    
    //! Get weight for quadrature point q
    virtual Real getQPweights(uint q) = 0;

    //! Get shape functions values at quadrature point q, node a
    virtual Real getN(uint q, uint a) = 0; 

    //! Get shape functions derivatives at quadrature point q, node a, direction i
    virtual Real getDN(uint q, uint a, uint i) = 0;

  protected:
    //! GeomElement ID
    int                    _elemID;

    //! List of nodes
    vector<int >           _nodesID;
    
  }; // GeomElement

}; // namespace voom

#endif
