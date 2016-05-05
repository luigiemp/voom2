//-*-C++-*-
/*!
  \file GelElement.h

  \brief Base class implementation for Geometry Element.
  The derived classes will compute shape functions and shape functions 
  derivatives with respect to the global coordinate system.
  There is no physics implemented in this and the derived classes. 
*/

#if !defined(__GelElement_h__)
#define __GelElement_h__

namespace voom {

  class GelElement {    
    
  public:
    /*!
      Constructor for parent geometry element.
      \param elemID is element ID
      \param nodesID is a vector of Node ID
    */
    GelElement(const int elemID, const vector<int > & nodesID, const vector<Vector3d> & nodesX):
      _elemID(elemID), _nodesID(nodesID),_nodesX(nodesX) {}

    ~GelElement(){};
    
    //! Get element ID
    int getGelElementID() {return _elemID; }

    //! Get nodes per element
    uint getNodesPerElement() {return _nodesID.size(); }

    //! Get node list
    const vector<int > & getNodesID() {return _nodesID; }

    //! Get node ref position
    const vector<Vector3d> getNodesX(){return _nodesX;}

  protected:
    //! GelElement ID
    int                    _elemID;

    //! List of nodes
    vector<int >           _nodesID;

    //! Nodes ref position
    vector<Vector3d>       _nodesX;
    
  }; // GelElement

}; // namespace voom

#endif
