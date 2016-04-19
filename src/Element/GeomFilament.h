//-*-C++-*-
/*!
  \file GeomFilament.h

  \brief Base class implementation for Geometry Filament.
  The derived classes will compute shape functions and shape functions 
  derivatives with respect to the global coordinate system.
  There is no physics implemented in this and the derived classes. 
*/

#if !defined(__GeomFilament_h__)
#define __GeomFilament_h__

namespace voom {

  class GeomFilament {    
    
  public:
    /*!
      Constructor for parent geometry element.
      \param elemID is element ID
      \param nodesID is a vector of Node ID
    */
    GeomFilament(const int elemID, const vector<int > & nodesID, const vector<Vector3d> & nodesX):
      _elemID(elemID), _nodesID(nodesID),_nodesX(nodesX) {}
    
    //! Get element ID
    int getGeomFilamentID() {return _elemID; }

    //! Get nodes per element
    uint getNodesPerFilament() {return _nodesID.size(); }

    //! Get node list
    const vector<int > & getNodesID() {return _nodesID; }

    //! Get node ref position
    const vector<Vector3d> getNodesX(){return _nodesX;}

  protected:
    //! GeomFilament ID
    int                    _elemID;

    //! List of nodes
    vector<int >           _nodesID;

    //! Nodes ref position
    vector<Vector3d>       _nodesX;
    
  }; // GeomFilament

}; // namespace voom

#endif
