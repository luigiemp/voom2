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
    GeomFilament(const int elemID, const vector<int > & nodesID):
      _elemID(elemID), _nodesID(nodesID) {}
    
    //! Get element ID
    int getGeomFilamentID() {return _elemID; }

    //! Get nodes per element
    uint getNodesPerFilament() {return _nodesID.size(); }

    //! Get node list
    const vector<int > & getNodesID() {return _nodesID; }


  protected:
    //! GeomFilament ID
    int                    _elemID;

    //! List of nodes
    vector<int >           _nodesID;

    
  }; // GeomFilament

}; // namespace voom

#endif
