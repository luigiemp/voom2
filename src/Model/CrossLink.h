//-*-C++-*-%
#ifndef __CrossLink_h__
#define __CrossLink_h__

#include "Model.h"
#include "FilamentMaterial.h"
#include "EigenResult.h"
#include "GelMesh.h"
//#include "GelModel.h"

namespace voom{

  // Model Results
  class CrossLink {
 
  public:

    //! Basic Constructor
    /*! Construct from basic data structures defining the mesh, materials, BCs. 
     */
    CrossLink(GelElement* aCrossLink, FilamentMaterial* spring);
		     
    //! Destructor
    
    ~CrossLink() {};
    
    // Get number of Nodes:
    int getNumberOfNodes(){return _myCrossLink->getNodesID().size();}

    const vector<int> & getNodesID(){return _myCrossLink->getNodesID();}

    //! Solve the system
    void compute(Result & R, vector<Vector3d> & d);


  protected:
    
    //void computeDeformation(vector<Vector3d > & dlist, GeomElement* geomEl);
    //! List of Material data at each element in the model
    FilamentMaterial*  _spring;
    GelElement* _myCrossLink;
    int _nodeNum;
    int _dim;
    vector<int> _NodesID;
    vector<Vector3d> _X;
  };

} // namespace voom

#endif
