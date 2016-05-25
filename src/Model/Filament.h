//-*-C++-*-%
#ifndef __Filament_h__
#define __Filament_h__

#include "Model.h"
#include "FilamentMaterial.h"
#include "EigenResult.h"
#include "GelMesh.h"
//#include "GelModel.h"
#include "PeriodicBox.h"

namespace voom{

  // Model Results
  class Filament {
 
  public:

    //! Basic Constructor
    /*! Construct from basic data structures defining the mesh, materials, BCs. 
     */
    Filament(GelElement* aFilament, FilamentMaterial* spring,FilamentMaterial* angleSpring, PeriodicBox* box);
		     
    //! Destructor
    
    ~Filament() {};
    
    // Get number of Nodes:
    int getNumberOfNodes(){return _myFilament->getNodesID().size();}

    const vector<int> & getNodesID(){return _myFilament->getNodesID();}

    //! Solve the system
    void compute(Result & R, vector<Vector3d> & d);


  protected:
    
    //void computeDeformation(vector<Vector3d > & dlist, GeomElement* geomEl);
    //! List of Material data at each element in the model
    FilamentMaterial*  _spring;
    FilamentMaterial*  _angleSpring;
    GelElement* _myFilament;
    int _nodeNum;
    int _dim;
    vector<int> _NodesID;
    vector<Vector3d> _X;
    PeriodicBox* _box;
  };

} // namespace voom

#endif
