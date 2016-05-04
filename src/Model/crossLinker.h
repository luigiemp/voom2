//-*-C++-*-
#ifndef __crossLinker_h__
#define __crossLinker_h__

#include "GelModel.h"
#include "FilamentMaterial.h"
#include "EigenResult.h"
#include "GelMesh.h"
#include "CrossLink.h"
#include "GelElement.h"

namespace voom{

  // Crosslinker class: find and update crosslinks with Kintetic Monte Carlo
  class crossLinker {

  public:

    //! Basic Constructor
    crossLinker(GelModel* aGelModel, FilamentMaterial *  clMat, Real maxClL);
    		     
    //! Destructor
    
    ~crossLinker() {
      
      set<CrossLink *> UNIQUEcrosslinks;
      for (uint i = 0; i < _crosslinks.size(); i++) 
	UNIQUEcrosslinks.insert(_crosslinks[i]);

      for (set<CrossLink *>::iterator it = UNIQUEcrosslinks.begin();
      	   it != UNIQUEcrosslinks.end(); it++) 
	delete (*it);
	
    };

    bool KMC();
    
    void addCrosslink(pair<int, Real> arg_dist, int i);

    void updateCrosslinks();

    void updateSingleCL(int i);

    void getCurrentConfig();

    pair<int,Real> findClosest(Vector3d  node , vector< double >  nodeList);
    

  protected:
    
    //! List of Material data at each element in the model
    FilamentMaterial *  _clMat;
    GelModel * _myGelModel;
    Real _maxClL;
    vector<CrossLink * > _crosslinks;
    VectorXd _X0;
    GelMesh* _myGelMesh;
    vector<bool> _CrossLinkedNodes;
    int _dim;
    int _NumNode;
    vector< double > _x;
   
  };

} // namespace voom

#endif
