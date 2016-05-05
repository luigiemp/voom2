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

  // Crosslink manager: initialize and update crosslinks based
  // on some update rule object (KMC for example)
  class crossLinker {

  public:

    //! Basic Constructor
    crossLinker(GelModel* aGelModel, FilamentMaterial *  clMat, Real maxClL);
    		     
    //! Destructor
    
    ~crossLinker(){}

    void computeRates(Real k0 = 1.0);

    void KMC();

    int findEvent();

    void addCrosslink(pair<int, Real> arg_dist, int i);

    void updateSingleCL(int i);

    void getCurrentConfig(){_myGelModel->getField(_x);};

    void ClUpdate(int i);

    pair<int,Real> findClosest(Vector3d  node , vector< double >  nodeList);
    

  protected:
    
    FilamentMaterial *  _clMat;
    GelModel * _myGelModel;
    Real _maxClL;
    VectorXd _X;
    GelMesh* _myGelMesh;
    vector<bool> _CrossLinkedNodes;
    int _dim;
    int _NumNode;
    vector< double > _x;
    vector< CrossLink *> _crossLinkMap;
    vector<Real> _rates; 
    vector<Real> _cumulativeRates; 
   
  };

} // namespace voom

#endif
