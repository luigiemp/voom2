//-*-C++-*-
#ifndef __CrossLinkManager_h__
#define __CrossLinkManager_h__

#include "GelModel.h"
#include "FilamentMaterial.h"
#include "EigenResult.h"
#include "GelMesh.h"
#include "CrossLink.h"
#include "GelElement.h"
#include "PeriodicBox.h"

namespace voom{

  // Crosslink manager: initialize and update crosslinks based
  // on some update rule object (KMC for example)
  class CrossLinkManager {

  public:

    //! Basic Constructor
    CrossLinkManager(GelModel* aGelModel, FilamentMaterial *  clMat, Real maxClL, PeriodicBox* box);
    		     
    //! Destructor
    
    ~CrossLinkManager(){}

    void computeRates(Real k0 = 1.0);

    void KMC(int kmcIt, Real maxDist);

    int findEvent();
    
    void addCrosslink(pair<int, int> nodePair);

    void getCurrentConfig(){_myGelModel->getField(_x);};

    void carryOutEvent(int i);

    pair<int,Real> findClosest(Vector3d  node , vector< double >  nodeList);
    
    void findNeighbours(Real maxDist);

    void findKMCPairsCandidate();

  protected:
    
    FilamentMaterial *  _clMat;
    GelModel * _myGelModel;
    Real _maxClL;
    VectorXd _X;
    GelMesh* _myGelMesh;
    
    int _dim;
    int _NumNode;
    vector< double > _x;
    vector< pair<bool,CrossLink * > > _crossLinkMap;
    vector<Real> _rates; 
    vector<Real> _cumulativeRates; 
    PeriodicBox* _box;
    vector<pair<int,int> > _KMCPairsCandidates;
    vector<vector<int> > _neighboursMap;
  };

} // namespace voom

#endif
