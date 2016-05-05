#include "crossLinker.h"
#include "CrossLink.h"
#include "FilamentMaterial.h"
#include <time.h>
namespace voom {

  // Constructor: generates initial crosslinks, based on a length constaint
  crossLinker::crossLinker(GelModel* aGelModel, FilamentMaterial *  clMat, Real maxClL):
    _myGelModel(aGelModel) , _clMat(clMat), _maxClL(maxClL)
  {
    // initialize random seed: 
    srand (time(NULL));
    // Ref config:
    _X = _myGelModel->getX();
    // Current config:
    _NumNode= _myGelModel->getNumberOfNodes();
    _dim = _myGelModel->getDimension(); 
    _x.resize(_dim*_NumNode);
    _rates.resize(_NumNode,0.0);
    _cumulativeRates.resize(_NumNode,0.0);
    this->getCurrentConfig();

    // Initialize crosslink map (which for each node, stores pointer to corresponding CL if exists):
    _crossLinkMap.resize(_NumNode, NULL);
    
    Real maxClL2 = pow(_maxClL,2);
    
    cout << "Starting crosslinking" << endl;
    // Initial crosslinking: Loop over nodes and create crosslink if 
    // mindist< distance < maxdist, using current config.
    // Each node can be crosslinked only once.
    for(int i = 0; i< _NumNode ; i++){
      if(_crossLinkMap[i] == NULL){
	Vector3d node1;
	node1 << _x[i*_dim+0],_x[i*_dim+1],_x[i*_dim+2]; 

	// Find closest node distance and ID:
	pair<int,Real> arg_dist = this->findClosest(node1, _x);
	
	// If it meets criteria, add to list of crosslinks
	if(_crossLinkMap[arg_dist.first]== NULL & arg_dist.second < maxClL2 ){
	  this->addCrosslink(arg_dist,i);
	}
      }

    }
  }

  void crossLinker::computeRates(Real k0){
    this->getCurrentConfig();
    
    for(int i = 0; i<_NumNode; i++){
      // For each node, check if a crosslink is already there
      // If yes, compute the energy of that crosslink
      // Otherwise, find the closest node and compute the energy of the would be crosslink:
      Vector3d x1;
      Vector3d x2;
      if(_crossLinkMap[i] != NULL){
	const vector<int> nodesID = _crossLinkMap[i]->getNodesID();
	// Get the current position for the 2 nodes
	x1 << _x[nodesID[0]*_dim] , _x[nodesID[0]*_dim+1], _x[nodesID[0]*_dim+2];
	x2 << _x[nodesID[1]*_dim] , _x[nodesID[1]*_dim+1], _x[nodesID[1]*_dim+2];
	
	vector<Vector3d> xlist;
	xlist.push_back(x1);
	xlist.push_back(x2);
      
	// Compute the energy:
	FilamentMaterial::Filresults Rf;
	_clMat->compute(Rf, xlist);
	_rates[i] = k0*exp(Rf.W);

      } else{
	
	x1 << _x[i*_dim+0],_x[i*_dim+1],_x[i*_dim+2]; 
	// Find closest node:
	pair<int,Real> arg_dist = this->findClosest(x1, _x);
	x2 << _x[arg_dist.first*_dim] , _x[arg_dist.first*_dim+1], _x[arg_dist.first*_dim+2];
	
	vector<Vector3d> xlist;
	xlist.push_back(x1);
	xlist.push_back(x2);
      
	// Compute the energy:
	FilamentMaterial::Filresults Rf;
	_clMat->compute(Rf, xlist);
	_rates[i] = k0*exp(-Rf.W);
	
      }
      
      
      for(int k = 0 ; k < i+1 ; k++){
	_cumulativeRates[i] += _rates[k];
      }
     }
  }

  int crossLinker::findEvent(){
    Real u = ((double) rand()/RAND_MAX);
    Real Q = _cumulativeRates[_NumNode];

    for(int i = 0; i < _NumNode; i++){
      if(_cumulativeRates[i] > u*Q){
	return i;
      }
    }
  }
  
  void crossLinker::KMC(){
    Real k0 = 1.0;
    this->computeRates(k0);
    int i = this->findEvent();
    
    this->ClUpdate(i);
  }
  
  void crossLinker::ClUpdate(int i){
    
    // Detach if linked, attach otherwise:
    if(_crossLinkMap[i] != NULL){
      const vector<int> NodesID = _crossLinkMap[i]->getNodesID();
      _myGelModel->deleteCrosslink(_crossLinkMap[i]);
      _crossLinkMap[NodesID[0]] = NULL;
      _crossLinkMap[NodesID[1]] = NULL;

    }else {

      this->getCurrentConfig();
      Vector3d node1;
      node1 << _x[i*_dim+0],_x[i*_dim+1],_x[i*_dim+2]; 
      
      // Find closest node and distance:
      pair<int,Real> arg_dist = this->findClosest(node1, _x);
      Real maxClL2 = pow(_maxClL,2);
      // Attach
      this->addCrosslink(arg_dist,i);
    }

  }

  void crossLinker::addCrosslink( pair<int, Real> arg_dist, int i){
    
    vector<int> nodesID;
    vector<Vector3d> NodeX;
    // Get the nodes ID for the cross link
    nodesID.push_back(i);
    nodesID.push_back(arg_dist.first);

    // Get the ref position for the 2 nodes
    Vector3d X1;
    X1 << _X(nodesID[0]*_dim) , _X(nodesID[0]*_dim+1), _X(nodesID[0]*_dim+2);
    Vector3d X2;
    X2 << _X(nodesID[1]*_dim) , _X(nodesID[1]*_dim+1), _X(nodesID[1]*_dim+2);

    NodeX.push_back(X1);
    NodeX.push_back(X2);
	
    // Construct cl element:
    GelElement* CLelem = new GelElement(0,nodesID, NodeX);
	
    // Construct cl object: gel element + cl  material prop:
    CrossLink* CL = new CrossLink(CLelem, _clMat);

    _myGelModel->addCrosslink(CL);
    _crossLinkMap[nodesID[0]] = CL;
    _crossLinkMap[nodesID[1]] = CL;
  }

 
  // find closest not attached node (use parsimoniously):
  pair<int, Real> crossLinker::findClosest(Vector3d  node , vector< double >  nodeList){
    
    Real minDist = 1.0;
    int minArg = 0;
    Real epsilon = 1e-6;
   
    for(int i = 0; i<_NumNode; i++){
      Real d2 = 0.0;
      for(int j = 0; j <_dim; j++){
	d2 += pow((node(j)-nodeList[i*_dim+j]),2);
      }
      
      if(d2 < minDist && d2> epsilon && _crossLinkMap[i] == NULL){
	minDist = d2;
      	minArg = i;
      }
    }
    return make_pair(minArg, minDist);
  }

} // namespace voom
