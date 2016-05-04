#include "crossLinker.h"
#include "CrossLink.h"
#include "FilamentMaterial.h"
#include <time.h>
namespace voom {

  // Constructor
  crossLinker::crossLinker(GelModel* aGelModel, FilamentMaterial *  clMat, Real maxClL):
    _myGelModel(aGelModel) , _clMat(clMat), _maxClL(maxClL)
  {
    /* initialize random seed: */
    srand (time(NULL));
    // Ref config:
    VectorXd _X0 = _myGelModel->getX();

    // Current config:
    int _NumNode= _myGelModel->getNumberOfNodes();
    int _dim = _myGelModel->getDimension();
    
    _x.resize(_dim*_NumNode);
    _CrossLinkedNodes.resize(_NumNode,false);

    _myGelModel->getField(_x);
    
    Real maxClL2 = pow(_maxClL,2);
    
    cout << "Starting crosslinking" << endl;
    // Loop over nodes and create crosslink if mindist< distance < maxdist, using current config
    // Each node can be crosslinked only once
    for(int i = 0; i< _NumNode ; i++){
      if(_CrossLinkedNodes[i] == false){

	Vector3d node1;
	node1 << _x[i*_dim+0],_x[i*_dim+1],_x[i*_dim+2]; 

	// Find closest node distance and ID:
	pair<int,Real> arg_dist = this->findClosest(node1, _x);

	// If it meets criteria, add to list of crosslinks
	if(_CrossLinkedNodes[arg_dist.first]== false & arg_dist.second < maxClL2 ){
	  
	  this->addCrosslink(arg_dist,i);
	}
      }

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
    X1 << _X0(nodesID[0]*_dim) , _X0(nodesID[0]*_dim+1), _X0(nodesID[0]*_dim+2);
    Vector3d X2;
    X2 << _X0(nodesID[1]*_dim) , _X0(nodesID[1]*_dim+1), _X0(nodesID[1]*_dim+2);

    NodeX.push_back(X1);
    NodeX.push_back(X2);
	
    // Construct cl element:
    GelElement* CLelem = new GelElement(0,nodesID, NodeX);
	
    // Construct cl object: cl element +  material prop:
    CrossLink* CL = new CrossLink(CLelem, _clMat);
	
    _myGelModel->addCrosslink(CL);
    _CrossLinkedNodes[nodesID[0]] = true;
    _CrossLinkedNodes[nodesID[1]] = true;
  }

  bool crossLinker::KMC(){
    Real p = ((double) rand() / (RAND_MAX));
    if(p<0.1){
      return true;
    }else{
      return false;
    }

    // // Compute enrergy of candidate crosslink:
    // FilamentMaterial::Filresults Rf;
    // _clMat->compute(Rf,xlist);
    // Real energy = Rf.W;
  }
  void crossLinker::getCurrentConfig(){
    _myGelModel->getField(_x);
  }

  void crossLinker::updateSingleCL(int i){
    
    this->getCurrentConfig();

    if(_CrossLinkedNodes[i] == false){
      vector<int> nodesID;
      vector<Vector3d> NodeX;

      Vector3d node1;
      node1 << _x[i*_dim+0],_x[i*_dim+1],_x[i*_dim+2]; 

      // Find closest node distance and ID:
      pair<int,Real> arg_dist = this->findClosest(node1, _x);
	
      Vector3d node2;
      node2 << _x[arg_dist.first*_dim+0],_x[arg_dist.first*_dim+1],_x[arg_dist.first*_dim+2]; 
	
      vector<Vector3d> xlist;
      xlist.push_back(node1);
      xlist.push_back(node2);
	
      bool link = this->KMC();
      if(link && _CrossLinkedNodes[i] == false){
	this->addCrosslink(arg_dist,i);
      }
    }else{
      

    }
    
  }

  // update CL funciton: takes as input the update rule
  void crossLinker::updateCrosslinks(){
    // First go over existing Crosslinks, decide if they stay
    int numCL = _crosslinks.size();
    int i = 0;
    while(i < _crosslinks.size()){
      bool link = this->KMC();
      if(link){
	const vector<int> NodesID = _crosslinks[i]->getNodesID();
	_CrossLinkedNodes[NodesID[0]] = false;
	_CrossLinkedNodes[NodesID[1]] = false;
	_myGelModel->deleteCrosslink(i);
	
      }else {
	++i;
      }
    }
    // Then go over unlinked nodes, decide if they should form crosslinks:
    
    // Current config:
    
    this->getCurrentConfig();
    
    for(int i = 0; i< _NumNode ; i++){
      if(_CrossLinkedNodes[i] == false){
	vector<int> nodesID;
	vector<Vector3d> NodeX;

	Vector3d node1;
	node1 << _x[i*_dim+0],_x[i*_dim+1],_x[i*_dim+2]; 

	// Find closest node distance and ID:
	pair<int,Real> arg_dist = this->findClosest(node1, _x);
	
	Vector3d node2;
	node2 << _x[arg_dist.first*_dim+0],_x[arg_dist.first*_dim+1],_x[arg_dist.first*_dim+2]; 
	
	vector<Vector3d> xlist;
	xlist.push_back(node1);
	xlist.push_back(node2);
	
	bool link = this->KMC();
	if(link && _CrossLinkedNodes[i] == false){
	  this->addCrosslink(arg_dist,i);
	}
      }
    }
  };


  // find closest node:
  pair<int, Real> crossLinker::findClosest(Vector3d  node , vector< double >  nodeList){
    
    Real minDist = 1.0;
    int minArg = 0;
    Real epsilon = 1e-6;
      
    for(int i = 0; i<nodeList.size()/_dim; i++){
      Real d2 = 0.0;
      for(int j = 0; j <_dim; j++){
	d2 += pow((node(j)-nodeList[i*_dim+j]),2);
      }
            
      if(d2 < minDist && d2> epsilon){
	minDist = d2;
      	minArg = i;
      }
    }
    return make_pair(minArg, minDist);
  }

} // namespace voom
