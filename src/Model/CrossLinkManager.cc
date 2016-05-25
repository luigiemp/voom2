#include "CrossLinkManager.h"
#include "CrossLink.h"
#include "FilamentMaterial.h"
#include <time.h>

namespace voom {

  // Constructor: generates initial crosslinks, based on a length constaint maxClL
  CrossLinkManager::CrossLinkManager(GelModel* aGelModel, FilamentMaterial *  clMat, Real maxClL, PeriodicBox* box):
    _myGelModel(aGelModel) , _clMat(clMat), _maxClL(maxClL),_box(box)
  {
    // initialize random seed: 
    srand (time(NULL));
    // Ref config:
    _X = _myGelModel->getX();
    // Current config:
    _NumNode= _myGelModel->getNumberOfNodes();
    _dim = _myGelModel->getDimension(); 
    _x.clear();
    _x.resize(_dim*_NumNode);
        
    this->getCurrentConfig();

    // Initialize crosslink map (which, for each node, stores pointer to corresponding CL if exists):
    _crossLinkMap.resize(_NumNode);
    
    Real maxClL2 = pow(_maxClL,2);
    
    cout << "Starting crosslinking" << endl;
    // Initial crosslinking: Loop over nodes and create crosslink if 
    // mindist< distance < maxdist, using current config.
    // Each node can be crosslinked only once.
    for(int i = 0; i< _NumNode ; i++){
      if(_crossLinkMap[i].first == false){
    	Vector3d node1;
    	node1 << _x[i*_dim+0],_x[i*_dim+1],_x[i*_dim+2]; 

    	// Find closest node distance and ID:
    	pair<int,Real> arg_dist = this->findClosest(node1, _x);
	
    	// If its not already crosslinked and it meets the length criteria, add to list of crosslinks
    	if(_crossLinkMap[arg_dist.first].first == false & arg_dist.second < maxClL2 ){
    	  this->addCrosslink(make_pair(arg_dist.first,i));
    	}
      }

    }
  }

  // find closest (and not attached node) (use parsimoniously):
  pair<int, Real> CrossLinkManager::findClosest(Vector3d  node , vector< double >  nodeList){
    
    Real minDist = 1.0;
    int minArg = 0;
    Real epsilon = 1e-6;
   
    for(int i = 0; i<_NumNode; i++){
      Real d2 = 0.0;
      Vector3d dx;
      for(int j = 0; j <_dim; j++){
	dx(j) = (node(j)-nodeList[i*_dim+j]);
      }
      // Periodic boudnary condition:
      _box->mapDistance(dx);
      d2 = dx.squaredNorm();
      
      if(d2 < minDist && d2> epsilon && _crossLinkMap[i].first == false){
	minDist = d2;
      	minArg = i;
      }
    }
    return make_pair(minArg, minDist);
  }

  void CrossLinkManager::findNeighbours(Real maxDist){
    /* Build a list of node neihbours, for each node.*/
    
    _neighboursMap.clear();

    Real epsilon = 1e-6;
    Real maxDist2 = pow(maxDist,2);

    //first lets build a vector of node index, on which we will loop and 
    //remove as we go to make sure we don't add nodes twice
    vector<int> nodelist;
    for(int i = 0; i< _NumNode; i++){
      nodelist.push_back(i);
    }
   
    /*Loop over the index vector:
      we add all the nodes that satisfy the max distance requirement.
      We then remove the initial node from the list of nodes so that each pair
      is added only once*/
   
    vector< int >::iterator iter1;
    //first node loop
    for(iter1 = nodelist.begin(); iter1 != nodelist.end();){
      vector<int> neighbours;
      vector< int >::iterator iter2;
      //second node loop
      for(iter2 = nodelist.begin(); iter2 != nodelist.end(); iter2++){
	Real d2 = 0.0;
	Vector3d dx;
	for(int j = 0; j <_dim; j++){
	  dx(j) = (_x[*iter1 * _dim + j]-_x[*iter2 * _dim +j]);
	}
	// Periodic boudnary condition:
	_box->mapDistance(dx);
	d2 = dx.squaredNorm();
	//add to the list if satisfies length criterion:
	if(d2 < maxDist2 && d2> epsilon){
	  neighbours.push_back(*iter2);
	}
      }//end second node loop
      _neighboursMap.push_back(neighbours);
      //remove the node from the list and increment iterator:
      iter1 = nodelist.erase(iter1);
    }//end first node loop
       
  }

  void CrossLinkManager::addCrosslink( pair<int, int> nodePair){
    
    int n1 = nodePair.first;
    int n2 = nodePair.second;
    vector<int> nodesID;
    vector<Vector3d> NodeX;
    // Get the nodes ID for the cross link
    nodesID.push_back(n1);
    nodesID.push_back(n2);

    // Get the ref position for the 2 nodes
    Vector3d X1;
    X1 << _X(n1*_dim) , _X(n1*_dim+1), _X(n1*_dim+2);
    Vector3d X2;
    X2 << _X(n2*_dim) , _X(n2*_dim+1), _X(n2*_dim+2);

    NodeX.push_back(X1);
    NodeX.push_back(X2);
	
    // Construct cl element:
    GelElement* CLelem = new GelElement(0,nodesID, NodeX);
	
    // Construct cl object: gel element + cl material prop:
    CrossLink* CL = new CrossLink(CLelem, _clMat,_box);
    // add the crosslink to the model:
    _myGelModel->addCrosslink(CL);
    // add crosslink pointer to the crosslinked node map:
    _crossLinkMap[n1].first = true;
    _crossLinkMap[n2].first = true;
    _crossLinkMap[n1].second = CL;
    _crossLinkMap[n2].second = CL;
    
    
  }
  

  // KMC STUFF:

  void CrossLinkManager::KMC(int kmcIt, Real maxDist){
    
    Real k0 = 0.001;
    this->getCurrentConfig();
    cout << "finding neighbours..." <<endl;
    this->findNeighbours(maxDist);
    
    for (int it = 0; it < kmcIt ; it++){
      cout << "finding par candidates..." <<endl;
      this->findKMCPairsCandidate();
      cout << "done" <<endl;
      cout << "computing rates..." <<endl;
      this->computeRates(k0);
      cout << "done" << endl;
      int i = this->findEvent();
      cout << "event is " << i << endl;
      this->carryOutEvent(i);
    }
    
  }

  // Use kinetic Monte Carlo to choose event to carry:
  int CrossLinkManager::findEvent(){
    Real u = ((double) rand()/RAND_MAX);
    int numRates = _cumulativeRates.size();
    Real Q = _cumulativeRates.back();
    cout << "numrates is "<< numRates << " Qu is " << Q*u << endl;
    for(int i = 0; i < _cumulativeRates.size(); i++){
      if(_cumulativeRates[i] > u*Q){
	return i;
      }
    }
  }

  void CrossLinkManager::computeRates(Real k0){
    
    int numRates = _KMCPairsCandidates.size();

    _cumulativeRates.clear();
    _rates.clear();
    _cumulativeRates.resize(numRates,0.0);
    _rates.resize(numRates,0.0);

    cout << "number of possible events is " << numRates << endl;
    cout<< _cumulativeRates.size() << endl;
    for(int i = 0; i< numRates; i++ ){
      int n1 = _KMCPairsCandidates[i].first;
      int n2 = _KMCPairsCandidates[i].second;
      
      Vector3d x1;
      Vector3d x2;
      
      // Get the current position for the 2 nodes
      x1 << _x[n1*_dim] , _x[n1*_dim+1], _x[n1*_dim+2];
      x2 << _x[n2*_dim] , _x[n2*_dim+1], _x[n2*_dim+2];
	
      Vector3d AB = x2 - x1;
      _box->mapDistance(AB);
      vector<Vector3d> dx;
      dx.push_back(AB);
	
      // Compute the energy:
      FilamentMaterial::Filresults Rf;
      Rf.request = 3;
      _clMat->compute(Rf, dx);
    
      if(_crossLinkMap[n1].first == true & _crossLinkMap[n2].first == true){
      	_rates[i] = k0*exp(Rf.W);
      } else if(_crossLinkMap[n1].first == false & _crossLinkMap[n2].first == false){
      	_rates[i] = k0*exp(-Rf.W);
      }else{
	cout<< "ERROR1" <<endl;
	_rates[i] = 0.0;
      }
      
      for(int k = 0 ; k < i+1 ; k++){
	_cumulativeRates[i] += _rates[k];
      }
    }
    
  }
 
  void CrossLinkManager::carryOutEvent(int i){
    pair<int,int> nodePair = _KMCPairsCandidates[i];
    int n1 = nodePair.first;
    int n2 = nodePair.second;

    assert((_crossLinkMap[n1].first == true & _crossLinkMap[n2].first == true)\
	   ||(_crossLinkMap[n1].first == false & _crossLinkMap[n2].first == false));

    if(_crossLinkMap[n1].first == true & _crossLinkMap[n2].first == true
       & _crossLinkMap[n1].second == _crossLinkMap[n2].second){
      cout << "DELETING" <<endl;
      _myGelModel->deleteCrosslink(_crossLinkMap[n1].second);
      const vector<int> NodesID1 = _crossLinkMap[n1].second->getNodesID();
      const vector<int> NodesID2 = _crossLinkMap[n2].second->getNodesID();
     
      _crossLinkMap[n1].first = false;
      _crossLinkMap[n2].first = false;
      _crossLinkMap[n1].second = NULL;
      _crossLinkMap[n2].second = NULL;
      
    }else if(_crossLinkMap[n1].first == false & _crossLinkMap[n2].first == false ) {
      cout << "connecting" <<endl;

      this->addCrosslink(nodePair);

    }
  }
    

  void CrossLinkManager::findKMCPairsCandidate(){
    /* Build a list of node pairs to be considered by the KMC algorithm for
       crosslinking or detachement.*/
    _KMCPairsCandidates.clear();
    vector <int > addedCLnodes(0);

    for(int node1 = 0; node1 < _NumNode; node1++){
      vector< int > localNeighbours = _neighboursMap[node1];
      if(_crossLinkMap[node1].first == false){
	for (int j = 0; j < localNeighbours.size(); j++){
	  int node2 = localNeighbours[j];
	  if(_crossLinkMap[node2].first == false){
	    _KMCPairsCandidates.push_back(make_pair(node1,node2));
	  }
	}
      }else{
	//find the other node it is crosslinked to:
	const vector<int> nodesID = _crossLinkMap[node1].second->getNodesID();
	//sanitary check:
	assert(nodesID[0] == node1 | nodesID[1] == node1);
	if(find(addedCLnodes.begin(), addedCLnodes.end(), nodesID[0]) != addedCLnodes.end() \
	   && find(addedCLnodes.begin(), addedCLnodes.end(), nodesID[1]) != addedCLnodes.end()){
	  //do nothing, the CL has already been added
	}else{
	  _KMCPairsCandidates.push_back(make_pair(nodesID[0],nodesID[1]));
	  addedCLnodes.push_back(nodesID[0]);
	  addedCLnodes.push_back(nodesID[1]);
	}
      }
    }
  }
  

  

} // namespace voom
