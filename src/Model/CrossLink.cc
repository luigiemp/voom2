#include "CrossLink.h"


namespace voom {
  //Constructor
  CrossLink::CrossLink(GelElement* aCrossLink, FilamentMaterial* spring,PeriodicBox* box):
    _spring(spring),_myCrossLink(aCrossLink),_box(box)
  {
    _X = _myCrossLink->getNodesX();
    _nodeNum = _myCrossLink->getNodesPerElement();
    _NodesID = _myCrossLink->getNodesID();
    _dim = 3;
  }

  
  void CrossLink::compute(Result & R,vector<Vector3d> & xlist)
  {
    // Compute single crosslink stretching energy
    
    FilamentMaterial::Filresults Rf_stretch;
    Rf_stretch.request = R.getRequest();
    
    for(int n = 0; n < _nodeNum-1; n++)
      {
	//Nodes are named as follows (see spring material class):
	//nodeA : n
	//nodeB : n+1

	vector<Vector3d> dx;
	Vector3d AB = xlist[n+1] - xlist[n];
  
	// Perdiodic boundary conditions:
	_box->mapDistance(AB);

	dx.push_back(AB);
	// stretching energy
	_spring->compute(Rf_stretch,dx);
	
	for(int i=0; i<_dim; i++) {
	  R.addResidual(_NodesID[n]*_dim + i,Rf_stretch.f1(i));
	  R.addResidual(_NodesID[n+1]*_dim + i,-Rf_stretch.f1(i));
	}
	
	R.addEnergy(Rf_stretch.W);
	
      }
    
  }// end compute function
}





