#include "CrossLink.h"


namespace voom {
  //Constructor
  CrossLink::CrossLink(GelElement* aCrossLink, FilamentMaterial* spring):
    _spring(spring),_myCrossLink(aCrossLink)
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
    FilamentMaterial::Filresults Rf_bend;
    Rf_stretch.request = R.getRequest();
    Rf_bend.request = R.getRequest();

    for(int n = 0; n < _nodeNum-1; n++)
      {
	vector<Vector3d> x;
	vector<Vector3d> X;

	x.push_back(xlist[n]);
	x.push_back(xlist[n+1]);
	X.push_back(_X[n]);
	X.push_back(_X[n+1]);
	
	// stretching energy
	_spring->compute(Rf_stretch,x);
	
	for(int i=0; i<_dim; i++) {
	  R.addResidual(_NodesID[n]*_dim + i,Rf_stretch.f1(i));
	  R.addResidual(_NodesID[n+1]*_dim + i,-Rf_stretch.f1(i));
	}
	
	R.addEnergy(Rf_stretch.W);
	
      }
    
  }// end compute function
}





