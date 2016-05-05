#include "Filament.h"


namespace voom {
  //Constructor
  Filament::Filament(GelElement* aFilament, FilamentMaterial* spring, FilamentMaterial* angleSpring):
    _spring(spring),_angleSpring(angleSpring),_myFilament(aFilament)
  {
    _X = _myFilament->getNodesX();
    _nodeNum = _myFilament->getNodesPerElement();
    _NodesID = _myFilament->getNodesID();
    _dim = 3;
  }

  
  void Filament::compute(Result & R,vector<Vector3d> & xlist)
  {
    // Compute single filament stretching and bending energy
    
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
	_spring->compute(Rf_stretch,x,X);
	
	for(int i=0; i<_dim; i++) {
	  R.addResidual(_NodesID[n]*_dim + i,Rf_stretch.f1(i));
	  R.addResidual(_NodesID[n+1]*_dim + i,-Rf_stretch.f1(i));
	}
	
	R.addEnergy(Rf_stretch.W);
	
	if(n>0){
	  x.push_back(xlist[n-1]);
	  X.push_back(_X[n-1]);
	  
	  //bending energy
	  _angleSpring->compute(Rf_bend,x,X);

	  for(int i=0; i<_dim; i++) {
	    R.addResidual(_NodesID[n-1]*_dim + i,Rf_bend.f1(i));
	    R.addResidual(_NodesID[n]*_dim + i,-Rf_bend.f1(i)-Rf_bend.f2(i));
	    R.addResidual(_NodesID[n+1]*_dim + i,Rf_bend.f2(i));
	  }	  	  
	}

	R.addEnergy(Rf_bend.W);
	
      }
    
  }// end compute function
}





