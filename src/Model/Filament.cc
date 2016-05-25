#include "Filament.h"


namespace voom {
  //Constructor
  Filament::Filament(GelElement* aFilament, FilamentMaterial* spring, FilamentMaterial* angleSpring, PeriodicBox* box):
    _spring(spring),_angleSpring(angleSpring),_myFilament(aFilament), _box(box)
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
	//Nodes are named as follows (see spring material class):
	//nodeA : n
	//nodeB : n+1

	vector<Vector3d> dx;
	Vector3d AB = xlist[n+1] - xlist[n];
	
	// Perdiodic boundary conditions:
	_box->mapDistance(AB);

	dx.push_back(AB);

	// //periodic BC test:
	// Vector3d A;
	// Vector3d B;
	// A<< 0.0,0.0,0.0;
	// B<<1.3,0.0,0.0;
	// Vector3d ab = B-A;
	// _box->mapDistance(ab);
	// cout << ab << endl;
	
	// stretching energy
	_spring->compute(Rf_stretch,dx);
	
	for(int i=0; i<_dim; i++) {
	  R.addResidual(_NodesID[n]*_dim + i,Rf_stretch.f1(i));
	  R.addResidual(_NodesID[n+1]*_dim + i,-Rf_stretch.f1(i));
	}
	
	R.addEnergy(Rf_stretch.W);
	
	//if more than 2 nodes, we compute the bending energy:
	if(n>0){
	  //In this case, nodes are named as follows (see anglespring mat class):
	  //nodeA : n-1
	  //nodeB : n
	  //nodeC : n+1
	  // so dx[0] is nodeC-nodeB = tBC, already in dx
	  // and tBA = nodeA-nodeB is added to dx:
	  Vector3d BA;
	  BA = xlist[n-1]-xlist[n];
	  _box->mapDistance(BA);
	  dx.push_back(BA);

	  //bending energy
	  _angleSpring->compute(Rf_bend,dx);

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





