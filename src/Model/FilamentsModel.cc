#include "FilamentModel.h"


namespace voom {
  //Constructor
  FilamentModel::FilamentModel(GeomFilament* aFilament, FilamentMaterial* material,
			       const uint NodeDoF,
			       Vector3d d0, Real k);
  Model(NodeDoF),_material(material),_nodalForcesFlag(NodalForcesFlag),_resetFlag(resetFlag), _d0(d0),_k(k)
    ,_myFilament(aFilament)
  {
    _material->setMaterialParameters(_k);
    _material->setInternalParameters(_d0);
    
  }

  void compute(Result & R,vector<Vector3d> & dlist)
  {
    const uint nodeNum = NodesID.size();
    vector<Vector3d > dlist(nodeNum-1,Vector3d::Zero());

    FilamentMaterial::Filresults Rf;
    Rf.request = 7;

    for(int n = 0; n<nodeNum-2; n++)
      {

	vector<int> bond;
	bond.push_back(NodesID[n]);
	bond.push_back(NodesID[n+1]);
	Vector3d d0;
	cout << "OK" << endl;
	Vector3d x1;
	Vector3d x2;

	x1 << getX(bond[0],0) , getX(bond[0],1) , getX(bond[0],2);
	x2 << getX(bond[1],0) , getX(bond[1],1) , getX(bond[1],2);

	d0 = x2-x1;

	_material->compute(Rf,dlist[n]);
	

      }
  }


  }

  
}


