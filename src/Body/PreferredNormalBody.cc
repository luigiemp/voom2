#include "PreferredNormalBody.h"

namespace voom {

  // Constructor
  PreferredNormalBody::PreferredNormalBody(Mesh* myMesh, State* myState, Real Krot):
    Body(myMesh, myState), _krot(Krot)
    {
      // Initialize _field vector for the part corresponding to this Body
      this->initializeField(1.0);
      this->initializeNormals();
    } // Constructor



  void PreferredNormalBody::initializeField(Real fact) { // This should not be necessary since PreferredNormalBody is always coupled to another mechanics body, which will be initialized
    vector<GeomElement* > Elements = _myMesh->getElements();
    for (int e = 0; e < _myMesh->getNumberOfElements(); e++) { // Loop Thhrough all elements - some nodes will be initiliazed multiple times - Otherwise we can first compile a set of unique nodal ids for this body
      vector<int > ElNodesID = Elements[e]->getNodesID();
      for (int n = 0; n < ElNodesID.size(); n++) {
	int nodeNum = ElNodesID[n];
	for (int i = 0; i < 3; i++) {
	  _myState->setPhi(nodeNum, i, _myState->getX(nodeNum)(i));
	}
      }
    }
  } // initializeField



  void PreferredNormalBody::initializeNormals() {
    
    // Compute the average normal in the reference configuration (_initNormals)
    vector<GeomElement*> elements = _myMesh->getElements();
    for(int e = 0; e < elements.size(); e++)
    {
      GeomElement* geomEl = elements[e];
      const vector<int  >& NodesID = geomEl->getNodesID();
      const int numQP    = geomEl->getNumberOfQuadPoints();
      const int numNodes = NodesID.size();

      // Loop over quadrature points
      for(int q = 0; q < numQP; q++) {

        // Compute normal based on current field - if used at the beginning, this is the same as using X
        Vector3d a1 = Vector3d::Zero(), a2 = Vector3d::Zero(), a3 = Vector3d::Zero();
        for (int a = 0; a < NodesID.size(); a++) {
          int nodeID = NodesID[a];
          Vector3d xa;
	  xa << _myState->getPhi(nodeID, 0), _myState->getPhi(nodeID, 1), _myState->getPhi(nodeID, 2);

          a1 += xa * geomEl->getDN(q, a, 0);
          a2 += xa * geomEl->getDN(q, a, 1);
        }
        a3 = a1.cross(a2);
	_initNormals.push_back(a3/a3.norm());
      } // end loop over quadrature points q
    } // end loop over elements e
  } // initializeNormals



  // Compute Function - Compute Energy, Force, Stiffness
  void PreferredNormalBody::compute(Result* R)
  {
    vector<GeomElement*> elements = _myMesh->getElements();
    
    for(int e = 0; e < elements.size(); e++) {
      GeomElement* geomEl = elements[e];
      const vector<int >& NodesID = geomEl->getNodesID();
      const int numQP    = geomEl->getNumberOfQuadPoints();
      const int numNodes = NodesID.size();

      // Loop over quadrature points
      for(int q = 0; q < numQP; q++) {
	Vector3d normalQ = _initNormals[e * numQP + q];
	// Compute the normal of the element at quadrature point
	Vector3d a1 = Vector3d::Zero(), a2 = Vector3d::Zero(), a3 = Vector3d::Zero();
        for (int a = 0; a < NodesID.size(); a++) {
          int nodeID = NodesID[a];
          Vector3d xa;
          xa << _myState->getPhi(nodeID, 0), _myState->getPhi(nodeID, 1), _myState->getPhi(nodeID, 2);
          a1 += xa * geomEl->getDN(q, a, 0);
          a2 += xa * geomEl->getDN(q, a, 1);
        }
        a3 = a1.cross(a2);
	Real a3norm = a3.norm();
	a3 = a3/a3norm;

	Vector3d normalDiff = a3 - normalQ;
	
	// Compute energy
	if (R->getRequest() & ENERGY) {
	  R->addEnergy(0.5 * _krot * square( normalDiff.norm() ) );
	}
    
	// Data structure used for force and stiffness
	vector<MatrixXd> dnedx(3, MatrixXd::Zero(numNodes, 3)); // Used for stiffness but computed here
	if (R->getRequest() & FORCE | STIFFNESS) {
	  for (int a = 0; a < numNodes; a++) {
	    for (int j = 0; j < 3; j++) {
	      int aj = NodesID[a] * 3 + j;
	      double tempResidual = 0.0;
	      for (int k = 0 ; k < 3; k++) {
		double temp1 = 0.0; double temp2 = 0.0;
		for (int m = 0; m < 3; m++) {
		  for (int b = 0; b < numNodes; b++) {
		    int bm = NodesID[b] * 3 + m;
		    temp1 += _myState->getPhi(bm) * (LeviCivita(k,j,m) * geomEl->getDN(q,a,0) * geomEl->getDN(q,b,1) + LeviCivita(k,m,j) * geomEl->getDN(q,b,0) * geomEl->getDN(q,a,1)) ;
		    for (int l = 0; l < 3; l++) {
		      temp2 += -1.0 * _myState->getPhi(bm) * a3(k) * a3(l) * (LeviCivita(l,j,m) * geomEl->getDN(q,a,0) * geomEl->getDN(q,b,1) + LeviCivita(l,m,j) * geomEl->getDN(q,b,0) * geomEl->getDN(q,a,1));
		    } // end l loop
		  } // end b loop
		} // end m loop
		dnedx[k](a,j) = 1.0/a3norm * (temp1 + temp2);
		tempResidual += dnedx[k](a,j) * normalDiff(k);
	      } // end k loop
	      if (R->getRequest() & FORCE) R->addResidual(aj, _krot * tempResidual);
	    } // end j loop
	  } // end a loop
	} // end FORCE | STIFFNESS request

	if (R->getRequest() & STIFFNESS) {
	  for (int a = 0; a < numNodes; a++) {
	    for (int j = 0; j < 3; j++) {
	      int aj = NodesID[a] * 3 + j;
	      for (int c = 0; c < numNodes; c++) {
		for (int n = 0; n < 3; n++) {
		  int cn = NodesID[c] * 3 + n;
		  double temp1 = 0.0;
		  double temp2 = 0.0;
		  double temp2_1a = 0.0;
		  double temp2_1b = 0.0;
		  double temp2_2 = 0.0;
		  double temp2_3 = 0.0;
		  for (int k = 0; k < 3; k++) {
		    temp1 += dnedx[k](c,n) * dnedx[k](a,j);
		    temp2_1a += normalDiff(k) * (LeviCivita(k,j,n) * geomEl->getDN(q,a,0) * geomEl->getDN(q,c,1) + LeviCivita(k,n,j) * geomEl->getDN(q,c,0) * geomEl->getDN(q,a,1));
		    for (int l = 0; l < 3; l++) {
		      temp2_1b -= normalDiff(k) * a3(k) * a3(l) * (LeviCivita(l,j,n) * geomEl->getDN(q,a,0) * geomEl->getDN(q,c,1) + LeviCivita(l,n,j) * geomEl->getDN(q,c,0) * geomEl->getDN(q,a,1));
		      for (int m = 0; m < 3; m++) {
			for (int b = 0; b < numNodes; b++) {
			  int bm = NodesID[b] * 3 + m;
			  temp2_2 += -1.0 * normalDiff(k) * _myState->getPhi(bm) * ((dnedx[k](c,n) * a3(l) + a3(k) * dnedx[l](c,n)) * (LeviCivita(l,j,m) * geomEl->getDN(q,a,0) * geomEl->getDN(q,b,1) + LeviCivita(l,m,j) * geomEl->getDN(q,b,0) * geomEl->getDN(q,a,1)));
			  temp2_3 += -1.0 * normalDiff(k) * _myState->getPhi(bm) * dnedx[k](a,j) * a3(l) * (LeviCivita(l,n,m) * geomEl->getDN(q,c,0) * geomEl->getDN(q,b,1) + LeviCivita(l,m,n) * geomEl->getDN(q,b,0) * geomEl->getDN(q,c,1));
			} // end b loop
		      } // end m loop
		    } // end l loop
		  } // end k loop
		  temp2 = temp2_1a + temp2_1b + temp2_2 + temp2_3;
		  temp2 *= 1.0/a3norm;
		  R->addStiffness(aj, cn, _krot*(temp1 + temp2) ); 
		} // end n loop
	      } // end c loop
	    } // end j loop
	  } // end a loop
	} // end STIFFNESS request
       
      } // end loop over quadrature points
    } // end loop over elements   
       

  } // Compute PreferredNormal Body

 

} // namespace voom

