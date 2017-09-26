#include "PressureBody.h"

namespace voom {

  // Constructor
  PressureBody::PressureBody(Mesh* myMesh, State* myState, Real Pressure):
    Body(myMesh, myState), _pressure(Pressure)
    {
      // Initialize _field vector for the part corresponding to this Body
      this->initializeField(1);
      // Initialize _prevPhi to initial nodal positions
      this->setPrevField();
    } // Constructor



  void PressureBody::initializeField(Real fact) {
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



  // Compute Function - Compute Energy, Force, Stiffness
  void PressureBody::compute(Result* R)
  {
    const vector<GeomElement* > elements = _myMesh->getElements();

    // Loop through elements
    for(int e = 0; e < elements.size(); e++)
    {
      GeomElement* geomEl = elements[e];
      const vector<int  >& NodesID = geomEl->getNodesID();
      const int numQP     = geomEl->getNumberOfQuadPoints();
      const int numNodes  = NodesID.size();

      // Loop over quadrature points
      for(int q = 0; q < numQP; q++) {

        // Compute normal based on _prevPhi and displacement
        Vector3d a1 = Vector3d::Zero(), a2 = Vector3d::Zero(), a3 = Vector3d::Zero(), u = Vector3d::Zero();
        for (int a = 0; a < NodesID.size(); a++) {
          Vector3d xa_prev, xa_curr;
	  int NodeIDGdof = _myState-> getGdof(NodesID[a]);
          xa_prev << _prevPhi[NodeIDGdof], _prevPhi[NodeIDGdof + 1], _prevPhi[NodeIDGdof + 2];
          xa_curr << _myState->getPhi(NodeIDGdof), _myState->getPhi(NodeIDGdof + 1), _myState->getPhi(NodeIDGdof + 2);
          a1 += xa_prev*geomEl->getDN(q, a, 0);
          a2 += xa_prev*geomEl->getDN(q, a, 1);

          u += ( xa_curr - _myState->getX(NodesID[a]) )*geomEl->getN(q, a);
        }
        a3 = a1.cross(a2);

        // Surface associated with QP q
        Real Area = a3.norm();
        a3 /= Area;
        Area *= geomEl->getQPweights(q);
        // cout << "Area = " << Area << endl;
	// cout << a3 << endl;

        // Compute energy
        if (R->getRequest() & ENERGY) {
          R->addEnergy( _pressure*Area*a3.dot(u) );
        }

        // Compute Residual
        if ( (R->getRequest() & FORCE) || (R->getRequest() & DMATPROP) ) {
          for(uint a = 0; a < NodesID.size(); a++) {
            for(uint i = 0; i < 3; i++) {
              R->addResidual(_myState->getGdof(NodesID[a])+i, _pressure * Area * a3(i) * geomEl->getN(q, a));
              // cout << NodesID[a] << " " << geomEl->getN(q, a) << endl;
            } // i loop
          } // a loop
        } // Internal force loop

      } // loop over QP
    } // loop over elements

  } // Compute Pressure Body

 

} // namespace voom
