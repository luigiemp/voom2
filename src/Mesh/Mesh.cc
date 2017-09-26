#include "Mesh.h"

namespace voom
{
  Mesh::Mesh(string Nodes, string ConnTable, State* myState, int dofPerNode, bool CheckOverlap)
    : _gState(myState) {
    ifstream inp;
    inp.open(Nodes.c_str(), ios::in);
    if (!inp) {
      cout << "Cannot open input file: " << Nodes << endl;
      exit(0);
    }

    int NumNodes = 0, dim = 0;
    // First line
    // First number is NumNodes, Second numbder is dim
    inp >> NumNodes >> dim;
    int StateSize = _gState->getXsize(); // Assume new mesh will not have overlapping nodes to be collapsed
    Real tol = 1.0e-10;
   
    int DuplicatedNodeCount = 0;
    _LocaltoGlobal.assign(NumNodes, -1);
    // All other lines are the nodal coordinates
    for (int i = 0; i < NumNodes; i++) {
      // Load current node coordinates
      Vector3d TempX;
      for (int j = 0; j < dim; j++) {
	inp >> TempX(j);
      }
      if (CheckOverlap && StateSize > 0) { // Check if node is already present in State
	int Inserted = 0;
	for (int k = 0; k < StateSize; k ++) {
	  Vector3d Diff = TempX - _gState->getX(k);
	  if (Diff.norm() < tol) {
	    _LocaltoGlobal[i] = k;
	    Inserted = 1;
	    DuplicatedNodeCount++;
	    // cout << " Local to Global (True) " << i << " " << k << endl;
	    // cout << TempX << endl;
	    // cout << _gState->getX(k) << endl;
	    break;
	  } 
	} // Loop over all nodes already in State
	if (Inserted == 0) {
	  _gState->insertNode(TempX, dofPerNode);
	  _LocaltoGlobal[i] = i + StateSize - DuplicatedNodeCount;
	  // cout << " Local to Global (True) " << i << " " << i+StateSize-DuplicatedNodeCount << endl;
	}
      } else {
	_gState->insertNode(TempX, dofPerNode);
	_LocaltoGlobal[i] = i + StateSize;
	// cout << " Local to Global (False) " << i << " " << i+StateSize << endl;
      }
    } // Loop over all new nodes
    // for (int i = 0; i< NumNodes; i++)
    //   cout << " Local to Global " << _LocaltoGlobal[i] << endl;
    // cout << endl << endl;
    
    inp.close();
  }



  Mesh::Mesh(const vector<VectorXd > &  Positions, State* myState, int dofPerNode, bool CheckOverlap)
    : _gState(myState) {
    
    int NumNodes = Positions.size(), dim = 3;
    int StateSize = _gState->getXsize(); // Assume new mesh will not have overlapping nodes to be collapsed
    Real tol = 1.0e-10;
   
    int DuplicatedNodeCount = 0;
    _LocaltoGlobal.assign(NumNodes, -1);
    // All other lines are the nodal coordinates
    for (int i = 0; i < NumNodes; i++) {
      // Load current node coordinates
      Vector3d TempX = Positions[i];
      
      if (CheckOverlap && StateSize > 0) { // Check if node is already present in State
	int Inserted = 0;
	for (int k = 0; k < StateSize; k ++) {
	  Vector3d Diff = TempX - _gState->getX(k);
	  if (Diff.norm() < tol) {
	    _LocaltoGlobal[i] = k;
	    Inserted = 1;
	    DuplicatedNodeCount++;
	    break;
	  } 
	} // Loop over all nodes already in State
	if (Inserted == 0) {
	  _gState->insertNode(TempX, dofPerNode);
	  _LocaltoGlobal[i] = i + StateSize - DuplicatedNodeCount;
	}
      } else {
	_gState->insertNode(TempX, dofPerNode);
	_LocaltoGlobal[i] = i + StateSize;
      }
    } // Loop over all new nodes

  }



} //namespace voom

