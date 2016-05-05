#include "GelMesh.h"
#include <iostream>

namespace voom {
  
  // Constructor from input file
  GelMesh::GelMesh(const string Nodes, const string ConnTable){
     
    ifstream inp;
    inp.open(Nodes.c_str(), ios::in);
    if (!inp) {
      cout << "Cannot open input file: " << Nodes << endl;
      exit(0);
    }
    
    uint NumNodes = 0, dim = 0, temp = 0;
    // First number is NumNodes, Second number is dim                                                         
    inp >> NumNodes >> dim;
    _X.resize(NumNodes, VectorXd::Zero(dim));

    // All other lines are the nodal coordinates
    for (uint i = 0; i < NumNodes; i++) {
                
      for (uint j = 0; j < dim; j++) {
        inp >> _X[i](j);
      }
    }

    inp.close();
    
    
    // Second we load and read the connectivity table: 1 element
    // per line with each node ID in order
    ifstream inpCon;
    inpCon.open(ConnTable.c_str(),ios::in);
    
    if (!inpCon) {
      cout << "Cannot open input file: " << ConnTable << endl;
      exit(0);
    }
    
    string line;
    uint NumEl = 0;
    if (inpCon.is_open())
      {
	while (getline (inpCon,line))
	  {
	    stringstream stream(line);
	    int n;
	    vector<int> fil;
	    vector<Vector3d> NodesX;
	    
	    while (stream >> n){
	      fil.push_back(n);
	      Vector3d NodeX;
	      NodeX << _X[n](0) , _X[n](1), _X[n](2);
	      NodesX.push_back(NodeX);
	    }
	    _elements.push_back(new GelElement(NumEl,fil,NodesX));
	    NumEl++;
	  }
	inpCon.close();
      }
    
   
  } // End constructor from input files
} // namespace voom
