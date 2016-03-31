#include "GelMesh.h"
#include <iostream>

namespace voom {
  
  // Constructor from input file
  GelMesh::GelMesh(const string Nodes, const string ConnTable){
    // First we load and read the connectivity table: 1 filament per line with each node ID in order
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
	    while (stream >> n){
	      fil.push_back(n);
	    }
	    _filaments.push_back(new GeomFilament(NumEl,fil));
	    NumEl++;
	  }
	inpCon.close();
      }
    
        
    ifstream inp;
    inp.open(Nodes.c_str(), ios::in);
    if (!inp) {
      cout << "Cannot open input file: " << Nodes << endl;
      exit(0);
    }
    
    uint NumNodes = 0, dim = 0, temp = 0;
    // First line                                                                                                                                       
    // First number is NumNodes, Second numbder is dim
    inp >> NumNodes >> dim;
    _X.resize(NumNodes, VectorXd::Zero(dim));

    // All other lines are the nodal coordinates
    for (uint i = 0; i < NumNodes; i++) {
      // inp >> temp; - if this is commented each line is just the [x,y,z] ccordinate, no node number
      for (uint j = 0; j < dim; j++) {
        inp >> _X[i](j);
      }
    }

    inp.close();
    
  } // End constructor from input files
} // namespace voom
