#include "Mesh.h"
#include "FEMesh.h"

namespace voom
{
  // Constructor. Creates position table only.
  Mesh::Mesh(const string Nodes, const string ConnTable) {
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
      // Real mag = 0;
      // inp >> temp; - if this is commented each line is just the [x,y,z] ccordinate, no node number
      for (uint j = 0; j < dim; j++) {
	inp >> _X[i](j);
	// mag += _X[i](j)*_X[i](j);
      }
      // mag = sqrt(mag); // To get an accurate sphere
      // mag = 1.0; //For any other surface
      // _X[i](0) = _X[i](0)/mag; _X[i](1) = _X[i](1)/mag; _X[i](2) = _X[i](2)/mag;
      //if (i<=12) _X[i] = (1-.000125)*_X[i];
    }

    inp.close();
  }
  
  // // Static function
  // Mesh* Mesh::New(const string inputFile) {
  //   Mesh* myMesh;
  //   // Parse inputfile to see element type and decide which element is used
  //   ifstream inp(inputFile.c_str());
  //   string line;
  //   if (!inp.is_open()){
  //     cerr << "** ERROR: Opening input file " << inputFile << endl;
  //     cerr << "Exiting...\n";
  //     exit(EXIT_FAILURE);
  //   }
  //   // Parse till we find element card
  //   while( getline(inp, line) ) 
  //     if ( line.find("*ELEMENT") == 0) break;
  //   // Read all element types in the model. Store in a set
  //   set<string> elementTypes;
  //   while( getline(inp, line) ) {
  //     if ( line.find("**") == 0) continue;
  //     if ( line.find("*") == 0) continue;
  //     vector<string> strs = splitString( line, " \t");
  //     elementTypes.insert( strs[0] );
  //   }
  //   inp.close();
  //   // Check element types to decide mesh type
  //   MESHTYPE type;
  //   for(set<string>::iterator it = elementTypes.begin();
  // 	it != elementTypes.end(); it++ ) 
  //     if ( it->find("MFRK") == 0) type = MESHFREERKPM;
  //     else if ( it->find("MFLM") == 0 ) type = MESHFREELME;
  //     else type = FINITEELEMENT;
    
  //   // Create derived class onject based on mesh type
  //   switch( type ) {
  //   case 0:
  //     {
  // 	// Finite element mesh
  // 	myMesh = new FEMesh(inputFile);
  // 	break;
  //     }
  //   case 1:
  //     {
  // 	// Meshfree RKPM
  // 	cout << "Meshfree RKPM not yet implemented\n";
  // 	break;
  //     }
  //   case 2:
  //     {
  // 	cout << "Meshfree LME not yet implemented\n";
  // 	// Meshfree LME
  // 	break;
  //     }
  //   default:
  //     {
  // 	cerr << "** ERROR: Unknown mesh type\n";
  // 	cerr << "Exiting...\n";
  // 	exit(EXIT_FAILURE);
  //     }
  //   } // Switch statement
  //   return myMesh;
  // }

} //namespace voom

