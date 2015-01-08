#include "Mesh.h"
#include "FEMesh.h"

namespace voom
{
  // Constructor. Creates position table only.
  Mesh::Mesh(const string inputFile) {
    ifstream inp(inputFile.c_str());
    string line;
    while( getline(inp, line) )
      if ( line.find("*NUMBEROFNODES") == 0) break;
    getline(inp, line); // Read number of Nodes
    const int nNodes = atoi( line.c_str() );
    _positions.resize( nNodes );
    while( getline(inp, line) )
      if ( line.find("*NODE") == 0) break;

    for(uint nodeId = 0; nodeId < nNodes; nodeId++) {    
      getline(inp, line);
      vector<string> strs = splitString( line, " \t");
      VectorXd data( strs.size() );
      for(int i = 0; i < strs.size(); i++) data(i) = atof( strs[i].c_str() );
      _positions[nodeId] = data;
    }
    inp.close();
  }
  
  // Static function
  Mesh* Mesh::New(const string inputFile) {
    Mesh* myMesh;
    // Parse inputfile to see element type and decide which element is used
    ifstream inp(inputFile.c_str());
    string line;
    if (!inp.is_open()){
      cerr << "** ERROR: Opening input file " << inputFile << endl;
      cerr << "Exiting...\n";
      exit(EXIT_FAILURE);
    }
    // Parse till we find element card
    while( getline(inp, line) ) 
      if ( line.find("*ELEMENT") == 0) break;
    // Read all element types in the model. Store in a set
    set<string> elementTypes;
    while( getline(inp, line) ) {
      if ( line.find("**") == 0) continue;
      if ( line.find("*") == 0) continue;
      vector<string> strs = splitString( line, " \t");
      elementTypes.insert( strs[0] );
    }
    inp.close();
    // Check element types to decide mesh type
    MESHTYPE type;
    for(set<string>::iterator it = elementTypes.begin();
	it != elementTypes.end(); it++ ) 
      if ( it->find("MFRK") == 0) type = MESHFREERKPM;
      else if ( it->find("MFLM") == 0 ) type = MESHFREELME;
      else type = FINITEELEMENT;
    
    // Create derived class onject based on mesh type
    switch( type ) {
    case 0:
      {
	// Finite element mesh
	myMesh = new FEMesh(inputFile);
	break;
      }
    case 1:
      {
	// Meshfree RKPM
	cout << "Meshfree RKPM not yet implemented\n";
	break;
      }
    case 2:
      {
	cout << "Meshfree LME not yet implemented\n";
	// Meshfree LME
	break;
      }
    default:
      {
	cerr << "** ERROR: Unknown mesh type\n";
	cerr << "Exiting...\n";
	exit(EXIT_FAILURE);
      }
    } // Switch statement
    return myMesh;
  }

} //namespace voom

