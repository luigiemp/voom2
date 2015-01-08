#include "FEMesh.h"

namespace voom {

  // Constructor from input file
  FEMesh::FEMesh(const string inputFile): Mesh(inputFile) {
    ifstream inp(inputFile.c_str());
    string line;

    // Look for number of elements
    while( getline(inp, line) )
      if (line.find( "*NUMBEROFELEMENTS") == 0) break;
    getline(inp, line); // Read number of Elements
    const int nElements = atoi( line.c_str() );
    // Resize element container
    _elements.resize( nElements );

    // Look for element types
    while( getline(inp, line) ) 
      if ( line.find("*ELEMENT") == 0) break;
    // List of Element types in the given mesh
    // Map between element type and type ID
    map<string, int> elemMap; 
    // Element connectivity
    vector<vector<int> > connectivity(nElements);
    // Element type vector
    vector<string> eTypes(nElements);
    for(uint i = 0; i < nElements; i++) {
      getline(inp, line);
      vector<string> strs = splitString( line, " \t");
      if ( elemMap.find(strs[0]) == elemMap.end() )
	elemMap[strs[0]] = elemMap.size();
      eTypes[i] = strs[0];
      // Fill in connectivity table
      vector<int> conn(strs.size() - 2);
      for(uint m = 1; m < strs.size()-1; m++) 
	conn[m-1] = atoi(strs[m].c_str());
      connectivity[i] = conn;
    } // loop over nElements

    // Created Required shapes and quadrature rules
    _quadrature.resize( elemMap.size() );
    // Loop through element map to find all the element types present in the mesh
    for(map<string, int>::iterator it = elemMap.begin(); it != elemMap.end(); it++) {
      if (it->first == "C3D8" ) {
	// Full integration hexahedral element
	HexQuadrature *qhex = new HexQuadrature(2); 
	vector<VectorXd > QuadPoints = qhex->getQuadPoints();
	uint NumQP = QuadPoints.size();
	vector<Shape* > elemShape(NumQP, NULL);
	for (uint s=0; s<NumQP; s++)
	  elemShape[s] = new HexShape(QuadPoints[s]);

	_quadrature[it->second] = qhex;
	_shapes.push_back(elemShape);
      } 
      else if ( it->first == "C3D8R") {
	// Reduced integration hexadehral element
	HexQuadrature *qhex = new HexQuadrature(1); 
	vector<VectorXd > QuadPoints = qhex->getQuadPoints();
	uint NumQP = QuadPoints.size();
	vector<Shape* > elemShape(NumQP, NULL);
	for (uint s=0; s<NumQP; s++)
	  elemShape[s] = new HexShape(QuadPoints[s]);

	_quadrature[it->second] = qhex;
	_shapes.push_back(elemShape);
      }
      else {
	cerr << "** ERROR: Unknown finite element " << it->first << endl;
	cerr << "Exiting...\n";
	exit(EXIT_FAILURE);
      }
    } // Loop over elemMap - different element in the mesh

    // Create list of elements
    for(uint i = 0; i < nElements; i++) {
      int lookUp = elemMap[ eTypes[i] ];
      vector< VectorXd > positionList( connectivity[i].size() );
      // Temporary list of postion of nodes which make the element
      for(uint m = 0; m < connectivity[i].size(); m++) 
	positionList[m] = _positions[ connectivity[i][m] ];

      FEgeomElement *myElement = new FEgeomElement(i, connectivity[i], 
						   positionList,
						   _shapes[lookUp],
						   _quadrature[lookUp]);
      _elements[i] = myElement;
    } // Loop over element list 

  } // Constructor from input file



  // Constructor from nodes and connectivities
  FEMesh::FEMesh(const vector<VectorXd > &  Positions,
		 const vector<vector<int > > & Connectivity, 
		 const vector<int > & LocalDoF,
		 const vector<int > & GhostDoF,
		 string ElementType,
		 uint QuadOrder): 
    Mesh(Positions, LocalDoF, GhostDoF)
  {
    // Created required shapes and quadrature rules
    if (ElementType == "Hexa") {
	// Hexahedral element
	HexQuadrature *qhex = new HexQuadrature(QuadOrder); 
	vector<VectorXd > QuadPoints = qhex->getQuadPoints();
	uint NumQP = QuadPoints.size();
	vector<Shape* > elemShape(NumQP, NULL);
	for (uint s = 0; s < NumQP; s++)
	  elemShape[s] = new HexShape(QuadPoints[s]);

	_quadrature.push_back(qhex);
	_shapes.push_back(elemShape);
    } 
    else {
      cerr << "** ERROR: Unknown finite element " <<  ElementType << endl;
      cerr << "Exiting...\n";
      exit(EXIT_FAILURE);
    }

    // Resize element container
    uint NumEl = Connectivity.size();
    _elements.resize(NumEl);
    
    // Create list of elements
    for(uint i = 0; i < NumEl; i++) {
      // Temporary list of nodal positions
      vector<VectorXd > positionList(Connectivity[i].size() );
      for(uint m = 0; m < Connectivity[i].size(); m++) 
	positionList[m] = _positions[Connectivity[i][m]];
      
      FEgeomElement *myElement = new FEgeomElement(i, Connectivity[i], 
						   positionList,
						   _shapes[0],
						   _quadrature[0]);
      _elements[i] = myElement;
    } // Loop over element list 

  } // Constructor from nodes and connectivities
  
} // namespace voom
