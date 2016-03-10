#include "FEMesh.h"

namespace voom {

  // Constructor from input file
  FEMesh::FEMesh(const string Nodes, const string ConnTable): Mesh(Nodes, ConnTable) {
    ifstream inp(ConnTable.c_str());
    
    Real radius = 1.0;
    uint NumEl = 0, temp = 0;
    string ElType;
    // First line
    // First number is NumNodes, Second is type of element
    inp >> NumEl >> ElType;
    _elements.resize(NumEl);
    
    uint NumNodesEl = 0;
    if (ElType == "C3D8") { 
      // Full integration hexahedral element
      _quadrature = new HexQuadrature(2); 
      vector<VectorXd > QuadPoints = _quadrature->getQuadPoints();

      // Fill quadrature and shapes
      for ( uint q = 0; q < QuadPoints.size(); q++ ) {
	_shapes.push_back( new HexShape(QuadPoints[q]) );
      }
      NumNodesEl = 8;
    } // end of C3D8
    else if (ElType == "L2") {
      // Line integral with Legendre-Gauss quadrature
      _quadrature = new LineQuadrature(2);
      vector<VectorXd > QuadPoints = _quadrature->getQuadPoints();
      
      // Fill quadrature and shapes    
      for ( uint q = 0; q < QuadPoints.size(); q++ ) {
        _shapes.push_back( new BarShape(QuadPoints[q]) );
      }
      NumNodesEl = 2;
    }
    else if (ElType == "C3D8R") {
      // Reduced integration hexahedral element
      _quadrature = new HexQuadrature(1); 
      vector<VectorXd > QuadPoints = _quadrature->getQuadPoints();
      
      // Fill quadrature and shapes
      for ( uint q = 0; q < QuadPoints.size(); q++ ) {
	_shapes.push_back( new HexShape(QuadPoints[q]) );
      }
      NumNodesEl = 8;
    }
    else if (ElType == "C3D4") {
      // Full integration linear tetrahedral element
      _quadrature = new TetQuadrature(1); 
      vector<VectorXd > QuadPoints = _quadrature->getQuadPoints();
      
      // Fill quadrature and shapes
      for ( uint q = 0; q < QuadPoints.size(); q++ ) {
	_shapes.push_back( new LinTetShape(QuadPoints[q]) );
      }
      NumNodesEl = 4;
    }
    else if (ElType == "C3D10") {
      // Full integration quadratic tetrahedral element
      _quadrature = new TetQuadrature(2); 
      vector<VectorXd > QuadPoints = _quadrature->getQuadPoints();
      
      // Fill quadrature and shapes
      for ( uint q = 0; q < QuadPoints.size(); q++ ) {
	_shapes.push_back( new QuadTetShape(QuadPoints[q]) );
      }
      NumNodesEl = 10;
    }
    else if (ElType == "TD3") {
      // Full integration linear triangular element
      _quadrature = new TriQuadrature(1); 
      vector<VectorXd > QuadPoints = _quadrature->getQuadPoints();
      
      // Fill quadrature and shapes
      for ( uint q = 0; q < QuadPoints.size(); q++ ) {
	_shapes.push_back( new LinTriShape(QuadPoints[q]) );
      }
      NumNodesEl = 3;
    }
    else if (ElType == "TD6") {
      // Full integration quadratic triangular element
      _quadrature = new TriQuadrature(2); 
      vector<VectorXd > QuadPoints = _quadrature->getQuadPoints();
      
      // Fill quadrature and shapes
      for ( uint q = 0; q < QuadPoints.size(); q++ ) {
	_shapes.push_back( new QuadTriShape(QuadPoints[q]) );
      }
      NumNodesEl = 6;
    }
    else {
      cerr << "** ERROR: Unknown finite element type: " << ElType << endl;
      cerr << "Exiting...\n";
      exit(EXIT_FAILURE);
    }
        
    if (ElType == "L2") {
      
      for (uint e = 0; e < NumEl; e++) {
	vector<int > ConnEl(NumNodesEl, 0);
	vector<VectorXd > Xel;
	for (uint n = 0; n < NumNodesEl; n++) {
	  inp >> ConnEl[n];
	  Xel.push_back(_X[ConnEl[n]]);
	  
	}
	
	
	_elements[e] = new FEgeomElement1D(e, ConnEl,
					   Xel,
					   _shapes,
					   _quadrature,
					   radius);
	
      } // End of FEgeom1D
      
    }                
    else {
      // Compute the geometric elements                                                                                                                                                                                                         
      for (uint e = 0; e < NumEl; e++) {
	vector<int > ConnEl(NumNodesEl, 0);
	vector<VectorXd > Xel;
	for (uint n = 0; n < NumNodesEl; n++) {
	  inp >> ConnEl[n];
	  Xel.push_back(_X[ConnEl[n]]);
	}

	_elements[e] = new FEgeomElement(e, ConnEl,
					 Xel,
					 _shapes,
					 _quadrature);
      } // End of FEgeom 


    }
    
  } // End constructor from input files



  // // Constructor from nodes and connectivities
  // FEMesh::FEMesh(const vector<VectorXd > &  Positions,
  // 		 const vector<vector<int > > & Connectivity, 
  // 		 string ElementType,
  // 		 uint QuadOrder): 
  //   Mesh(Positions)
  // {
  //   // Created required shapes and quadrature rules
  //   if (ElementType == "Hexa") {
  // 	// Hexahedral element
  // 	HexQuadrature *qhex = new HexQuadrature(QuadOrder); 
  // 	vector<VectorXd > QuadPoints = qhex->getQuadPoints();
  // 	uint NumQP = QuadPoints.size();
  // 	vector<Shape* > elemShape(NumQP, NULL);
  // 	for (uint s = 0; s < NumQP; s++)
  // 	  elemShape[s] = new HexShape(QuadPoints[s]);

  // 	_quadrature.push_back(qhex);
  // 	_shapes.push_back(elemShape);
  //   } 
  //   else {
  //     cerr << "** ERROR: Unknown finite element " <<  ElementType << endl;
  //     cerr << "Exiting...\n";
  //     exit(EXIT_FAILURE);
  //   }

  //   // Resize element container
  //   uint NumEl = Connectivity.size();
  //   _elements.resize(NumEl);
    
  //   // Create list of elements
  //   for(uint i = 0; i < NumEl; i++) {
  //     // Temporary list of nodal positions
  //     vector<VectorXd > positionList(Connectivity[i].size() );
  //     for(uint m = 0; m < Connectivity[i].size(); m++) 
  // 	positionList[m] = _positions[Connectivity[i][m]];
      
  //     FEgeomElement *myElement = new FEgeomElement(i, Connectivity[i], 
  // 						   positionList,
  // 						   _shapes[0],
  // 						   _quadrature[0]);
  //     _elements[i] = myElement;
  //   } // Loop over element list 

  // } // Constructor from nodes and connectivities
  
} // namespace voom
