#include "FEMesh.h"

namespace voom {

  // Constructor from input file
  FEMesh::FEMesh(const string Nodes, const string ConnTable, State* myState, int dofPerNode, bool CheckOverlap)
    : Mesh(Nodes, ConnTable, myState, dofPerNode, CheckOverlap) {
    ifstream inp(ConnTable.c_str());

    int NumEl = 0, temp = 0;
    string ElType;
    // First line
    // First number is NumNodes, Second is type of element
    inp >> NumEl >> ElType;
    _elementType = ElType;
    _elements.resize(NumEl);

    int NumNodesEl = this->createElementShapeAndQuadrature(ElType);
    if (NumNodesEl == 0) {
      cout << "Failed reading " << ConnTable << endl;
      exit(EXIT_FAILURE);
    }

    // Compute the geometric elements
    for (int e = 0; e < NumEl; e++) {
      vector<int > ConnEl(NumNodesEl, 0);
      vector<VectorXd > Xel;
      for (int n = 0; n < NumNodesEl; n++) {
        inp >> temp;
	ConnEl[n] = _LocaltoGlobal[temp];
        Xel.push_back(_gState->getX(ConnEl[n]));
      }
      _elements[e] = new FEgeomElement(e, ConnEl, Xel, _shapes, _quadrature);
    } // End of FEgeom
  } // End constructor from input files



  // Constructor from nodes and connectivities
  FEMesh::FEMesh(const vector<VectorXd > &  Positions, State* myState, int dofPerNode, bool CheckOverlap,
    const vector<vector<int > > & Connectivity, string ElementType):
    Mesh(Positions, myState, dofPerNode, CheckOverlap)
  {
    _elementType = ElementType;
    int NumNodesEl = this->createElementShapeAndQuadrature(ElementType);

    // Resize element container
    int NumEl = Connectivity.size();
    _elements.resize(NumEl);

    // Create list of elements
    for(int e = 0; e < NumEl; e++) {
      // Temporary list of nodal positions
      vector<VectorXd > Xel;
      vector<int > ConnEl(NumNodesEl, 0);
      for(int n = 0; n < Connectivity[e].size(); n++) {
	ConnEl[n] = _LocaltoGlobal[Connectivity[e][n]];
        Xel.push_back(_gState->getX(ConnEl[n]));
      }

      _elements[e] = new FEgeomElement(e, ConnEl, Xel, _shapes, _quadrature);
    } // Loop over element list
  } // Constructor from nodes and connectivities



  int FEMesh::createElementShapeAndQuadrature(const string ElType) {
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
    else if (ElType == "Q4") {
      // Full integration linear quadrilateral element
      _quadrature = new QuadQuadrature(1);
      vector<VectorXd> QuadPoints = _quadrature->getQuadPoints();

      // Fill quadrature and shapes
      for ( uint q = 0; q < QuadPoints.size(); q++ ) {
        _shapes.push_back( new LinQuadShape(QuadPoints[q]) );
      }
      NumNodesEl = 4;
    }
    else {
      cerr << "** ERROR: Unknown finite element type: " << ElType << endl;
      cerr << "Exiting...\n";
    }
    return NumNodesEl;
  } // End createElementShapeAndQuadrature
} // namespace voom
