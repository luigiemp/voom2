#include "FEgeomElement.h"
#include "HexShape.h"
#include "HexQuadrature.h"
#include "LinTetShape.h"
#include "QuadTetShape.h"
#include "TetQuadrature.h"
#include "QuadTriShape.h"
#include "TriQuadrature.h"

using namespace voom;

Vector3d getRand() {
  Vector3d rand3D;
    rand3D << 0.12*(Real(rand())/RAND_MAX - 0.5), 
              0.18*(Real(rand())/RAND_MAX - 0.5), 
              0.06*(Real(rand())/RAND_MAX - 0.5);
  return rand3D;
}

Vector2d getRand2d() {
  Vector2d rand2D;
    rand2D << 0.05*(Real(rand())/RAND_MAX - 0.5), 
              0.1*(Real(rand())/RAND_MAX - 0.5);
  return rand2D;
}

int main() 
{
  {
    cout << "Inspecting Hex Element" << endl;
    // Create a Hex Element
    int id = 0;                // Element ID
    vector<int > nodesID(8,0); // Connectivity table
    for (uint i=0; i<8; i++)
      nodesID[i] = i;
    
    // Nodes
    vector<VectorXd > nodesX;
    srand( time(NULL) );
    Vector3d X; 
    X << 0.0, 0.0, 0.0; X = X + getRand();
    nodesX.push_back(X);
    X << 1., 0., 0.; X = X + getRand();
    nodesX.push_back(X);
    X << 1., 1., 0.; X = X + getRand();
    nodesX.push_back(X);
    X << 0., 1., 0.; X = X + getRand();
    nodesX.push_back(X);
    X << 0., 0., 1.; X = X + getRand();
    nodesX.push_back(X);
    X << 1., 0., 1.; X = X + getRand();
    nodesX.push_back(X);
    X << 1., 1., 1.; X = X + getRand();
    nodesX.push_back(X);
    X << 0., 1., 1.; X = X + getRand();
    nodesX.push_back(X);
    
    // Hex Shape and Quadrature
    HexQuadrature QuadRule(1);
    const vector<VectorXd > QuadPoints = QuadRule.getQuadPoints();
    
    vector<Shape *> shapePointers;
    for (uint q = 0; q < QuadPoints.size(); q++) {
      Shape * shapeP = new HexShape(QuadPoints[q]);
      shapePointers.push_back(shapeP);
    }
    
    // Create Geom Element
    FEgeomElement HexElement(id, nodesID, nodesX, shapePointers, &QuadRule);
    
    // Print element members
    cout << "ElemID             = " << HexElement.getGeomElementID()   << endl;
    cout << "NodesNumber        = " << HexElement.getNodesPerElement() << endl;
    cout << "Connectivity table = " ;
    vector<int > ConnTable = HexElement.getNodesID();
    for (uint a = 0; a < 8; a++)
      cout << ConnTable[a] << " ";
    cout << endl;
    
    for (uint q = 0; q < QuadPoints.size(); q++) {
      cout << "QP weight = " << endl ;
      cout << HexElement.getQPweights(q) << " ";
      cout << endl;
      
      cout << "N = " << endl ;
      for (uint a = 0; a < 8; a++)
	cout << HexElement.getN(q, a) << " ";
      cout << endl;
      
      cout << "DN = " << endl;
      for (uint a = 0; a < 8; a++) {
	for (uint i = 0; i < 3; i++)
	  cout << HexElement.getDN(q, a, i) << " ";
	cout << endl;
      }
      cout << endl;
    } // q loop

    // Delete shape function objects
    for (uint q = 0; q < QuadPoints.size(); q++) {
      delete shapePointers[q];
    }
    
  }// End of test for hexa element



  {
    cout << "Inspecting Lin Tet Element" << endl;
    // Create a Lin Tet Element
    int id = 0;                // Element ID
    vector<int > nodesID(4,0); // Connectivity table
    for (uint i=0; i<4; i++)
      nodesID[i] = i;
    
    // Nodes
    vector<VectorXd > nodesX;
    srand( time(NULL) );
    Vector3d X; 
    X << 0., 0., 1.; X = X + getRand();
    nodesX.push_back(X);
    X << 1., 0., 0.; X = X + getRand();
    nodesX.push_back(X);
    X << 0., 1., 0.; X = X + getRand();
    nodesX.push_back(X);
    X << 0., 0., 0.; X = X + getRand();
    nodesX.push_back(X);
    
    // Tet Shape and Quadrature
    TetQuadrature QuadRule(1);
    const vector<VectorXd > QuadPoints = QuadRule.getQuadPoints();
    
    vector<Shape *> shapePointers;
    for (uint q = 0; q < QuadPoints.size(); q++) {
      Shape * shapeP = new LinTetShape(QuadPoints[q]);
      shapePointers.push_back(shapeP);
    }
    
    // Create Geom Element
    FEgeomElement TetElement(id, nodesID, nodesX, shapePointers, &QuadRule);
    
    // Print element members
    cout << "ElemID             = " << TetElement.getGeomElementID()   << endl;
    cout << "NodesNumber        = " << TetElement.getNodesPerElement() << endl;
    cout << "Connectivity table = " ;
    vector<int > ConnTable = TetElement.getNodesID();
    for (uint a = 0; a < 4; a++)
      cout << ConnTable[a] << " ";
    cout << endl;
    
    for (uint q = 0; q < QuadPoints.size(); q++) {
      cout << "QP weight = " << endl ;
      cout << TetElement.getQPweights(q) << " ";
      cout << endl;
      
      cout << "N = " << endl ;
      for (uint a = 0; a < 4; a++)
	cout << TetElement.getN(q, a) << " ";
      cout << endl;
      
      cout << "DN = " << endl;
      for (uint a = 0; a < 4; a++) {
	for (uint i = 0; i < 3; i++)
	  cout << TetElement.getDN(q, a, i) << " ";
	cout << endl;
      }
      cout << endl;
    } // q loop

    // Delete shape function objects
    for (uint q = 0; q < QuadPoints.size(); q++) {
      delete shapePointers[q];
    }
    
  }// End of test for linear tet element



 {
    cout << "Inspecting Quad Tet Element" << endl;
    // Create a Quad Tet Element
    int id = 0;                // Element ID
    vector<int > nodesID(10,0); // Connectivity table
    for (uint i=0; i<10; i++)
      nodesID[i] = i;
    
    // Nodes
    vector<VectorXd > nodesX;
    srand( time(NULL) );
    Vector3d X; 
    X << 0., 0., 0.; X = X + getRand();
    nodesX.push_back(X);
    X << 1., 0., 0.; X = X + getRand();
    nodesX.push_back(X);
    X << 0., 1., 0.; X = X + getRand();
    nodesX.push_back(X);
    X << 0., 0., 1.; X = X + getRand();
    nodesX.push_back(X);
    X << 0.5, 0., 0.; X = X + getRand();
    nodesX.push_back(X);
    X << 0.5, 0.5, 0.; X = X + getRand();
    nodesX.push_back(X);
    X << 0., 0.5, 0.0; X = X + getRand();
    nodesX.push_back(X);
    X << 0., 0., 0.5; X = X + getRand();
    nodesX.push_back(X);
    X << 0.5, 0., 0.5; X = X + getRand();
    nodesX.push_back(X);
    X << 0., 0.5, 0.5; X = X + getRand();
    nodesX.push_back(X);
    
    // Tet Shape and Quadrature
    TetQuadrature QuadRule(1);
    const vector<VectorXd > QuadPoints = QuadRule.getQuadPoints();
    
    vector<Shape *> shapePointers;
    for (uint q = 0; q < QuadPoints.size(); q++) {
      Shape * shapeP = new QuadTetShape(QuadPoints[q]);
      shapePointers.push_back(shapeP);
    }
    
    // Create Geom Element
    FEgeomElement TetElement(id, nodesID, nodesX, shapePointers, &QuadRule);
    
    // Print element members
    cout << "ElemID             = " << TetElement.getGeomElementID()   << endl;
    cout << "NodesNumber        = " << TetElement.getNodesPerElement() << endl;
    cout << "Connectivity table = " ;
    vector<int > ConnTable = TetElement.getNodesID();
    for (uint a = 0; a < 10; a++)
      cout << ConnTable[a] << " ";
    cout << endl;
    
    for (uint q = 0; q < QuadPoints.size(); q++) {
      cout << "QP weight = " << endl ;
      cout << TetElement.getQPweights(q) << " ";
      cout << endl;
      
      cout << "N = " << endl ;
      for (uint a = 0; a < 10; a++)
	cout << TetElement.getN(q, a) << " ";
      cout << endl;
      
      cout << "DN = " << endl;
      for (uint a = 0; a < 10; a++) {
	for (uint i = 0; i < 3; i++)
	  cout << TetElement.getDN(q, a, i) << " ";
	cout << endl;
      }
      cout << endl;
    } // q loop

    // Delete shape function objects
    for (uint q = 0; q < QuadPoints.size(); q++) {
      delete shapePointers[q];
    }
    
  }// End of test for quadratic tet element



 {
    cout << "Inspecting Quad Tri Element" << endl;
    // Create a Quad Tri Element
    int id = 0;                // Element ID
    vector<int > nodesID(6,0); // Connectivity table
    for (uint i=0; i<6; i++)
      nodesID[i] = i;
    
    // Nodes
    vector<VectorXd > nodesX;
    srand( time(NULL) );
    Vector2d X; 
    X << 0., 0.; X = X + getRand2d();
    nodesX.push_back(X);
    X << 1., 0.; X = X + getRand2d();
    nodesX.push_back(X);
    X << 0., 1.; X = X + getRand2d();
    nodesX.push_back(X);
    X << 0.5, 0.; X = X + getRand2d();
    nodesX.push_back(X);
    X << 0.5, 0.5; X = X + getRand2d();
    nodesX.push_back(X);
    X << 0., 0.5; X = X + getRand2d();
    nodesX.push_back(X);
    
    // Tri Shape and Quadrature
    TriQuadrature QuadRule(1);
    const vector<VectorXd > QuadPoints = QuadRule.getQuadPoints();
    
    vector<Shape *> shapePointers;
    for (uint q = 0; q < QuadPoints.size(); q++) {
      Shape * shapeP = new QuadTriShape(QuadPoints[q]);
      shapePointers.push_back(shapeP);
    }
    
    // Create Geom Element
    FEgeomElement TriElement(id, nodesID, nodesX, shapePointers, &QuadRule);
    
    // Print element members
    cout << "ElemID             = " << TriElement.getGeomElementID()   << endl;
    cout << "NodesNumber        = " << TriElement.getNodesPerElement() << endl;
    cout << "Connectivity table = " ;
    vector<int > ConnTable = TriElement.getNodesID();
    for (uint a = 0; a < 6; a++)
      cout << ConnTable[a] << " ";
    cout << endl;
    
    for (uint q = 0; q < QuadPoints.size(); q++) {
      cout << "QP weight = " << endl ;
      cout << TriElement.getQPweights(q) << " ";
      cout << endl;
      
      cout << "N = " << endl ;
      for (uint a = 0; a < 6; a++)
	cout << TriElement.getN(q, a) << " ";
      cout << endl;
      
      cout << "DN = " << endl;
      for (uint a = 0; a < 6; a++) {
	for (uint i = 0; i < 2; i++)
	  cout << TriElement.getDN(q, a, i) << " ";
	cout << endl;
      }
      cout << endl;
    } // q loop

    // Delete shape function objects
    for (uint q = 0; q < QuadPoints.size(); q++) {
      delete shapePointers[q];
    }
    
  }// End of test for linear tet element

}
