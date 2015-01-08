#include "FEgeomElement.h"
#include "HexShape.h"
#include "HexQuadrature.h"

using namespace voom;

Vector3d getRand() {
  Vector3d rand3D;
    rand3D << 0.2*(Real(rand())/RAND_MAX - 0.5), 
              0.3*(Real(rand())/RAND_MAX - 0.5), 
              0.4*(Real(rand())/RAND_MAX - 0.5);
  return rand3D;
}

int main() 
{
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

}
