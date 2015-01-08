#include "LMEgeomElement.h"
#include "LMEShape.h"

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
  // Create a Mesh Free Element with 9 nodes
  int elemID = 0;            // Element ID
  vector<int > nodesID(9,0); // Connectivity table
  for (uint i=0; i<9; i++)
    nodesID[i] = i;

  // Nodes
  vector<VectorXd > nodesX;
  srand( time(NULL) );
  Vector3d X; 
  X << 0.0, 0.0, 0.0; X = X + getRand();
  nodesX.push_back(X);
  X << 1.0, 0.0, 0.0; X = X + getRand();
  nodesX.push_back(X);
  X << 1.0, 1.0, 0.0; X = X + getRand();
  nodesX.push_back(X);
  X << 0.0, 1.0, 0.0; X = X + getRand();
  nodesX.push_back(X);
  X << 0.0, 0.0, 1.0; X = X + getRand();
  nodesX.push_back(X);
  X << 1.0, 0.0, 1.0; X = X + getRand();
  nodesX.push_back(X);
  X << 1.0, 1.0, 1.0; X = X + getRand();
  nodesX.push_back(X);
  X << 0.0, 1.0, 1.0; X = X + getRand();
  nodesX.push_back(X);
  X << 0.5, 0.5, 0.5; X = X + getRand();
  nodesX.push_back(X);

  // Quadrature points
  vector<VectorXd > MaterialPoints(2, Vector3d::Zero(3));
  MaterialPoints[0] << 0.4, 0.4, 0.4;
  MaterialPoints[1] << 0.3, 0.5, 0.6;

  // Volume per quadrature point
  vector<Real > Weights(2, 1.0);
  Weights[1] += 0.2;

  // Create LME Geom Element
  Real beta = 0.8, tol = 1.0e-12;
  uint MaxIter = 50;
  LMEgeomElement LMEelement(elemID, nodesID, nodesX, MaterialPoints, Weights, beta, tol, MaxIter);

  // Print element members
  cout << "ElemID             = " << LMEelement.getGeomElementID()   << endl;
  cout << "NodesNumber        = " << LMEelement.getNodesPerElement() << endl;
  cout << "Connectivity table = " ;
  vector<int > ConnTable = LMEelement.getNodesID();
  for (uint a = 0; a < 9; a++)
    cout << ConnTable[a] << " ";
  cout << endl;

  for (uint q = 0; q < MaterialPoints.size(); q++) {
    cout << "QP weight = " << endl ;
    cout << LMEelement.getQPweights(q) << " ";
    cout << endl;

    cout << "N = " << endl ;
    for (uint a = 0; a < 9; a++)
      cout << LMEelement.getN(q, a) << " ";
    cout << endl;

    cout << "DN = " << endl;
    for (uint a = 0; a < 9; a++) {
      for (uint i = 0; i < 3; i++)
	cout << LMEelement.getDN(q, a, i) << " ";
      cout << endl;
    }
    cout << endl;
  } // q loop

  cout << "All done " << endl;
} // main
