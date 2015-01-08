#include "RKPMgeomElement.h"
#include "MRKPMShape.h"

using namespace voom;

Vector2d getRand() {
  Vector2d rand2D;
    rand2D << 0.2*(Real(rand())/RAND_MAX - 0.5), 
      0.3*(Real(rand())/RAND_MAX - 0.5);
  return rand2D;
}

int main() 
{
  // Create a Mesh Free Element with 9 nodes
  int elemID = 0;            // Element ID
  vector<int > nodesID(4,0); // Connectivity table
  for (uint i=0; i<4; i++)
    nodesID[i] = i;

  // Nodes
  vector<VectorXd > nodesX;
  srand( time(NULL) );
  Vector2d X; 
  X << 0.0, 0.0; X = X + getRand();
  nodesX.push_back(X);
  X << 1.0, 0.0; X = X + getRand();
  nodesX.push_back(X);
  X << 1.0, 1.0; X = X + getRand();
  nodesX.push_back(X);
  X << 0.0, 1.0; X = X + getRand();
  nodesX.push_back(X);

  // Quadrature points
  vector<VectorXd > MaterialPoints(2, Vector2d::Zero(2));
  MaterialPoints[0] << 0.5, 0.5;
  MaterialPoints[1] << 0.3, 0.5;

  // Volume per quadrature point
  vector<Real > Weights(2, 1.0);

  // Create RKPM Geom Element
  const Real support = 1.414, supportHat = 1.0, radius = 0.001;
  RKPMgeomElement RKPMelement(elemID, nodesID, nodesX, MaterialPoints,
			      Weights, support, radius, supportHat);

  // Print element members
  cout << "ElemID             = " << RKPMelement.getGeomElementID()   << endl;
  cout << "NodesNumber        = " << RKPMelement.getNodesPerElement() << endl;
  cout << "Connectivity table = " ;
  vector<int > ConnTable = RKPMelement.getNodesID();
  for (uint a = 0; a < 4; a++)
    cout << ConnTable[a] << " ";
  cout << endl;

  for (uint q = 0; q < MaterialPoints.size(); q++) {
    cout << "QP weight = " << endl ;
    cout << RKPMelement.getQPweights(q) << " ";
    cout << endl;

    cout << "N = " << endl ;
    for (uint a = 0; a < 4; a++)
      cout << RKPMelement.getN(q, a) << " ";
    cout << endl;

    cout << "DN = " << endl;
    for (uint a = 0; a < 4; a++) {
      for (uint i = 0; i < 2; i++)
	cout << RKPMelement.getDN(q, a, i) << " ";
      cout << endl;
    }
    cout << endl;
  } // q loop

  cout << "All done " << endl;
} // main
