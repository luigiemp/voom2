#include "FEgeomElement1D.h"
#include "BarShape.h"
#include "LineQuadrature.h"

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
  // Create a 1D Element 
  unsigned int id = 0;
  const Real radius = 1.0;
  vector<int > nodesID(2, 0.0);
  nodesID[1] = 1;
  vector<Vector3d > nodesX(2, Vector3d::Zero());
  nodesX[0] += getRand();
  nodesX[1](0) = 1.0;
  nodesX[1] += getRand();

  // Bar Shape and Quadrature
  LineQuadrature QuadRule(2);
  const vector<VectorXd > QuadPoints = QuadRule.getQuadPoints();
  vector<Shape *> shapePointers;
  for (uint q = 0; q < QuadPoints.size(); q++) {
    Shape * shapeP = new BarShape(QuadPoints[q]);
    shapePointers.push_back(shapeP);
  }
  
  // Create Geom Element
  FEgeomElement1D BarElement(id, nodesID, nodesX, shapePointers, &QuadRule, radius);  
  
  // Print element members
  cout << "ElemID             = " << BarElement.getGeomElementID()   << endl;
  cout << "NodesNumber        = " << BarElement.getNodesPerElement() << endl;
  cout << "Connectivity table = " ;
  vector<int > ConnTable = BarElement.getNodesID();
  for (int a = 0; a < 2; a++)
    cout << ConnTable[a] << " ";
  cout << endl;

  for (uint q = 0; q < QuadPoints.size(); q++)
  {
    cout << "QP weight = " << endl ;
    cout << BarElement.getQPweights(q) << " ";
    cout << endl;

    cout << "N = " << endl;
    for (int a = 0; a < 2; a++)
      cout << BarElement.getN(q, a) << " ";
    cout << endl;

    cout << "DN = " << endl ;
    for (int a = 0; a < 2; a++) {
      for (int i = 0; i < 3; i++)
	cout << BarElement.getDN(q, a, i) << " ";
      cout << endl;
    }
    cout << endl;
  } // q loop

  for (uint q = 0; q < QuadPoints.size(); q++) {
    delete shapePointers[q];
  }

}
