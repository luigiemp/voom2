#include <vector>
#include <iostream>

#include "ScalarNode.h"
#include "MechanicsNode.h"
#include "VectorNode.h"

using namespace voom;

int main()
{
  cout << endl << "Testing voom node class ... " << endl;
  cout << "................................... " << endl << endl;

  // Scalar node
  int ID = 1;
  vector<int > DOF(3, 0); DOF[1] = 1; DOF[2] = 2;
  Vector3d X; X << 0.0, 1.0, -1.0;
  ScalarNode SN(ID, DOF, X);
  cout << "Scalar node testing" << endl;
  cout << "ID should be "  << ID  << " and is "  << SN.getId()  << endl;
  cout << "DOF should be " << DOF[0] << " " <<  DOF[1] << " " <<  DOF[2] << " and are " << SN.getDoF()[0] << " "  << SN.getDoF()[1] << " "  << SN.getDoF()[2] << endl;
  cout << "X should be "   << X   << " and is "  << SN.getX()   << endl;

  // Mechanics node
  ID = 2;
  DOF[1] = 3; DOF[2] = 5;
  X << 10.0, 9.0, -1.0;
  MechanicsNode MN(ID, DOF, X);
  cout << endl << endl;
  cout << "Mechanics node testing" << endl;
  cout << "ID should be "  << ID  << " and is "  << MN.getId()  << endl;
  cout << "DOF should be " << DOF[0] << " " <<  DOF[1] << " " <<  DOF[2] << " and are " << MN.getDoF()[0] << " "  << MN.getDoF()[1] << " "  << MN.getDoF()[2] << endl;
  cout << "X should be "   << X   << " and is "  << MN.getX()   << endl;
  cout << "Field should be " << X   << " and is "  << MN.getField()   << endl;
  Vector3d FieldMN; FieldMN << 1.0, 1.5, 2.0;
  MN.setField(FieldMN);
  cout << "Field should be " << FieldMN   << " and is "  << MN.getField()   << endl;

  // Vector node
  ID = 3;
  DOF[1] = 6; DOF[2] = 7;
  X << 1.0, 29.0, -5.0;
  VectorXd FieldVN(4); FieldVN << 1.0, 2.0, 3.0, 7.0;
  VectorNode VN(ID, DOF, X, FieldVN);
  cout << endl << endl;
  cout << "Mechanics node testing" << endl;
  cout << "ID should be "  << ID  << " and is "  << VN.getId()  << endl;
  cout << "DOF should be " << DOF[0] << " " <<  DOF[1] << " " <<  DOF[2] << " and are " << VN.getDoF()[0] << " "  << VN.getDoF()[1] << " "  << VN.getDoF()[2] << endl;
  cout << "X should be "   << X   << " and is "  << VN.getX()   << endl;
  cout << "Field should be " << FieldVN   << " and is "  << VN.getField()   << endl;

  return 0;
}
