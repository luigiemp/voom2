#include <vector>
#include <iostream>

#include "../MechanicsNode.h"

using namespace voom;

int main()
{
  cout << endl << "Testing voom node class ... " << endl;
  cout << "................................... " << endl << endl;

  VectorXd a(3); a << 1.0,0.0,0.0;
  MechanicsNode A(0,3, a); 
  cout << "Node ID  : " << A.getId() << endl;
  cout << "Node DOF : " << A.getDoF() << endl;
  vector<string> property = A.getPropertyList();
  Real *data = NULL; uint length;
  for(vector<string>::iterator it=property.begin();it != property.end();++it){
    cout << *it << endl;
    A.getProperty(*it, &data, length);
    for(unsigned int i = 0; i < length; i++) cout << data[i] << " ";
    cout << endl << endl;
  }

  // Linear Update
  cout << "Updaing all internal variables" << endl;
  VectorXd deltaU = VectorXd::Random(3);
  A.linearizedUpdate(deltaU);
  cout << property[1] << endl;
  A.getProperty(property[1], &data, length);
  for(unsigned int i = 0; i < length; i++) cout << data[i] << " ";
  cout << endl << endl;

  cout << "Updating only one internal variable" << endl;
  A.linearizedUpdate(2, -0.5);
  A.getProperty(property[1], &data, length);
  cout << property[1] << endl;
  for(unsigned int i = 0; i < length; i++) cout << data[i] << " ";
  cout << endl << endl;

  delete data;
  cout << endl << "........................... " << endl;
  cout << "Test of voom node class completed" << endl;
  return 0;
}
