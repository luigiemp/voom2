#include "../VoomMath.h"

using namespace voom;

int main()
{
  cout << endl << "Testing VoomMath class ... " << endl;
  cout << ".......................... " << endl << endl;
  

  Real r = 3.0; 
  cout << r << " square    = " << square(r)         << endl << endl;
  cout << r << " factorial = " << factorial(int(r)) << endl << endl;

  VectorXd a( VectorXd::Zero(2)), b(VectorXd::Zero(2));
  a(0) = 1.0; a(1) = 2.0;
  b(0) =-1.0; b(1) = 3.0;
  cout << a << "\n dot \n" << b   << "\n = \n" << a.dot(b) << endl << endl;

  MatrixXd Ma( MatrixXd::Zero(2,2)), Mb( MatrixXd::Zero(2,2)), 
    Mc( MatrixXd::Zero(2,2));
  Ma << 1.0,-3.0,
        2.0, 5.0;
  Mb << 0.0, 1.0,
        1.0, 3.0;
  
  b = Ma*a;
  cout << Ma << "\n dot \n" << a  << "\n = \n" << b << endl << endl;

  Mc = Ma*Mb;
  cout << Ma << "\n dot \n" << Mb << "\n = \n" << Mc << endl << endl;
  
  cout << Ma << " \n determinant = \n " << det(Ma) << endl << endl;
  
  inv(Ma, Mc);
  cout << Ma <<  " \n inverse \n  = " << Mc << endl << endl;

  cout << "Ma*inv(Ma) = \n " << Ma*Mc << endl;


  FourthOrderTensor K(3,3,3,3);
  K.set(0,1,2,0, 5.0);
  cout << endl << "Entry at [0,1,2,0] is " << K.get(0,1,2,0) << endl;
  K.incrementIterator();
  K.sequentialSet(10.0);
  cout << "Entry at [1,0,0,0] is " << K.get(1,0,0,0) << endl;
  K.resetIterator();
  cout << "Printing the transpose " << endl;
  for (uint k=0; k<3; k++) {
    for (uint i=0; i<3; i++) {
      for (uint j=0; j<3; j++) {
	cout << K.sequentialGet() << " ";
	K.incrementIterator();
      }
      cout << endl;
    }
    cout << endl;
  }
  
  cout << endl << "....................... " << endl;
  cout << "Test of VoomMath class completed" << endl;
  return 0;
}
