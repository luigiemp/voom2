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
  {
    ThirdOrderTensor K(3,2,2);
    K.set(0,1,1, 5.0);
    cout << endl << "Entry at [0,1,1] is " << K.get(0,1,1) << endl;
    K.incrementIterator();
    K.sequentialSet(10.0);
    cout << "Entry at [1,0,0] is " << K.get(1,0,0) << endl;
    K.resetIterator();
    cout << "Printing the transpose " << endl;
    for (uint k=0; k<3; k++) {
      for (uint i=0; i<2; i++) {
	for (uint j=0; j<2; j++) {
	  cout << K.sequentialGet() << " ";
	  K.incrementIterator();
	}
	cout << endl;
      }
      cout << endl;
    }
  }

  {
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
  }

  // Testing matrix exponential
  {
    Matrix3d A;
    A << 0, -M_PI/4, 0,
       M_PI/4, 0,     0,
       0,    0,     0;
    std::cout << "The matrix A is:\n" << A << "\n\n";
    std::cout << "The matrix exponential of A is:\n" << A.exp() << "\n\n";

    // With symmetric matrices
     A << 0, M_PI/4, 0,
          M_PI/4, 0,     0,
          0,    0,     0;
    std::cout << "The matrix A is:\n" << A << "\n\n";
    std::cout << "The matrix exponential of A is:\n" << A.exp() << "\n\n";
    std::cout << "OR the matrix exponential of A is:\n" << VoomExpSymmMatrix(A) << "\n\n";

    // Correctness test
    for (int i = 0; i < 3; i++) {
      Matrix3d X = Matrix3d::Random();
      A = X + X.transpose();
      cout << "Eigen - Voom = " <<  A.exp() - VoomExpSymmMatrix(A) << endl;
    }
  
    // Speed test
    clock_t t;
    t = clock();
    for (int i = 0; i < 10000000; i++) {
      Matrix3d X = Matrix3d::Random();
      A = X + X.transpose();
      Matrix3d B = A.exp();
    }
    t = clock() - t;
    std::cout << "Time using Eigen exp = " << ((float)t)/CLOCKS_PER_SEC << endl;

    t = clock();
    for (int i = 0; i < 10000000; i++) {
      Matrix3d X = Matrix3d::Random();
      A = X + X.transpose();
      Matrix3d B = VoomExpSymmMatrix(A);
    }
    t = clock() - t;
    std::cout << "Time using VoomExpSymmMatrix = " << ((float)t)/CLOCKS_PER_SEC << endl;
  }

  // Testing Levi-Civita 3D Symbol
  {
    Matrix3d IdentityMat = Matrix3d::Identity();
    std::cout << "Testing Levi Civita 3D Symbol." << endl;
    std::cout << "Test 1. \\delta_{ij} \\epsilon_{ijk} = {0 0 0} " << endl;
    Vector3d test1(0.0, 0.0, 0.0);
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
	  test1(i) += IdentityMat(i,j) * LeviCivita(i,j,k);
	}
      }
    }
    if (test1.norm() == 0) cout << "Test 1 PASSED." << endl;
    else cout << "Test 1 FAILED." << endl;

    cout << endl << "Test 2. \\epsilon_{ipq} \\epsilon_{jpq} = 2 \\delta_{ij} " << endl;
    Matrix3d test2 = Matrix3d::Zero();
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
	for (int p = 0; p < 3; p++) {
	  for (int q = 0; q < 3; q++) {
	    test2(i,j) += LeviCivita(i,p,q) * LeviCivita(j,p,q);
	  }
	}
      }
    }
    test2 = test2 - 2 * IdentityMat;
    if (test2.norm() == 0) cout << "Test 2 PASSED." << endl;
    else cout << "Test 2 FAILED." << endl;

    cout << endl << "Test 3. \\epsilon_{ijk} \\epsilon_{ijk} = 6" << endl;
    int test3 = 0;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          test3 += LeviCivita(i,j,k) * LeviCivita(i,j,k);
        }
      }
    }
    if (test3 == 6) cout << "Test 3 PASSED." << endl;
    else cout << "Test 3 FAILED." << endl;
  }


  cout << endl << "....................... " << endl;
  cout << "Test of VoomMath class completed" << endl;
  return 0;
}
