#include "VoomMath.h"

namespace voom {
      	
  //! recursive funtion to compute factorial
  unsigned int factorial(unsigned int n)
  {
    unsigned int temp;
    if(n <= 1) return 1;
    temp = n * factorial(n - 1);
    return temp;
  };



  Real det(const MatrixXd& A)
  {
    if (A.size() == 4)
      return ( A(0,0)*A(1,1) - A(0,1)*A(1,0) );
    else if (A.size() == 9)
      return ( A(0,0)*A(1,1)*A(2,2) +
	       A(1,0)*A(2,1)*A(0,2) +
	       A(0,1)*A(1,2)*A(2,0) -
	       A(0,2)*A(1,1)*A(2,0) -
	       A(0,1)*A(1,0)*A(2,2) -
	       A(1,2)*A(2,1)*A(0,0) );
    else if (A.size() == 16)
      return ( A(0,0)*A(1,1)*A(2,2)*A(3,3) 
	       - A(0,0)*A(1,1)*A(2,3)*A(3,2) 
	       - A(0,0)*A(1,2)*A(2,1)*A(3,3) 
	       + A(0,0)*A(1,2)*A(2,3)*A(3,1) 
	       + A(0,0)*A(1,3)*A(2,1)*A(3,2) 
	       - A(0,0)*A(1,3)*A(2,2)*A(3,1) 
	       - A(0,1)*A(1,0)*A(2,2)*A(3,3) 
	       + A(0,1)*A(1,0)*A(2,3)*A(3,2) 
	       + A(0,1)*A(1,2)*A(2,0)*A(3,3) 
	       - A(0,1)*A(1,2)*A(2,3)*A(3,0) 
	       - A(0,1)*A(1,3)*A(2,0)*A(3,2) 
	       + A(0,1)*A(1,3)*A(2,2)*A(3,0) 
	       + A(0,2)*A(1,0)*A(2,1)*A(3,3) 
	       - A(0,2)*A(1,0)*A(2,3)*A(3,1) 
	       - A(0,2)*A(1,1)*A(2,0)*A(3,3) 
	       + A(0,2)*A(1,1)*A(2,3)*A(3,0) 
	       + A(0,2)*A(1,3)*A(2,0)*A(3,1) 
	       - A(0,2)*A(1,3)*A(2,1)*A(3,0)
	       - A(0,3)*A(1,0)*A(2,1)*A(3,2)  
	       + A(0,3)*A(1,0)*A(2,2)*A(3,1) 
	       + A(0,3)*A(1,1)*A(2,0)*A(3,2) 
	       - A(0,3)*A(1,1)*A(2,2)*A(3,0) 
	       - A(0,3)*A(1,2)*A(2,0)*A(3,1) 
	       + A(0,3)*A(1,2)*A(2,1)*A(3,0) );

    else {
      cerr << "** Matrix dimension not 2x2, 3x3 or 4x4. Cannot compute "
	   << "determinant" << endl;
      return -1;
    }
  };



  void inv(const  MatrixXd & A, MatrixXd & B)
  {
    const double detA = det(A);
    if ( A.size() == 4) {
      B(0,0) =  A(1,1)/detA;
      B(0,1) = -A(0,1)/detA;
      B(1,0) = -A(1,0)/detA;
      B(1,1) =  A(0,0)/detA;
    } else if (A.size() == 9 ) {
      B(0,0) = ( A(1,1)*A(2,2) - A(2,1)*A(1,2) )/detA;
      B(1,0) = ( A(2,0)*A(1,2) - A(1,0)*A(2,2) )/detA;
      B(2,0) = ( A(1,0)*A(2,1) - A(1,1)*A(2,0) )/detA;
      B(0,1) = ( A(2,1)*A(0,2) - A(0,1)*A(2,2) )/detA;
      B(1,1) = ( A(0,0)*A(2,2) - A(2,0)*A(0,2) )/detA;
      B(2,1) = ( A(0,1)*A(2,0) - A(0,0)*A(2,1) )/detA;
      B(0,2) = ( A(0,1)*A(1,2) - A(1,1)*A(0,2) )/detA;
      B(1,2) = ( A(0,2)*A(1,0) - A(0,0)*A(1,2) )/detA;
      B(2,2) = ( A(0,0)*A(1,1) - A(1,0)*A(0,1) )/detA;
    } else if (A.size() == 16 ) {
      B(0,0) = A(1,1)* A(2,2)* A(3,3) - A(1,1) *A(2,3)* A(3,2) 
	- A(1,2)* A(2,1)* A(3,3) + A(1,2)* A(2,3)* A(3,1) 
	+ A(1,3)* A(2,1)* A(3,2) - A(1,3)* A(2,2)* A(3,1);
      B(0,1) = A(0,1)* A(2,3)* A(3,2) - A(0,1) *A(2,2)* A(3,3)
	+ A(0,2)* A(2,1)* A(3,3) - A(0,2)* A(2,3) *A(3,1) 
	- A(0,3)* A(2,1)* A(3,2) + A(0,3) *A(2,2)* A(3,1);
      B(0,2) = A(0,1)* A(1,2)* A(3,3) 
	- A(0,1)* A(1,3) *A(3,2) - A(0,2)* A(1,1)* A(3,3) 
	+ A(0,2) *A(1,3) *A(3,1) + A(0,3)* A(1,1) *A(3,2) 
	- A(0,3)* A(1,2)* A(3,1);
      B(0,3) = A(0,1)* A(1,3)* A(2,2) - A(0,1)* A(1,2) *A(2,3) 
	+ A(0,2) *A(1,1)* A(2,3) 
	- A(0,2)* A(1,3) *A(2,1) - A(0,3)* A(1,1) *A(2,2) 
	+ A(0,3)* A(1,2) *A(2,1);
      B(1,0) = A(1,0) *A(2,3)* A(3,2) - A(1,0) *A(2,2)* A(3,3) 
	+ A(1,2)* A(2,0)* A(3,3) 
	- A(1,2)* A(2,3)*A(3,0) - A(1,3)* A(2,0)* A(3,2) 
	+ A(1,3)* A(2,2)* A(3,0);
      B(1,1) = A(0,0) *A(2,2)* A(3,3) - A(0,0)* A(2,3)* A(3,2) 
	- A(0,2)* A(2,0)*A(3,3) + A(0,2)* A(2,3)* A(3,0) 
	+ A(0,3)* A(2,0)* A(3,2) - A(0,3)* A(2,2)*A(3,0);
      B(1,2) = A(0,0)* A(1,3) *A(3,2) - A(0,0)* A(1,2)* A(3,3) 
	+ A(0,2)* A(1,0)*A(3,3) - A(0,2) *A(1,3) *A(3,0) 
	- A(0,3)* A(1,0) *A(3,2) + A(0,3)* A(1,2)*A(3,0);
      B(1,3) = A(0,0)* A(1,2)* A(2,3) - A(0,0) *A(1,3)* A(2,2) 
	- A(0,2) *A(1,0)*A(2,3) + A(0,2)* A(1,3)* A(2,0) 
	+ A(0,3)* A(1,0)* A(2,2) - A(0,3) *A(1,2)*A(2,0);
      B(2,0) = A(1,0) *A(2,1)* A(3,3) - A(1,0)* A(2,3)* A(3,1) 
	- A(1,1)* A(2,0)* A(3,3) + A(1,1)* A(2,3)* A(3,0) 
	+ A(1,3)* A(2,0) *A(3,1) - A(1,3)* A(2,1)*A(3,0);
      B(2,1) = A(0,0) *A(2,3)* A(3,1) - A(0,0)* A(2,1) *A(3,3) 
	+ A(0,1)* A(2,0)*A(3,3) - A(0,1) *A(2,3) *A(3,0) 
	- A(0,3)* A(2,0)* A(3,1) + A(0,3) *A(2,1)*A(3,0);
      B(2,2) =  A(0,0)* A(1,1) *A(3,3) - A(0,0) *A(1,3) *A(3,1) 
	- A(0,1) *A(1,0)*A(3,3) + A(0,1) *A(1,3)* A(3,0) 
	+ A(0,3) *A(1,0) *A(3,1) - A(0,3) *A(1,1)*A(3,0);
      B(2,3)= A(0,0) *A(1,3)* A(2,1) - A(0,0) *A(1,1) *A(2,3) 
	+ A(0,1) *A(1,0)*A(2,3) - A(0,1) *A(1,3) *A(2,0) 
	- A(0,3) *A(1,0) *A(2,1) + A(0,3) *A(1,1)*A(2,0);
      B(3,0) = A(1,0) *A(2,2) *A(3,1) - A(1,0) *A(2,1)* A(3,2) 
	+ A(1,1)* A(2,0) *A(3,2) - A(1,1) *A(2,2)* A(3,0) 
	- A(1,2)* A(2,0) *A(3,1) + A(1,2)* A(2,1)*A(3,0);
      B(3,1) = A(0,0)* A(2,1)* A(3,2) - A(0,0) *A(2,2)* A(3,1) 
	- A(0,1)* A(2,0)*A(3,2) + A(0,1)* A(2,2)* A(3,0) 
	+ A(0,2)* A(2,0)* A(3,1) - A(0,2) *A(2,1)*A(3,0);
      B(3,2) = A(0,0) *A(1,2) *A(3,1) - A(0,0)* A(1,1)* A(3,2) 
	+ A(0,1)* A(1,0)*A(3,2) - A(0,1)* A(1,2)* A(3,0) 
	- A(0,2)* A(1,0)* A(3,1) + A(0,2)* A(1,1)*A(3,0);
      B(3,3) = A(0,0) *A(1,1) *A(2,2) - A(0,0)* A(1,2)* A(2,1) 
	- A(0,1)* A(1,0)*A(2,2) + A(0,1) *A(1,2)* A(2,0) 
	+ A(0,2) *A(1,0)*A(2,1) - A(0,2)*A(1,1)*A(2,0) ;
    B = B/detA;

    } else {
      cerr << "** Unknown matrix dimension. Cannot compute inverse " << endl;
    }	
  }



  Matrix3d VoomExpSymmMatrix(const Matrix3d &A) {
    SelfAdjointEigenSolver<Matrix3d> SAES;
    SAES.compute(A);
    Vector3d v0 = SAES.eigenvectors().col(0);
    Vector3d v1 = SAES.eigenvectors().col(1);
    Vector3d v2 = SAES.eigenvectors().col(2);

    return exp(SAES.eigenvalues()[0])*( v0*v0.transpose() ) +
           exp(SAES.eigenvalues()[1])*( v1*v1.transpose() ) +
           exp(SAES.eigenvalues()[2])*( v2*v2.transpose() );
  }

  int LeviCivitaSymbol3(int i, int j, int k) {
    if (i == j || j == k || k == i) return 0;
    else if ((i == 0 && j == 1 && k == 2) || (i == 1 && j == 2 && k == 0) || (i == 2 && j == 0 && k == 1)) return 1;
    else if ((i == 0 && j == 2 && k == 1) || (i == 1 && j == 0 && k == 2) || (i == 2 && j == 1 && k == 0)) return -1;
    else {
      // For debugging
      cout << "Something went wrong with LeviCivitaSymbol3 Function" << endl;
      exit(EXIT_FAILURE);
    }
  }

};
