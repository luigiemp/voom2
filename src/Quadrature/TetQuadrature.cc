#include "TetQuadrature.h"

namespace voom {
  
  TetQuadrature::TetQuadrature(unsigned int order)
  {
    Vector3d tempPoint(0.0, 0.0, 0.0);

    switch(order) {  
    case 1 : // 1 point - degree of precision 1
	_quadWeights.push_back(1.0/6.0);
	_quadPoints.push_back(tempPoint);
	_quadPoints[0](0) = 0.25; _quadPoints[0](1) = 0.25; _quadPoints[0](2) = 0.25;
	break;

    case 2 : // 4 points - degree of precision 2
	_quadWeights.assign(4, 1.0/24.0);
	_quadPoints.resize(4, tempPoint);
	_quadPoints[0](0) = 0.5854101966249685; _quadPoints[0](1) = 0.1381966011250105; _quadPoints[0](2) = 0.1381966011250105;
	_quadPoints[1](0) = 0.1381966011250105; _quadPoints[1](1) = 0.5854101966249685; _quadPoints[1](2) = 0.1381966011250105;
	_quadPoints[2](0) = 0.1381966011250105; _quadPoints[2](1) = 0.1381966011250105; _quadPoints[2](2) = 0.5854101966249685;
	_quadPoints[3](0) = 0.1381966011250105; _quadPoints[3](1) = 0.1381966011250105; _quadPoints[3](2) = 0.1381966011250105;
	break;

    case 3 : // 5 points - degree of precision 3
	_quadWeights.assign(5, 0.45/6.0);
	_quadWeights[0] = -0.8/6.0;
	_quadPoints.resize(5, tempPoint);
	_quadPoints[0](0) = 0.25;     _quadPoints[0](1) = 0.25;      _quadPoints[0](2) = 0.25;
	_quadPoints[1](0) = 0.5;      _quadPoints[1](1) = 1.0/6.0;   _quadPoints[1](2) = 1.0/6.0;
	_quadPoints[2](0) = 1.0/6.0;  _quadPoints[2](1) = 0.5;       _quadPoints[2](2) = 1.0/6.0;
	_quadPoints[3](0) = 1.0/6.0;  _quadPoints[3](1) = 1.0/6.0;   _quadPoints[3](2) = 0.5;
	_quadPoints[4](0) = 1.0/6.0;  _quadPoints[4](1) = 1.0/6.0;   _quadPoints[4](2) = 1.0/6.0;

	break;

    // case 4 : // 5 points - degree of precision 2
    // 	_quadWeights.assign(5, 0.45/6.0);
    // 	_quadWeights[0] = -0.8/6.0;
    // 	_quadPoints.resize(5, tempPoint);
    // 	_quadPoints[0](0) = 0.25;     _quadPoints[0](1) = 0.25;      _quadPoints[0](2) = 0.25;
    // 	_quadPoints[1](0) = 0.5;      _quadPoints[1](1) = 1.0/6.0;   _quadPoints[1](2) = 1.0/6.0;
    // 	_quadPoints[2](0) = 1.0/6.0;  _quadPoints[2](1) = 0.5;       _quadPoints[2](2) = 1.0/6.0;
    // 	_quadPoints[3](0) = 1.0/6.0;  _quadPoints[3](1) = 1.0/6.0;   _quadPoints[3](2) = 0.5;
    // 	_quadPoints[4](0) = 1.0/6.0;  _quadPoints[4](1) = 1.0/6.0;   _quadPoints[4](2) = 1.0/6.0;

    // 	break;


   
    // case 3   %  K4  N=11
    //   xa= [0.2500000000000000, 0.7857142857142857, 0.0714285714285714, 0.0714285714285714, 0.0714285714285714, ...
    // 	   0.1005964238332008, 0.3994035761667992, 0.3994035761667992, 0.3994035761667992, 0.1005964238332008, 0.1005964238332008];
    // ya= [0.2500000000000000, 0.0714285714285714, 0.0714285714285714, 0.0714285714285714, 0.7857142857142857, ...
    // 	 0.3994035761667992, 0.1005964238332008, 0.3994035761667992, 0.1005964238332008, 0.3994035761667992, 0.1005964238332008];
    // za= [0.2500000000000000, 0.0714285714285714, 0.0714285714285714, 0.7857142857142857, 0.0714285714285714, ...
    // 	 0.3994035761667992, 0.3994035761667992, 0.1005964238332008, 0.1005964238332008, 0.1005964238332008, 0.3994035761667992];
    // wt=[-0.0789333333333333, 0.0457333333333333, 0.0457333333333333, 0.0457333333333333, 0.0457333333333333, ...
    // 	0.1493333333333333, 0.1493333333333333, 0.1493333333333333, 0.1493333333333333, 0.1493333333333333, 0.1493333333333333]/6;
    
    // case 4  %K6  N=15
    //   xa=[0.2500000000000000, 0.0000000000000000, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333, ...
    // 	  0.7272727272727273, 0.0909090909090909, 0.0909090909090909, 0.0909090909090909, 0.4334498464263357, ...
    // 	  0.0665501535736643, 0.0665501535736643, 0.0665501535736643, 0.4334498464263357, 0.4334498464263357];
    // ya=[0.2500000000000000, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333, 0.0000000000000000, ...
    // 	0.0909090909090909, 0.0909090909090909, 0.0909090909090909, 0.7272727272727273, 0.0665501535736643, ...
    // 	0.4334498464263357, 0.0665501535736643, 0.4334498464263357, 0.0665501535736643, 0.4334498464263357];
    // za=[0.2500000000000000, 0.3333333333333333, 0.3333333333333333, 0.0000000000000000, 0.3333333333333333, ...
    // 	0.0909090909090909, 0.0909090909090909, 0.7272727272727273, 0.0909090909090909, 0.0665501535736643, ...
    // 	0.0665501535736643, 0.4334498464263357, 0.4334498464263357, 0.4334498464263357, 0.0665501535736643];
    // wt=[0.1817020685825351, 0.0361607142857143, 0.0361607142857143, 0.0361607142857143, 0.0361607142857143, ...
    // 	0.0698714945161738, 0.0698714945161738, 0.0698714945161738, 0.0698714945161738, 0.0656948493683187, ...
    // 	0.0656948493683187, 0.0656948493683187, 0.0656948493683187, 0.0656948493683187, 0.0656948493683187]/6;
    
    default :
      std::cout << "TetQuadrature::::TetQuadrature(unsigned int order): No quadrature rule is implmented for order " 
		<< order << ". Terminating program execution." << std::endl;
      
    }; // switch loop

  }; // TetQuadrature::TetQuadrature(unsigned int order) 



  bool TetQuadrature::check(unsigned int degree) const
  {
    // Adapted from voom
    // Create a random polynomial of requested degree in x,y,z
    // \sum_{i+j+k<=d} a_{ijk} s_1^i s_2^j s_3^k
    srand(time(NULL));
    std::vector<double> a;
    for(int i = 0, q = 0; i <= degree; i++) {
      for(int j = 0; i+j <= degree; j++) {
        for(int k = 0; i+j+k <= degree; k++, q++) {
          a.push_back( 100.0*(double(rand())/double(RAND_MAX) - 0.33) );
        }
      }
    }

    double Iexact = 0.0;
    for(int i = 0, q = 0; i <= degree; i++) {
      for(int j = 0; i+j <= degree; j++) {
        for(int k = 0; i+j+k <= degree; k++, q++) {
          Iexact += a[q]*factorial(i)*factorial(j)*factorial(k)/factorial(3+i+j+k);
        }
      }
    }
    
   
    // Integrate the polynomial P(s) numerically
    //  \sum_p P(s_p) w_p A
    double Inumerical = 0.0;
    double s1 = 0.0, s2 = 0.0, s3 = 0.0;
    for (int p = 0; p < _quadWeights.size(); p++) {
      s1 = _quadPoints[p](0);
      s2 = _quadPoints[p](1);
      s3 = _quadPoints[p](2);
      for(int i = 0, q = 0; i <= degree; i++) {
        for(int j = 0; i+j <= degree; j++) {
          for(int k = 0; i+j+k <= degree; k++, q++) {
            Inumerical += a[q]*pow(s1,i)*pow(s2,j)*pow(s3,k)*(_quadWeights[p]);
          }
        }
      }
    }

     
    // Print test results to terminal
    double Error = abs(Iexact-Inumerical)/abs(Iexact);
    cout << "TetQuadrature::check for polynomial of degree " << degree << endl
	 << "I_exact     = " << Iexact     << endl
	 << "I_numerical = " << Inumerical << endl
	 << "Error       = " << Error      << endl;

    // Check if test has been passed
    double tol = 1.0e-7;
    if ( Error < tol*abs(Iexact) ) {
      cout << "TetQuadrature::check( " << degree << " ) PASSED!" << endl << endl;
      return true;
    }

    cout << "TetQuadrature::check( " << degree << " ) FAILED!" << endl << endl;

    return false;
  }; // TetQuadrature::check(unsigned int degree)

} // namespace voom


      
      
 
