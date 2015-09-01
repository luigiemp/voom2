
#include "TriQuadrature.h"

namespace voom {
  
  TriQuadrature::TriQuadrature(unsigned int order)
  {
    Vector2d tempPoint(0.0, 0.0);

    switch(order) {  
    case 1 : // 1 point - degree of precision 1
	_quadWeights.push_back(0.5);
	_quadPoints.push_back(tempPoint);
	_quadPoints[0](0) = _quadPoints[0](1) = 1.0/3.0;
	break;

    case 2 : // 3 points - degree of precision 2
      // _quadWeights.assign(3, 1.0/6.0);
      // _quadPoints.resize(3, tempPoint);
      // _quadPoints[0](0) = 2.0/3.0; _quadPoints[0](1) = 1.0/6.0; 
      // _quadPoints[1](0) = 1.0/6.0; _quadPoints[1](1) = 2.0/3.0; 
      // _quadPoints[2](0) = 1.0/6.0; _quadPoints[2](1) = 1.0/6.0; 
      _quadWeights.assign(3, 0.5*1.0/3.0);
      _quadPoints.resize(3, tempPoint);
      _quadPoints[0](0) = 0.5; _quadPoints[0](1) = 0.5;
      _quadPoints[1](0) = 0.0; _quadPoints[1](1) = 0.5; 
      _quadPoints[2](0) = 0.5; _quadPoints[2](1) = 0.0;
      
      break;

    case 3 : // 6 points - degree of precision 3
	_quadWeights.assign(6, 1/12.0);
	_quadPoints.resize(6, tempPoint);
	_quadPoints[0](0) = 0.659027622374092;     _quadPoints[0](1) = 0.231933368553031;
	_quadPoints[1](0) = 0.659027622374092;     _quadPoints[1](1) = 0.109039009072877;
	_quadPoints[2](0) = 0.231933368553031;     _quadPoints[2](1) = 0.659027622374092;
	_quadPoints[3](0) = 0.231933368553031;     _quadPoints[3](1) = 0.109039009072877;
	_quadPoints[4](0) = 0.109039009072877;     _quadPoints[4](1) = 0.659027622374092;
	_quadPoints[5](0) = 0.109039009072877;     _quadPoints[5](1) = 0.231933368553031;

	break;
    
    default :
      std::cout << "TriQuadrature::::TriQuadrature(unsigned int order): No quadrature rule is implmented for order " 
		<< order << ". Terminating program execution." << std::endl;
      
    }; // switch loop

  }; // TriQuadrature::TriQuadrature(unsigned int order) 



  bool TriQuadrature::check(unsigned int degree) const
  {
    // Adapted from voom
    // Create a random polynomial of requested degree in x,y
    // create a random polynomial of degree=d, i.e., 
    // \sum_{i+j<=d} a_{ij} s_1^i s_2^j
    srand(time(NULL));
    std::vector<double> a;
    for(int i=0; i<=degree; i++) {
      for(int j=0; i+j<=degree; j++) {
        a.push_back( 100.0*(static_cast<double>(rand())/RAND_MAX - 0.33) );
      }
    }

    double Iexact = 0.0;
    for(int i=0, k=0; i<=degree; i++) {
      for(int j=0; i+j<=degree; j++, k++) {
        Iexact += a[k]*factorial(i)*factorial(j)/factorial(2+i+j);
      }
    }

    // Integrate the polynomial P(s) numerically
    //  \sum_p P(s_p) w_p A
    double Inumerical = 0.0, s1 = 0.0, s2 = 0.0;
    for (int p = 0; p < _quadWeights.size(); p++) {
      s1 = _quadPoints[p](0);
      s2 = _quadPoints[p](1);
      for(int i=0, k=0; i<=degree; i++) {
        for(int j=0; i+j<=degree; j++, k++) {
          Inumerical += a[k]*pow(s1,i)*pow(s2,j)*(_quadWeights[p]);
        }
      }
    }



    // Print test results to terminal
    double Error = abs(Iexact-Inumerical)/abs(Iexact);
    cout << "TriQuadrature::check for polynomial of degree " << degree << endl
	 << "I_exact     = " << Iexact     << endl
	 << "I_numerical = " << Inumerical << endl
	 << "Error       = " << Error      << endl;

    // Check if test has been passed
    double tol = 1.0e-7;
    if ( Error < tol*abs(Iexact) ) {
      cout << "TriQuadrature::check( " << degree << " ) PASSED!" << endl << endl;
      return true;
    }

    cout << "TriQuadrature::check( " << degree << " ) FAILED!" << endl << endl;

    return false;
  }; // TriQuadrature::check(unsigned int degree)

} // namespace voom
