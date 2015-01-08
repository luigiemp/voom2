
#include "HexQuadrature.h"

namespace voom {
  
  HexQuadrature::HexQuadrature(unsigned int order)
  {
    Vector3d tempPoint(0.0, 0.0, 0.0);

    switch(order) {  
      case 1 : // 1 point 
	_quadWeights.push_back(8.0);
	_quadPoints.push_back(tempPoint);
	break;

      case 2 : // 8 points
	_quadWeights.assign(8, 1.0);
	_quadPoints.resize(8, tempPoint);
	_quadPoints[0](0)=-0.577350269189626; _quadPoints[0](1)=-0.577350269189626; _quadPoints[0](2)=-0.577350269189626;
	_quadPoints[1](0)= 0.577350269189626; _quadPoints[1](1)=-0.577350269189626; _quadPoints[1](2)=-0.577350269189626;
	_quadPoints[2](0)=-0.577350269189626; _quadPoints[2](1)= 0.577350269189626; _quadPoints[2](2)=-0.577350269189626;
	_quadPoints[3](0)= 0.577350269189626; _quadPoints[3](1)= 0.577350269189626; _quadPoints[3](2)=-0.577350269189626;
	_quadPoints[4](0)=-0.577350269189626; _quadPoints[4](1)=-0.577350269189626; _quadPoints[4](2)= 0.577350269189626;
	_quadPoints[5](0)= 0.577350269189626; _quadPoints[5](1)=-0.577350269189626; _quadPoints[5](2)= 0.577350269189626;
	_quadPoints[6](0)=-0.577350269189626; _quadPoints[6](1)= 0.577350269189626; _quadPoints[6](2)= 0.577350269189626;
	_quadPoints[7](0)= 0.577350269189626; _quadPoints[7](1)= 0.577350269189626; _quadPoints[7](2)= 0.577350269189626;
	break;

    case 3: {
      _quadPoints.resize(27, tempPoint);
      _quadWeights.resize(27);
      Real qp[] = { sqrt(3./5.), 0., -sqrt(3./5.) };
      Real w[] = { 5./9., 8./9., 5./9.};
      for(uint i = 0; i < 3; i++)
	for(uint j = 0; j < 3; j++)
	  for(uint k = 0; k < 3; k++) {
	    _quadPoints[ (i*3 + j)*3 + k ](0) = qp[i];
	    _quadPoints[ (i*3 + j)*3 + k ](1) = qp[j];
	    _quadPoints[ (i*3 + j)*3 + k ](2) = qp[k];
	    _quadWeights[ 9*i + 3*j +k] = w[i]*w[j]*w[k];
	  }
      break;
    }
    default :
      std::cout << "HexQuadrature::::HexQuadrature(unsigned int order): No quadrature rule is implmented for order " 
		<< order << ". Terminating program execution." << std::endl;
      

    }; // switch loop

  }; // HexQuadrature::HexQuadrature(unsigned int order) 



  bool HexQuadrature::check(unsigned int degree) const
  {
    // Create a random polynomial of requested degree in x,y,z
    srand(time(NULL));
    vector<Real > coeffx, coeffy, coeffz;
    int i = 0, q = 0;
    for(i = 0; i <= degree; i++) {
      coeffx.push_back( 100.0*( (double(rand())/RAND_MAX) - 0.5) );
      coeffy.push_back( 100.0*( (double(rand())/RAND_MAX) - 0.5) );
      coeffz.push_back( 100.0*( (double(rand())/RAND_MAX) - 0.5) );
    }

    // Compute the polynomial analitically in the [-1,1] reference domain
    double Iexact = 0.0, Ie_x = 0.0, Ie_y = 0.0, Ie_z = 0.0;
    for(i = 0; i <= degree; i++) {
      Ie_x += (1.0+pow(-1.0,i))*coeffx[i]/(i+1);
      Ie_y += (1.0+pow(-1.0,i))*coeffy[i]/(i+1);
      Ie_z += (1.0+pow(-1.0,i))*coeffz[i]/(i+1);
    }
    Iexact = Ie_x*Ie_y*Ie_z;

    // Compute the polynomial numerically in the [-1,1] reference domain
    double Inumerical = 0.0, In_x = 0.0, In_y = 0.0, In_z = 0.0, weight = 0.0 ,x = 0.0, y = 0.0, z = 0.0;
    for (q = 0; q < _quadWeights.size(); q++) {
      weight = _quadWeights[q];
      x = _quadPoints[q](0);
      y = _quadPoints[q](1);
      z = _quadPoints[q](2);

      In_x = 0.0; In_y = 0.0; In_z = 0.0;
      for(i = 0; i <= degree; i++) {
	In_x += coeffx[i]*pow(x,i);
	In_y += coeffy[i]*pow(y,i);
	In_z += coeffz[i]*pow(z,i);
      }
      Inumerical += In_x*In_y*In_z*weight;
    }
      
    // Print test results to terminal
    double Error = abs(Iexact-Inumerical)/abs(Iexact);
    cout << "HexQuadrature::check for polynomial of degree " << degree << endl
	 << "I_exact     = " << Iexact     << endl
	 << "I_numerical = " << Inumerical << endl
	 << "Error       = " << Error      << endl;

    // Check if test has been passed
    double tol = 1.0e-7;
    if ( Error < tol*abs(Iexact) ) {
      cout << "HexQuadrature::check( " << degree << " ) PASSED!" << endl << endl;
      return true;
    }

    cout << "HexQuadrature::check( " << degree << " ) FAILED!" << endl << endl;

    return false;
  }; // HexQuadrature::check(unsigned int degree)

} // namespace voom


      
      
 
