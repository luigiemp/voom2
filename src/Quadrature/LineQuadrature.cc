#include "LineQuadrature.h"

namespace voom {
  LineQuadrature::LineQuadrature(unsigned int order)
  {
    VectorXd tempPoint(1);
    tempPoint(0) = 0.0;

    switch(order) {
      case 1: // 1 point
	_quadWeights.push_back(2.0);
	_quadPoints.push_back(tempPoint); 
	
	break;

      case 2: // 2 points
	_quadWeights.assign(2, 1.0);
	_quadPoints.resize(2, tempPoint);
	_quadPoints[0](0) = -0.577350269189626;				
	_quadPoints[1](0) =  0.577350269189626;
	
	break;

      case 3: // 3 point
	_quadWeights.resize(3);
	_quadWeights[0] = _quadWeights[2] = 5./9.;
	_quadWeights[1] = 8.0/9.0;

	_quadPoints.resize(3, tempPoint); 
	_quadPoints[0](0) = -0.774596669241483;
	_quadPoints[1](0) =  0.0;
	_quadPoints[2](0) =  0.774596669241483;
 
	break;

      case 4: // 4 points
	_quadWeights.resize(4);
	_quadWeights[0] = _quadWeights[1] = 0.652145154862546;
	_quadWeights[2] = _quadWeights[3] = 0.347854845137454;

	_quadPoints.resize(4, tempPoint); 
	_quadPoints[0](0) =  0.339981043584856;
	_quadPoints[1](0) = -0.339981043584856;
	_quadPoints[2](0) =  0.861136311594053;
	_quadPoints[3](0) = -0.861136311594053;
	
	break;

      case 5:  // 5points
	_quadWeights.resize(5);
	_quadWeights[0] = 0.568888888888889;
	_quadWeights[1] = _quadWeights[2] = 0.478628670499366;
	_quadWeights[3] = _quadWeights[4] = 0.236926885056189;
	
	_quadPoints.resize(5, tempPoint); 
	_quadPoints[0](0) = 0.;
	_quadPoints[1](0) =  0.538469310105683;
	_quadPoints[2](0) = -0.538469310105683;
	_quadPoints[3](0) =  0.906179854938664;
	_quadPoints[4](0) = -0.906179854938664;
	
	break;

      case 6: // 6 points
	_quadWeights.resize(6);
	_quadWeights[0] = _quadWeights[1] = 0.467913934572691;
	_quadWeights[2] = _quadWeights[3] = 0.360761573048139;
	_quadWeights[4] = _quadWeights[5] = 0.171324492379170;

	_quadPoints.resize(6, tempPoint); 
	_quadPoints[0](0) =  0.238619186083197;
	_quadPoints[1](0) = -0.238619186083197;
	_quadPoints[2](0) =  0.661209386466265;
	_quadPoints[3](0) = -0.661209386466265;
	_quadPoints[4](0) =  0.932469514203152;
	_quadPoints[5](0) = -0.932469514203152;

	break;

      case 7: // 7 points	
	_quadWeights.resize(7);
	_quadWeights[0] = 0.417959183673469;
	_quadWeights[1] = _quadWeights[2] = 0.381830050505119;
	_quadWeights[3] = _quadWeights[4] = 0.279705391489277;
	_quadWeights[5] = _quadWeights[6] = 0.129484966168870;

	_quadPoints.resize(7, tempPoint); 
	_quadPoints[0](0) =  0.;
	_quadPoints[1](0) =  0.405845151377397;
	_quadPoints[2](0) = -0.405845151377397;
	_quadPoints[3](0) =  0.741531185599394;
	_quadPoints[4](0) = -0.741531185599394;
	_quadPoints[5](0) =  0.949107912342759;
	_quadPoints[6](0) = -0.949107912342759;

	break;

      case 8: // 8 points
	_quadWeights.resize(8);
	_quadWeights[0] = _quadWeights[1] = 0.362683783378362;
	_quadWeights[2] = _quadWeights[3] = 0.313706645877887;
	_quadWeights[4] = _quadWeights[5] = 0.222381034453374;
	_quadWeights[6] = _quadWeights[7] = 0.101228536290376;

	_quadPoints.resize(8, tempPoint); 
	_quadPoints[0](0) =  0.183434642495690;
	_quadPoints[1](0) = -0.183434642495690;
	_quadPoints[2](0) =  0.525532409916329;
	_quadPoints[3](0) = -0.525532409916329;
	_quadPoints[4](0) =  0.796666477413627;
	_quadPoints[5](0) = -0.796666477413627;
	_quadPoints[6](0) =  0.960289856497536;
	_quadPoints[7](0) = -0.960289856497536;

	break;
	
      case 9: // 9 points
	_quadWeights.resize(9);
	_quadWeights[0] = 0.330239355001260;
	_quadWeights[1] = _quadWeights[2] = 0.180648160694857;
	_quadWeights[3] = _quadWeights[4] = 0.081274388361574;
	_quadWeights[5] = _quadWeights[6] = 0.312347077040003;
	_quadWeights[7] = _quadWeights[8] = 0.260610696402935;

	_quadPoints.resize(9, tempPoint); 
	_quadPoints[0](0) = 0.;
	_quadPoints[1](0) =  0.836031107326636;
	_quadPoints[2](0) = -0.836031107326636;
	_quadPoints[3](0) =  0.968160239507626;
	_quadPoints[4](0) = -0.968160239507626;
	_quadPoints[5](0) =  0.3242534234038089;
	_quadPoints[6](0) = -0.3242534234038089;
	_quadPoints[7](0) =  0.613371432700590;
	_quadPoints[8](0) = -0.613371432700590;
	
	break;

      case 10: // 10 points
	_quadWeights.resize(10);
	_quadWeights[0] = _quadWeights[1] = 0.295524224714753;
	_quadWeights[2] = _quadWeights[3] = 0.269266719309996;
	_quadWeights[4] = _quadWeights[5] = 0.219086362515982;
	_quadWeights[6] = _quadWeights[7] = 0.149451349150581;
	_quadWeights[8] = _quadWeights[9] = 0.066671344308688 ;

	_quadPoints.resize(10, tempPoint); 
	_quadPoints[0](0) =  0.148874338981631;
	_quadPoints[1](0) = -0.148874338981631;
	_quadPoints[2](0) =  0.433395394129247;
	_quadPoints[3](0) = -0.433395394129247;
	_quadPoints[4](0) =  0.679409568299024;
	_quadPoints[5](0) = -0.679409568299024;
	_quadPoints[6](0) =  0.865063366688985;
	_quadPoints[7](0) = -0.865063366688985;
	_quadPoints[8](0) =  0.973906528517172;
	_quadPoints[9](0) = -0.973906528517172;
       
	break;

      case 11: // 11 points
	_quadWeights.resize(11);
	_quadWeights[0] = 0.272925086777901;
	_quadWeights[1] = _quadWeights[2] = 0.262804544510247;
	_quadWeights[3] = _quadWeights[4] = 0.233193764591990;
	_quadWeights[5] = _quadWeights[6] = 0.186290210927734;
	_quadWeights[7] = _quadWeights[8] = 0.125580369464905;
	_quadWeights[9] = _quadWeights[10] = 0.055668567116174; 

	_quadPoints.resize(11, tempPoint);
	_quadPoints[0](0) =  0.;
	_quadPoints[1](0) =  0.269543155952345;
	_quadPoints[2](0) = -0.269543155952345;
	_quadPoints[3](0) =  0.519096129110681;
	_quadPoints[4](0) = -0.519096129110681;
	_quadPoints[5](0) =  0.730152005574049;
	_quadPoints[6](0) = -0.730152005574049;
	_quadPoints[7](0) =  0.887062599768095;
	_quadPoints[8](0) = -0.887062599768095;
	_quadPoints[9](0) =  0.978228658146057;
	_quadPoints[10](0) = -0.978228658146057;

	break;

      default:
	std::cout << "LineQuadrature::_initialize(): No quadrature rule is implmented for order " 
		<< order << ".  Quadrature uninitialized." << std::endl;
    } // Switch loop

  } // Line quadrature constructor



  bool LineQuadrature::check(unsigned int degree) const {

    // Create a random polynomial of requested degree in x,y,z
    srand(time(NULL));
    vector<Real > coeffx;
    int i = 0, q = 0;
    for(i = 0; i <= degree; i++) {
      coeffx.push_back( 100.0*( (double(rand())/RAND_MAX) - 0.5) );
    }

    // Compute the polynomial analytically in the [-1,1] reference domain
    double Iexact = 0.0, Ie_x = 0.0;
    for(i = 0; i <= degree; i++) {
      Ie_x += (1.0+pow(-1.0,i))*coeffx[i]/(i+1);
    }
    Iexact = Ie_x;

    // Compute the polynomial numerically in the [-1,1] reference domain
    double Inumerical = 0.0, In_x = 0.0, weight = 0.0 ,x = 0.0;
    for (q = 0; q < _quadWeights.size(); q++) {
      weight = _quadWeights[q];
      x = _quadPoints[q](0);

      In_x = 0.0;
      for(i = 0; i <= degree; i++) {
	In_x += coeffx[i]*pow(x,i);
      }
      Inumerical += In_x*weight;
    }
      
    // Print test results to terminal
    double Error = abs(Iexact-Inumerical)/abs(Iexact);
    cout << "LineQuadrature::check for polynomial of degree " << degree << endl
	 << "I_exact     = " << Iexact     << endl
	 << "I_numerical = " << Inumerical << endl
	 << "Error       = " << Error      << endl;

    // Check if test has been passed
    double tol = 1.0e-7;
    if ( Error < tol*abs(Iexact) ) {
      cout << "LineQuadrature::check( " << degree << " ) PASSED!" << endl << endl;
      return true;
    }

    cout << "LineQuadrature::check( " << degree << " ) FAILED!" << endl << endl;

    return false;
  } // CheckConsistency


} // namespace voom
