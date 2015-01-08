#include "LineQuadrature.h"
#include "QuadQuadrature.h"

namespace voom {
  QuadQuadrature::QuadQuadrature(unsigned int order)
  {
    Vector2d tempPoint(0.0, 0.0);
    LineQuadrature quad(order);
    _quadWeights.resize(order*order, 0.);
    _quadPoints.resize(order*order, tempPoint );
    const vector<Real>& weight    = quad.getQuadWeights();
    const vector<VectorXd>& point = quad.getQuadPoints();
    for(uint i = 0; i < order; i++)
      for(uint j = 0; j < order; j++) {
	_quadWeights[i*order + j] = weight[i]*weight[j];
	_quadPoints[i*order + j](0)= point[i](0);
	_quadPoints[i*order + j](1)= point[j](0);
      }
    return;
  }

  bool QuadQuadrature::check(unsigned int d) const {
    // create a random polynomial of degree=d, i.e., 
    // \sum_{i+j<=d} a_{ij} s_1^i s_2^j
    srand(time(0));
    vector<double> a;
    for(int i=0; i<=d; i++) {
      for(int j=0; i+j<=d; j++) {
	a.push_back( 100.0*( double(rand())/RAND_MAX - 0.5) );
      }
    }
    
    // Integrate the polynomial exactly using the formula from Cook's text:
    //
    //   \int_A s_1^i s_2^j (1-s_1-s_2)^k dA = 2A\frac{i!j!m!}{(2+i+j+k)!}
    // 
    // or with k=0
    //
    //   \int_A s_1^i s_2^j dA = 2A\frac{i!j!}{(2+i+j)!}
    //
    // For the standard quadrilateral A=1.0
    
    double I_exact=0.0;
    for(int i=0, k=0; i<=d; i++) 
      for(int j=0; i+j<=d; j++, k++) {
	if ((i%2 != 0) || (j%2 !=0)) continue;
	I_exact += a[k]*(2./(i+1))*(2./(j+1));
      }
    

    // Integrate the polynomial P(s) numerically
    //  \sum_p P(s_p) w_p A
    double I_numerical=0.0;
    for(uint qp = 0; qp < this->_quadWeights.size(); qp++) {
      double s1 = _quadPoints[qp](0);
      double s2 = _quadPoints[qp](1);
      for(int i=0, k=0; i<=d; i++) {
	for(int j=0; i+j<=d; j++, k++) {
	  I_numerical += a[k]*pow(s1,i)*pow(s2,j)*(_quadWeights[qp]);
	}
      }
    }
    
    cout << "QuadQuadrature::check for polynomial of degree "<< d << endl
	 << "I_exact     = " << I_exact << endl
	 << "I_numerical = " << I_numerical << endl
	 << "Error       = " << abs(I_exact-I_numerical)/abs(I_exact) 
	 << endl;

    double tol = 1.0e-8;
    if ( abs(I_numerical-I_exact) <= tol * abs(I_exact) ) {  
      cout << "QuadQuadrature::check( "<< d << ") PASSED!" << endl;
      return true;
    }
    cout << "QuadQuadrature::check( " << d <<  " ) FAILED!"
	      << endl;
    return false;
  }
} // namespace voom
