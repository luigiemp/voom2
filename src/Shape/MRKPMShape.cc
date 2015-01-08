//-*-C++-*-
#include "MRKPMShape.h"
#include "LineQuadrature.h"
#include "QuadQuadrature.h"

namespace voom{
  // Constructor
  MRKPMShape::MRKPMShape(const vector<VectorXd> & Nodes, const VectorXd & Point,
			 const Real support, const Real radius, 
			 const Real suppHat):
    MFShape(Nodes), _support(support), _radius(radius), _supportHat(suppHat) 
  {
    update(Point);
  }

  // Compute Kernel function
  Real MRKPMShape::computePhi(Real z, const kernelType kernel) {
    switch (kernel) {
    case CONSTANT:
      if ( z <= 1.) return 1;
      break;
    case LINEAR:
      if (z <= 1.) return (1 - z);
      break;
    case CUBIC:
      if ( z > 0.5 && z <= 1) return (4./3. - 4.*z + 4.*z*z - 4./3.*z*z*z);
      if ( z <= 0.5) return (2./3. - 4.*z*z + 4.*z*z*z);
      break;
    case QUARTIC:
      if (z <=1 ) return (1. - 6.*z*z + 8.*z*z*z - 3.*z*z*z*z);
      break;
    case EXPONENTIAL:
      if (z <=1 ) return exp(-z*z/0.09);
      break;
    default:
	cerr << "** ERROR: Unknown kernel function provided" << endl;
	cerr << "** Exiting...\n";
	exit(1);
    } // switch statement
    return 0.;
  }

  // Compute shape functions and derivatives at the given Point
  void MRKPMShape::update(const VectorXd & Point)
  {
    uint dim = _nodes[0].size();
    assert(dim == Point.size());
    // Reset derivatives value to 0
    for(uint m = 0; m < _shapeNum; m++) _DN[m] = VectorXd::Zero(dim);
    // Base class quadrature pointer
    Quadrature* quad;
    if ( dim == 2) quad = new LineQuadrature(10);
    else quad = new QuadQuadrature(10);
    const vector<VectorXd>& quadPoints = quad->getQuadPoints();
    const vector<Real>&     quadWeight = quad->getQuadWeights();
    const int numQP                    = quadPoints.size();
    // Initialize points of the contour (circle or sphere)
    vector< vector<Real> > functionValues;

    for(uint i = 0; i < numQP; i++) {
      // Set up the Point position
      VectorXd offset = Point;
      if (dim == 2) {
	offset(0) += _radius*cos(M_PI*quadPoints[i](0));
	offset(1) += _radius*sin(M_PI*quadPoints[i](0));
      } else {
	Real phi = M_PI*(quadPoints[i](0) + 1.); // phi 0 to 2 pi
	Real theta = M_PI/2.*(quadPoints[i](1) + 1.); // theta 0 to pi
	// dA = r^2 sin(theta) dTheta dPhi
	offset(0) += _radius*sin(theta)*cos(phi);
	offset(1) += _radius*sin(theta)*sin(phi);
	offset(2) += _radius*cos(theta);
      }
      this->computeFunction( offset );
      functionValues.push_back( _N );
    } // Loop over numQP

    // Compute Derivatives Now
    const Real scalar  = (dim == 2) ? 1./_radius : 3.*M_PI/(8.*_radius);
    vector<Real> factor(dim, 0.);
    for(uint i = 0; i < numQP; i++) {
      if (dim == 2) {
	factor[0] = cos(M_PI*quadPoints[i](0)) * quadWeight[i];
	factor[1] = sin(M_PI*quadPoints[i](0)) * quadWeight[i];
      } else {
	const Real phi = M_PI*(quadPoints[i](0) + 1.);
	const Real theta = M_PI/2*(quadPoints[i](1) + 1.);
	factor[0] = pow(sin(theta),2)*cos(phi) * quadWeight[i];
        factor[1] = pow(sin(theta),2)*sin(phi) * quadWeight[i];
        factor[2] = sin(theta)*cos(theta) * quadWeight[i];
      }
      for(uint j = 0; j < _nodes.size(); j++) 
	for(uint k = 0; k < dim ; k++) 
	  _DN[j](k) += functionValues[i][j]*factor[k]*scalar;

    } // Loop over numQp

    delete quad;

    // Get Shape Functions at required point
    computeFunction(Point);
  }
  
  // Compute Shape Functions
  void MRKPMShape::computeFunction(const VectorXd& Point) {
    const int dim = Point.size();
    MatrixXd M = MatrixXd::Zero(dim + 1, dim + 1); // Moment Matrix
    VectorXd Fhat = VectorXd::Zero(dim + 1);
    
    // Computing Moment Matrix
    for(uint a = 0; a < _shapeNum; a++){
      VectorXd dX = Point - _nodes[a];
      const Real distance = dX.norm();
      const Real phi = computePhi(distance/_support);
      const Real phiHat = computePhi(distance/_supportHat);
      VectorXd Ht(dim + 1); Ht(0) = 1.;
      for(uint j = 0; j < dim; j++) Ht(j+1) = dX[j];
      M = M + phi * Ht * Ht.transpose();
      Fhat = Fhat + phiHat * Ht;
    } // Loop over number of nodes in support
    MatrixXd Minv = MatrixXd::Zero(dim + 1, dim + 1);
    inv(M, Minv);
    VectorXd H0 = VectorXd::Zero(dim + 1); H0[0] = 1.;
    Fhat = H0 - Fhat;
    // Loop over all nodes in support again
    for(uint a = 0; a < _shapeNum; a++) {
      VectorXd dX = Point - _nodes[a];
      const Real distance = dX.norm();
      const Real phi = computePhi(distance/_support);
      const Real phiHat = computePhi(distance/_supportHat);
      VectorXd H(dim + 1); H[0] = 1.;
      for(uint j = 0; j < dim; j++) H[j+1] = dX[j];
      _N[a] = (Minv*Fhat).dot(H)*phi + phiHat;
    }
  }

} // namespace voom

