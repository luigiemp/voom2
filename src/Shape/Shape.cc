#include "Shape.h"

namespace voom 
{
  
  bool Shape::checkConsistency(VectorXd Point, const Real eps, const Real tol)
  {
    // Test dimension (1D, 2D, 3D)
    const uint dim = Point.size();
    
    const uint NumF = this->getShapeFunctionNum();
    vector<Real > N(NumF, 0.0);
    vector<VectorXd > DNnumerical(NumF, VectorXd::Zero(dim));
   
    Real error = 0.0, norm = 0.0;

    // Compute numerical derivatives
    for(uint i = 0; i < dim; i++)
    {
      Point(i) += eps;
      this->update(Point);
      for(uint a = 0; a < NumF; a++) 
	DNnumerical[a](i) = this->getN(a);
      
      Point(i) -= 2*eps;
      this->update(Point);
      for(uint a = 0; a < NumF; a++)
	DNnumerical[a](i)-= this->getN(a);
      
      Point(i) += eps;
      for(uint a = 0; a < NumF; a++) 
	DNnumerical[a](i) /= (2.0*eps);    
    }

    // Compute error norm
    for(uint a = 0; a < NumF; a++) {
      for(uint i = 0; i < dim; i++) {
	error += pow(this->getDN(a,i) - DNnumerical[a](i), 2);
	norm  += pow(this->getDN(a,i), 2);
      }
    }
    norm = sqrt(norm);
    error = sqrt(error);
    
    cout << "Error = " << error << " Norm = " << norm << endl;
    if ( abs(error) < norm * tol) {
      cout << "Shape consistency check passed" << endl;
      return true;
    }
    else {
      cout << "Shape consistency check failed" << endl
	   << setw(10) << "n"
	   << setw(24) << "analytical"
	   << setw(24) << "numerical" << endl;   
      for(uint a = 0; a < NumF; a++)
	for(uint i = 0; i < dim; i++)
	  cout << setw(10) << a << setw(24) << this->getDN(a,i)
	       << setw(24) << DNnumerical[a](i) << endl;
    }

    return false;
  } // consistency check



  bool Shape::checkPartitionUnity(const VectorXd Point, const Real tol) 
  {
    // Number of shape functions
    const uint NumF = this->getShapeFunctionNum();
    Real Unity = 0.0;
    VectorXd Zero;
    Zero = VectorXd::Zero( Point.size() );
   
    // Compute sum of shape functions at Point
    this->update(Point);
    for(uint a = 0; a < NumF; a++)
      Unity += this->getN(a); 
    
    // Sum of shape function derivatives at a Point
    for(uint a = 0; a < NumF; a++) 
      for(uint b = 0; b < Point.size(); b++)
	Zero[b] += this->getDN(a,b);

    if ( abs(Unity - 1.0) < tol && Zero.norm() < tol) {
      cout << "Shape partition of unity check passed" << endl << endl;
      return true;
    }
    else {
      cout << "Shape partition of unity check failed." << endl
	   << "Shape function sum = " << Unity << endl
	   << "Shape Derivative sum = ";
      for(uint m = 0; m < Point.size(); m++) cout << Zero(m) << " " ;
      cout << endl;
    }

    cout << endl;
    return false;

  } // partition of unity check

} // namespace voom
