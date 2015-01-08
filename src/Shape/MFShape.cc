#include "MFShape.h"

namespace voom {
  //Constructor
  MFShape::MFShape(const vector<VectorXd>& nodes): _nodes(nodes) {
    _shapeNum = _nodes.size();
    _N.resize(_shapeNum, 0.0);
    _DN.resize(_shapeNum, VectorXd::Zero(_nodes[0].size()) );
  }

  // Overloaded partition of unity
  bool MFShape::checkReproducingCondition(VectorXd Point, const Real tol)
  {  
    // Compute Shape functions at the point
    this->update(Point);

    // Checking reproducing conditions
    VectorXd pos = VectorXd::Zero(Point.size());
    for(uint a = 0; a < _nodes.size(); a++)
      pos += _nodes[a] * this->getN(a);

    VectorXd Unity = VectorXd::Zero(Point.size()),
      One = VectorXd::Ones(Point.size());
    for(uint a = 0; a < _nodes.size(); a++)
      for(uint b = 0; b < Point.size(); b++)
	Unity[b] += this->getDN(a,b)*_nodes[a][b];

    if( (pos - Point).norm() < tol && (Unity - One).norm() < tol) {
      cout << "Shape reproducing condition check passed" << endl << endl;
      return true;
    } 
    else {
      cout << "Shape reproducing condition check failed" << endl;
      cout << "Sum N_iX_i\tX(s)\n";
      for(uint m = 0; m < pos.size(); m++) 
	cout << pos(m) << "\t" << Point(m) << endl;
      cout << "Error: " << (pos - Point).norm() << endl;
    
      cout << endl << "Sum X_i DN_i = ";
      for(uint m = 0; m < pos.size(); m++)
	cout << Unity(m) << " " ;
      cout << endl << "Error: " << (Unity - One).norm() << endl;
      return false;
    }
  } // End of function
}
