//-*-C++-*-
#include "LMEShape.h"

namespace voom{

  LMEShape::LMEShape(const vector<VectorXd > & Nodes, const VectorXd & Point, 
		     const Real beta, const Real tol, const uint MaxIter):
    MFShape(Nodes), _beta(beta), _tol(tol), _maxIter(MaxIter)
  {
    update(Point);
  }


  void LMEShape::update(const VectorXd & Point)
  {
    uint dim = _nodes[0].size();
    assert(dim==Point.size());

    // Compute LME shape functions and shape functions derivatives
    VectorXd lambda = VectorXd::Zero(dim);
    VectorXd r = VectorXd::Zero(dim), dN = VectorXd::Zero(dim);
    MatrixXd J = MatrixXd::Zero(dim,dim), Jinv = MatrixXd::Zero(dim,dim);

    // Compute error
    r = LME_r(Point, lambda);    
    Real res = pow(r.dot(r), 0.5);
   
    // Start minimization problem to determine LME functions
    uint iter = 0;
    while (res>_tol && iter<_maxIter)
    {
      J = LME_J(Point, lambda);
      Jinv = J.inverse();
      lambda -= Jinv*r;

      r = LME_r(Point, lambda);
      res = pow(r.dot(r), 0.5);
      iter++;
    }

    for(uint i=0; i<_nodes.size(); i++)
    {
      _N[i]  = LME_pa(Point, lambda, i);
      _DN[i] = _N[i]*Jinv*(_nodes[i]-Point);
    }
  }



  Real LMEShape::LME_pa(const VectorXd & Point,
			const VectorXd & lambda,
			uint j)
  {
    Real Z = 0.0;
    VectorXd temp(_nodes[0].size());

    for(int i=0; i<_nodes.size(); i++)
    { 
      temp = Point - _nodes[i];
      Z += exp( - _beta*temp.dot(temp) + temp.dot(lambda)  );
    }
    
    temp = Point - _nodes[j];
    return exp( - _beta*temp.dot(temp) + temp.dot(lambda) )/Z;
  }



  VectorXd LMEShape::LME_r(const VectorXd & Point, 
			   const VectorXd & lambda)
  {
    VectorXd r = VectorXd::Zero(Point.size());
    
    for(int i=0; i<_nodes.size(); i++){
      r += LME_pa(Point, lambda, i) * (Point-_nodes[i]);
    }

    return r;
  }



  MatrixXd LMEShape::LME_J(const VectorXd & Point,
			   const VectorXd & lambda)
  {
    uint dim = Point.size();
    MatrixXd J = MatrixXd::Zero(dim,dim);
    VectorXd temp = VectorXd::Zero(dim);

    for(int i=0; i<_nodes.size(); i++){  
      temp = Point - _nodes[i];
      J += LME_pa(Point, lambda, i) * (temp*temp.transpose() );
    }
    
    temp = LME_r(Point, lambda);

    J -= temp*temp.transpose();
    return J;
  }



} // namespace voom

