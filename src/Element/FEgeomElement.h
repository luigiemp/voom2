//-*-C++-*-
/*!
  \file FEgeomElement.h

  \brief Base class implementation for FE geometry element class. This class
  performs the operation of computing the shape function derivatives 
  w.r.t. global CSYS using the shape-function derivatives in isoparametric
  CSYS. There is no physics implemented in this class. 
*/

#if !defined(__FEgeomElement_h__)
#define __FEgeomElement_h__

#include "GeomElement.h"

namespace voom {

  class FEgeomElement: public GeomElement {    
  public:
    /*!
      Specialized constructors that computes shape functions and
      shape functions derivatives.
      We need \f$\frac{\partial N}{\partial x} \f$. Using the chain rule
      we can compute this as
      \f[
      \frac{\partial N}{\partial \xi} \frac{\partial \xi}{\partial x}
      \f]
      The second term is the inverse of the jacobian. The jacobian is 
      computed as
      \f[
      \frac{\partial x_a}{\partial \xi_b} = \sum_{m=1}^{3} 
      \frac{\partial N_m}{\partial \xi_j} x_m^a
      \f]
    */

    FEgeomElement(const int elemID, const vector<int > & nodesID, 
		  const vector<VectorXd > & nodesX, 
		  vector<Shape* > shape, Quadrature* quadrature);

    //! Get number of quadrature points
    uint getNumberOfQuadPoints() { return _QPweights.size(); }

    //! Get weight for quadrature point q
    Real getQPweights(uint q) {return _QPweights[q]; }

    //! Get shape functions values at quadrature point q, node a
    Real getN(uint q, uint a) {
      return _N[q*_nodesID.size() + a]; 
    }

    //! Get shape functions derivatives at quadrature point q, node a, direction i
    Real getDN(uint q, uint a, uint i) {
      return _DN[q*_nodesID.size() + a](i); 
    }

  protected:
    //! Quadrature weights
    vector<Real >          _QPweights;

    /*! 
      Shape functions values at quadrature points (Each entry corresponds to a 
      node, different quadrature points are stored subsequently)
    */
    vector<Real >          _N;

    /*! 
      Shape functions derivatives values at quadrature points (Each entry 
      corresponds to a node, different quadrature points are stored 
      subsequently)
    */
    vector<VectorXd >      _DN;

  }; // GeomElement

} // namespace voom

#endif
