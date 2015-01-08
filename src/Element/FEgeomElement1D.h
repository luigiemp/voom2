//-*-C++-*-
/*!
  \file FEgeomElement1D.h

  \brief A derived class from GeomElement. This is for handling 1D elements
  which would be used in a 3D context.

*/
#if !defined(__FEgeomElement1D_h__)
#define __FEgeomElement1D_h__

#include  "GeomElement.h"

namespace voom{
  class FEgeomElement1D: public GeomElement {
  public:
   
    FEgeomElement1D(const int elemID, const vector<int > & nodesID, 
		    const vector<Vector3d > & nodesX, 
		    vector<Shape* > shape, Quadrature* quadrature,
		    const Real radius);

    //! Get number of quadrature points
    uint getNumberOfQuadPoints() {return _QPweights.size(); }

    //! Get weight for quadrature point q
    Real getQPweights(uint q) {return _QPweights[q]; }

    //! Get shape functions values at quadrature point q, node a
    Real getN(uint q, uint a) {
      return _N[q*_nodesID.size() + a]; 
    }

    //! Get shape functions derivatives at quadrature point q, node a, direction i
    Real getDN(uint q, uint a, uint i) {
      return (_DN[q*_nodesID.size() + a])(i); 
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
    vector<Vector3d >      _DN;

  }; // FEgeomElement1D

} // namespace voom

#endif
  
