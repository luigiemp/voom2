//-*-C++-*-
/*!
  \file Shape.h
  
  \brief Base class implementation for FE and meshfree shape functions.
  
*/

#if !defined(__Shape_h__)
#define __Shape_h__

#include "voom.h"
#include "VoomMath.h"

namespace voom {

  class Shape {
  public:
    //! Constructor.
    Shape() {;}

    //! Destructor
    ~Shape() {;}

    //! Update recomputes N and DN at new Point
    virtual void update(const VectorXd & Point) = 0;

    //! Get number of shape functions
    virtual uint getShapeFunctionNum() = 0;

    //! getN returns shape function value at a given node a.
    virtual Real getN(const uint a) = 0;
    
    //! GetDN returns shape function derivatives at node a in direction i.
    virtual Real getDN(const uint a, const uint i) = 0;

    /*! 
      Consistency check. The shape function derivatives are computed 
      numerically and compared against analytical derivatives. Numerical 
      derivatives are computed as
      \f[
      \frac{\partial N(\xi)}{\partial \xi_i} \approx \frac{ N(\xi + \Delta
      \xi_i) - N(\xi - \Delta \xi_i) }{2\Delta \xi_i}
      \f]
      The \f$L_2\f$ error norm between numerical and analytical derivatives
      is computed.
    */

    bool checkConsistency(VectorXd Point, 
			  const Real eps = 1.0e-8, const Real tol = 1.0e-7);

    bool checkPartitionUnity(VectorXd Point, const Real tol = 1.0e-7);

  }; // class Shape

} // namespace voom

#endif
