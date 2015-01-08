//-*-C++-*-
/*!
  \file MRKPMShape.h
  
  \brief Mesh free Modified Reproducing Kernel Particle Method
*/

#if !defined(__MRKPMShape_h__)
#define __MRKPMShape_h__

#include "MFShape.h"

namespace voom {
  enum kernelType { CONSTANT, LINEAR, CUBIC, QUARTIC, EXPONENTIAL };
  class MRKPMShape : public MFShape {
  public:
    //! Constructor compute N and DN
    MRKPMShape(const vector<VectorXd> & Nodes, const VectorXd & Point, 
	       const Real support, const Real radius, const Real supportHat);

    //! Update recomputes N and DN at new Point
    void update(const VectorXd & Point);

    // Accessors/Mutators to specific MRKPM class members
    Real getSupport() {return _support; }
    Real getRadius() {return _radius; }
    Real getSupportHat() {return _supportHat; }

    void setSupport(Real support) {_support = support;}
    void setRadius(const Real radius) { _radius = radius; }
    void setSupportHat(const Real supportHat) {_supportHat = supportHat; }

  private:
    // private functions for MRKPM calculation

    //! Compute kernel function value
    Real computePhi(Real z, const kernelType kernel = CUBIC);

    //! Compute Shapefuction
    void computeFunction(const VectorXd& Point);

    // extra members wrt classic FE shape functions
    //! Support value associated with the nodes
    Real _support;

    //! Radius
    Real _radius;

    //! Support Hat value
    uint _supportHat;
  };

} // namespace voom

#endif
