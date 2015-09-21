//-*-C++-*-
#ifndef __ViscousPotential_h__
#define __ViscousPotential_h__

#include "VoomMath.h"


namespace voom
{
  class ViscousPotential
  {
  public:
    ViscousPotential(){;}
    //! Destructor
    ~ViscousPotential(){;}

    virtual double phi(Matrix3d Fn, Matrix3d Fnp1, double deltaT) = 0;
    virtual Matrix3d dphidF(Matrix3d Fn, Matrix3d Fnp1, double deltaT) = 0;
    virtual FourthOrderTensor d2phidF2(Matrix3d Fn, Matrix3d Fnp1, double deltaT) = 0;

    virtual void setViscosity(double eta){;}

  }; // class Viscous Potential

} // namespace voom

#endif
