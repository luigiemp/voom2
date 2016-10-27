//-*-C++-*-
#ifndef __BlankViscousPotential_h__
#define __BlankViscousPotential_h__

#include "VoomMath.h"
#include "ViscousPotential.h"

namespace voom
{
  class BlankViscousPotential : public ViscousPotential
  {
  public:
    BlankViscousPotential(){;}
    //! Destructor
    ~BlankViscousPotential(){;}

    virtual double phi(Matrix3d Fn, Matrix3d Fnp1, double deltaT){return 0.0;}
    virtual Matrix3d dphidF(Matrix3d Fn, Matrix3d Fnp1, double deltaT){return Matrix3d::Zero();}
    virtual FourthOrderTensor d2phidF2(Matrix3d Fn, Matrix3d Fnp1, double deltaT){return FourthOrderTensor(3,3,3,3);}

  }; // class BlankViscousPotential

} // namespace voom

#endif
