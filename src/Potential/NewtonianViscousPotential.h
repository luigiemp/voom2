//-*-C++-*-
#ifndef __NewtonianViscousPotential_h__
#define __NewtonianViscousPotential_h__

#include "VoomMath.h"
#include "ViscousPotential.h"

namespace voom
{
  class NewtonianViscousPotential : public ViscousPotential
  {
  public:
    NewtonianViscousPotential(double eta, double alpha):_eta(eta), _alpha(alpha){;}
    //! Destructor
    ~NewtonianViscousPotential(){;}

    virtual double phi(Matrix3d Fn, Matrix3d Fnp1, double deltaT);
    virtual Matrix3d dphidF(Matrix3d Fn, Matrix3d Fnp1, double deltaT);
    virtual FourthOrderTensor d2phidF2(Matrix3d Fn, Matrix3d Fnp1, double deltaT);

    virtual void setViscosity(double eta){_eta = eta;}

  private:
    double _eta;
    double _alpha;

  }; // class NewtonianViscousPotential

} // namespace voom

#endif
