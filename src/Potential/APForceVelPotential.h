//-*-C++-*-
#ifndef __APForceVelPotential_h__
#define __APForceVelPotential_h__

#include "VoomMath.h"
#include "Potential.h"

namespace voom
{
  class APForceVelPotential : public Potential
  {
  public:
    APForceVelPotential(double a, double S0, double v0=3.0):Potential(),_a(a), _S0(S0),_v0(v0){;}
    //! Destructor
    ~APForceVelPotential(){;}

    double Psi(Vector3d deltaQ, double deltaT);
    Vector3d DPsiDQ(Vector3d deltaQ, double deltaT);
    Vector3d DPsiDQdot(Vector3d deltaQ, double deltaT);
    Matrix3d D2PsiDQDQ(Vector3d deltaQ, double deltaT);
    Matrix3d D2PsiDQDQdot(Vector3d deltaQ, double deltaT);
    Matrix3d D2PsiDQdotDQdot(Vector3d deltaQ, double deltaT);
  
  private:
    double _a;
    double _S0;
    double _v0;

  }; // class APForceVelPotential

} // namespace voom

#endif
