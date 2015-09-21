//-*-C++-*-
#ifndef __HillForceVelPotential_h__
#define __HillForceVelPotential_h__

#include "VoomMath.h"
#include "Potential.h"

namespace voom
{
  class HillForceVelPotential : public Potential
  {
  public:
    HillForceVelPotential(double a, double b, double S0):Potential(),_a(a), _b(b), _S0(S0){;}
    //! Destructor
    ~HillForceVelPotential(){;}

    double Psi(Vector3d deltaQ, double deltaT);
    Vector3d DPsiDQ(Vector3d deltaQ, double deltaT);
    Vector3d DPsiDQdot(Vector3d deltaQ, double deltaT);
    Matrix3d D2PsiDQDQ(Vector3d deltaQ, double deltaT);
    Matrix3d D2PsiDQDQdot(Vector3d deltaQ, double deltaT);
    Matrix3d D2PsiDQdotDQdot(Vector3d deltaQ, double deltaT);
  
  private:
    double _a;
    double _b;
    double _S0;

  }; // class HillForceVelPotential

} // namespace voom

#endif
