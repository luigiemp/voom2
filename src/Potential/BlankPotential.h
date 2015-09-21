//-*-C++-*-
#ifndef __BlankPotential_h__
#define __BlankPotential_h__

#include "VoomMath.h"
#include "Potential.h"

namespace voom
{
  class BlankPotential : public Potential
  {
  public:
    BlankPotential():Potential(){;}
    //! Destructor
    ~BlankPotential(){;}

    double Psi(Vector3d deltaQ, double deltaT){return 0;}
    Vector3d DPsiDQ(Vector3d deltaQ, double deltaT){return Vector3d::Zero();}
    Vector3d DPsiDQdot(Vector3d deltaQ, double deltaT){return Vector3d::Zero();}
    Matrix3d D2PsiDQDQ(Vector3d deltaQ, double deltaT){return Matrix3d::Zero();}
    Matrix3d D2PsiDQDQdot(Vector3d deltaQ, double deltaT){return Matrix3d::Zero();}
    Matrix3d D2PsiDQdotDQdot(Vector3d deltaQ, double deltaT){return Matrix3d::Zero();}


  }; // class BlankPotential

} // namespace voom

#endif
