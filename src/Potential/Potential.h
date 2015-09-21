//-*-C++-*-
#ifndef __Potential_h__
#define __Potential_h__

#include "VoomMath.h"


namespace voom
{
  class Potential
  {
  public:
    Potential(){;}
    //! Destructor
    ~Potential(){;}

    virtual double Psi(Vector3d deltaQ, double deltaT) = 0;
    virtual Vector3d DPsiDQ(Vector3d deltaQ, double deltaT) = 0;
    virtual Vector3d DPsiDQdot(Vector3d deltaQ, double deltaT) = 0;
    virtual Matrix3d D2PsiDQDQ(Vector3d deltaQ, double deltaT) = 0;
    virtual Matrix3d D2PsiDQDQdot(Vector3d deltaQ, double deltaT) = 0;
    virtual Matrix3d D2PsiDQdotDQdot(Vector3d deltaQ, double deltaT) = 0;


  }; // class Potential

} // namespace voom

#endif
