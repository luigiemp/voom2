#if !defined(__LB_h__)
#define __LB_h__

#include "VoomMath.h"
#include "ShellGeometry.h"

namespace voom
{
  
  class LandauBrazovskii
  {

  public:
    struct Scalarresults {
      Scalarresults(){
	W = 0.0;
	n1 = 0.0;
	n2 = Vector2d::Zero();
	n3 = 0.0;
	request = 0;
      };
      
      double W; //Energy Density
      double n1; //coefficient of variation in conc
      Vector2d n2; //coefficient of variation in conc_alpha
      double n3; //coefficient of variation in Lap(conc)
      int request;
    };

  LandauBrazovskii(){;};
  LandauBrazovskii(int MatID): _matID(MatID) {;}
  LandauBrazovskii(Real c, Real k0, Real r, Real alpha, Real beta): _c(c), _k0(k0), _r(r), _alpha(alpha), _beta(beta) {;}
  void compute(Scalarresults & R, Real psi, 
	       Vector2d psiPartials, Matrix2d psi2Partials, ShellGeometry& geom);

  void setR(Real r) {_r = r;}
  void setC(Real c) {_c = c;}
  double UU(Real psi){ return _r*pow(psi,2.0)-_alpha*pow(psi,3.0)+_beta*pow(psi,4.0); }
  double DUU(Real psi){ return 2* _r*psi - 3*_alpha*pow(psi,2.0)+4*_beta*pow(psi,3.0);}

  protected:

    Real _k0;
    Real _c;
    Real _r;
    Real _alpha;
    Real _beta;
    int _matID;


  };
  
} //namespace voom

#endif 
