#if !defined(__CahnHilliard_h__)
#define __CahnHilliard_h__

#include "VoomMath.h"
#include "ShellGeometry.h"

namespace voom
{
  
  class CahnHilliard
  {

  public:
    struct Scalarresults {
      Scalarresults(){
	W = 0.0;
	n = 0.0;
	m.resize(2,0);
	constr = 0;
	request = 0;
      };
      
      double W; //Energy Density
      double n; //coefficient of variation in conc
      vector<double> m; //coefficient of variation in conc_alpha
      Real constr;
      int request;
    };

  CahnHilliard(){;};
  CahnHilliard(int MatID): _matID(MatID) {;}
  CahnHilliard(Real eps): _eps(eps) {;}
  void compute(Scalarresults & R, Real phi, 
	       Vector2d phiPartials, Real lagMult, ShellGeometry& geom);

  void setEpsilon(double spe) {_eps = 4.0/(spe*(spe+1));}
  double UU(Real phi){ return pow(1-phi*phi,2.0); }
  double DUU(Real phi){ return -4*phi*(1-phi*phi); }

  protected:

    Real _eps;
    int _matID;


  };
  
} //namespace voom

#endif //  !defined(__SCElastic_h__)
