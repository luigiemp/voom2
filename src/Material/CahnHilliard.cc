#include "CahnHilliard.h"
#include <iostream>

namespace voom 
{
  void CahnHilliard::compute(Scalarresults & R, Real phi, 
			     Vector2d phiPartials, Real lagMult, ShellGeometry & geom)
  {
    Matrix2d gab = geom.metricTensorInverse();
    Real mu = 0;
    if( R.request & ENERGY ) {
      R.W = 0.5 * _eps * (gab(0,0)*phiPartials(0)*phiPartials(0) +
			  2*gab(0,1)*phiPartials(0)*phiPartials(1)+
			  gab(1,1)*phiPartials(1)*phiPartials(1) )+ UU(phi) + lagMult*(phi-mu);
    }
		
    if( R.request & FORCE ) {
      R.n = DUU(phi) + lagMult; //Coefficient of delta_phi 
      R.m[0] = 2* _eps*(gab(0,0)*phiPartials(0) + gab(0,1)*phiPartials(1) ); //Coeff of delta_phi_0
      R.m[1] = 2* _eps*(gab(1,0)*phiPartials(0) + gab(1,1)*phiPartials(1) ); //Coeff of delta_phi_1
      R.constr = (phi-mu);
    return;
  }
}

}
