#include "LandauBrazovskii.h"
#include <iostream>

namespace voom 
{
  void LandauBrazovskii::compute(Scalarresults & R, Real psi, 
				 Vector2d psiPartials, Matrix2d psi2Partials, ShellGeometry & geom)
  {
    Matrix2d gab = geom.metricTensorInverse();
    Real LapPsi = (gab*psi2Partials).trace() - psiPartials.dot(geom.gklGammaI_kl());
    Real lambda = 0;
    Real c = _c;
    
    //if (abs(LapPsi/psi+20) > 3) cout <<"LapPsi/psi : " << LapPsi/psi << endl << "---------------------" << endl;
    if( R.request & ENERGY ) {
      
      R.W =  c*(LapPsi*LapPsi +2*_k0*_k0*LapPsi*psi 
		+ pow(_k0*_k0*psi,2.0) )+ UU(psi)+ lambda*psi ; 
      
      /*
      // ------------ without IBP ----------------
      R.W =  LapPsi*LapPsi - 2*_k0*_k0*psiPartials.dot(gab*psiPartials) 
	+ pow(_k0*_k0*psi,2.0)+ UU(psi); 
      */
    }
		
    if( R.request & FORCE ) {
      
      R.n1 = c*(2*pow(_k0,4.0)*psi + 2*pow(_k0,2.0)*LapPsi) + DUU(psi) + lambda; //Coefficient of delta_psi 

      R.n2 =  Vector2d::Zero();

      R.n3 = c*(2*LapPsi + 2*pow(_k0,2.0)*psi);
      
      
      /*
      // ------------ without IBP ----------------
      R.n1 = 2*pow(_k0,4.0)*psi + DUU(psi) ; //Coefficient of delta_psi 

      R.n2 = -4*_k0*_k0*gab*psiPartials;

      R.n3 = 2*LapPsi;
      */
    return;
  }
}

}
