#include <iostream>
#include "SCElastic.h"

namespace voom 
{
  void SCElastic::compute(Shellresults & R, const ShellGeometry & defGeom)
  {

    typedef vector<Vector3d> BasisVectors;

    const BasisVectors & basis = defGeom.a();
    const BasisVectors & dual = defGeom.aDual();
    const BasisVectors & dPartials = defGeom.dPartials();

    _H = - 0.5 * ( dual[0].dot( dPartials[0] ) + \
		   dual[1].dot( dPartials[1]) );
    
    const Matrix2d & aInv = defGeom.metricTensorInverse();
		
    // 		Tensor2D curvatureTensor;
    // 		curvatureTensor(0,0) = -dot(basis(0), dPartials(0));
    // 		curvatureTensor(0,1) = -dot(basis(0), dPartials(1));
    // 		curvatureTensor(1,0) = -dot(basis(1), dPartials(0));
    // 		curvatureTensor(1,1) = -dot(basis(1), dPartials(1));	       

    Matrix2d mixedCurvatureTensor;
    // b^alpha_beta = mixedCurvatureTensor(alpha,beta)
    // first index is upper, second is lower
    mixedCurvatureTensor(0,0) = -dual[0].dot( dPartials[0] );
    mixedCurvatureTensor(0,1) = -dual[0].dot( dPartials[1] );
    mixedCurvatureTensor(1,0) = -dual[1].dot( dPartials[0] );
    mixedCurvatureTensor(1,1) = -dual[1].dot( dPartials[1] );

    _K =  (  mixedCurvatureTensor(0,0) * mixedCurvatureTensor(1,1) 
	     - mixedCurvatureTensor(0,1) * mixedCurvatureTensor(1,0) );
    
    //
    //  edited by Feng on Dec. 7th, 2004
    //  Original code:
    //  double twoHminusC0 = 2.0*_H - _C0;
    //  modified code:
    //  Reason: we are using positive spontaneous curvatures (_C0), but a sphere has
    //          negative mean curvature (_H), -1/R, I believe that the operatore should be "+"
    double twoHminusC0 = 2.0*_H /* + */ - _C0;

    if( R.request & ENERGY ) {
      // compute strain energy
      R.W =  0.5 * _kC * twoHminusC0 * twoHminusC0;
      R.W += _kG * _K ;
      //R.W = 1.0;
      //cout << "H : " << _H << " K : "<< _K <<endl;
    }
		
    if( R.request & FORCE ) {
      //cout << "H : " << _H << " K : "<< _K <<endl;
      // -------------- stress resultants --------------
      R.n[0] = _kC * twoHminusC0 * 
	( aInv(0,0) * dPartials[0] + 
	  aInv(0,1) * dPartials[1] )
	+ 0.5 * _kC * twoHminusC0 * twoHminusC0 * dual[0] ;
      R.n[0] += _kG*_K*dual[0];

      R.n[1] = _kC * twoHminusC0 * 
	( aInv(1,0) * dPartials[0] + 
	  aInv(1,1) * dPartials[1] )
	+ 0.5 * _kC * twoHminusC0 * twoHminusC0 * dual[1] ;
      R.n[1] += _kG*_K*dual[1]; 

      for(int beta=0; beta<2; beta++) {
	for(int alpha=0; alpha<2; alpha++) {
	  R.n[beta] += _kG*2.0*_H*aInv(beta,alpha)*dPartials[alpha];
	  for(int nu=0; nu<2; nu++) {
	    R.n[beta] -= _kG*dual[beta].dot(dPartials[nu]) * 
	      mixedCurvatureTensor(nu,alpha)*dual[alpha];
	  }
	}
      }
      // 
      //R.n[2] << 0.0, 0.0, 0.0;
      //
      // ------------------- moment resultants ------------------
      R.m[0] = - _kC * twoHminusC0 * dual[0];
      R.m[1] = - _kC * twoHminusC0 * dual[1];

      for(int beta=0; beta<2; beta++) {
	R.m[beta] -= _kG*2.0*_H*dual[beta];
	for(int alpha=0; alpha<2; alpha++) {
	  R.m[beta] += _kG*mixedCurvatureTensor(beta,alpha)*dual[alpha];
	}
      }
    }

    return;
  }
}
