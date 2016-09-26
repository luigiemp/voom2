#include "Guccione.h"

namespace voom {

  // Operators
  void Guccione::compute(FKresults & R, const Matrix3d & F)
  {
    // Needed for all requests
    Vector3d f, c, r;
    f = _fibers[0];
    c = _fibers[1];
    r = _fibers[2];
    Matrix3d E, ID, Fff, Fcc, Frr, Fcr, Frc, Ffc, Fcf, Ffr, Frf;
    ID << 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0;
    E = 0.5*(F.transpose()*F - ID);
    Real Eff  = f.dot(E*f);
    Real Ecc  = c.dot(E*c);
    Real Err  = r.dot(E*r);
    Real Ecr  = c.dot(E*r);
    Real Efc  = f.dot(E*c);
    Real Efr  = f.dot(E*r);

    Fff  = F*f*(f.transpose());
    Fcc  = F*c*(c.transpose());
    Frr  = F*r*(r.transpose());
    Fcr  = F*c*(r.transpose());; Frc  = F*r*(c.transpose());;
    Ffc  = F*f*(c.transpose());; Fcf  = F*c*(f.transpose());;
    Ffr  = F*f*(r.transpose());; Frf  = F*r*(f.transpose());;

    // ENERGY
    R.W = _a1*exp( _b1*Eff*Eff + _b2*(Ecc*Ecc + Err*Err + 2*Ecr*Ecr) + 2*_b3*(Efc*Efc + Efr*Efr) );

    if( (R.request & FORCE) || (R.request & DMATPROP) || (R.request & STIFFNESS) ) 
    {
      R.P = R.W * ( 2.0*_b1*Eff*Fff + 2.0*_b2*( Ecc*Fcc + Err*Frr + Ecr*(Fcr+Frc) ) + _b3*( 2.0*Efc*(Fcf+Ffc) + 2.0*Efr*(Ffr+Frf) ) );
    }

    if( R.request & STIFFNESS )
    {
      R.K.resetIterator();
      for (unsigned int L = 0; L<3; L++) {
	for (unsigned int k = 0; k<3; k++) {
	  for (unsigned int J = 0; J<3; J++) {
	    for (unsigned int i = 0; i<3; i++) {
	      R.K.sequentialSet( R.P(k,L)*R.P(i,J)/R.W + 
				 R.W*( 2.0*_b1*( Fff(k,L)*Fff(i,J) + Eff*ID(i,k)*f(L)*f(J) ) +
				       2.0*_b2*( Fcc(k,L)*Fcc(i,J) + Ecc*ID(i,k)*c(L)*c(J) +
                                                 Frr(k,L)*Frr(i,J) + Err*ID(i,k)*r(L)*r(J) +
                                                 0.5*( Fcr(k,L) + Frc(k,L) ) * ( Fcr(i,J) + Frc(i,J) ) +
						 Ecr*ID(i,k)*( c(L)*r(J) + r(L)*c(J) ) ) +
				       _b3*( ( Fcf(k,L) + Ffc(k,L) )*( Fcf(i,J) + Ffc(i,J) ) +
                                             2.0*Efc*ID(i,k)*( f(L)*c(J) + c(L)*f(J) ) +
                                             ( Ffr(k,L) + Frf(k,L) )*( Ffr(i,J) + Frf(i,J) ) +
                                             2.0*Efr*ID(i,k)*( f(L)*r(J) + r(L)*f(J) ) ) )
				 );
	      R.K.incrementIterator();
	    } // L
	  } // k
	} // J
      } // i
    } // STIFFNESS



    if( R.request & DMATPROP ) 
    {
      R.Dmat.resize(1, 3, 3);     // Already initialized to zero
      R.DDmat.resize(1, 1, 3, 3); // Already initialized to zero 
      Matrix3d Pa = (1.0/_a1)*R.P;
      
      for (unsigned int i = 0; i<3; i++) {
	for (unsigned int J = 0; J<3; J++) {
	  (R.Dmat).set( 0, i, J, Pa(i,J) );
	} // i
      } // J
    } // DMATPROP

  } // Guccione::compute

} // namespace voom


