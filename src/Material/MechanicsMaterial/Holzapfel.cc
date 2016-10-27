#include "Holzapfel.h"

namespace voom {

  // Operators
  void Holzapfel::compute(FKresults & R, const Matrix3d & F)
  {
    // Needed for all requests
    Vector3d f, s;
    f = _fibers[0];
    s = _fibers[1];
    Matrix3d FinvT, Cbar, Mff, Mss, Mfs;
    FinvT = (F.inverse()).transpose();
    Real I3 = pow(F.determinant(), 2.0);
    Real I3oneThird = pow(I3, -1.0/3.0);
    Cbar = F.transpose()*F*I3oneThird;
    Real I1bar   = Cbar.trace();
    Real I4fbar  = f.dot(Cbar*f);
    Real I4sbar  = s.dot(Cbar*s);
    Real I8fsbar = f.dot(Cbar*s);   
    Mff = f*(f.transpose());
    Mss = s*(s.transpose());
    Mfs = f*(s.transpose()) + s*(f.transpose());

    if( R.request & ENERGY )
    {
      R.W = 0.5*( (_a1/_b1)*exp(_b1*(I1bar-3.0)) +
		  (_a2/_b2)*(exp(_b2*pow(I4fbar-1.0, 2.0))-1.0) +
		  (_a3/_b3)*(exp(_b3*pow(I4sbar-1.0, 2.0))-1.0) +
		  (_a4/_b4)*(exp(_b4*pow(I8fsbar,    2.0))-1.0) );
    }
    if( (R.request & FORCE) || (R.request & DMATPROP) ) 
    {
      R.P = _a1*exp(_b1*(I1bar-3.0))*(I3oneThird*F - (I1bar/3.0)*FinvT) + 
	2.0*( _a2*(I4fbar-1.0)*exp(_b2*pow(I4fbar-1.0, 2.0))*(I3oneThird*F*Mff - (I4fbar/3.0)*FinvT)   + 
	      _a3*(I4sbar-1.0)*exp(_b3*pow(I4sbar-1.0, 2.0))*(I3oneThird*F*Mss - (I4sbar/3.0)*FinvT) ) +
	_a4*I8fsbar*exp(_b4*pow(I8fsbar, 2.0))*(I3oneThird*F*Mfs - (2.0*I8fsbar/3.0)*FinvT);
    }
    if( R.request & STIFFNESS )
    {
      cout << " ***** To be coded ***** " << endl;
    } // STIFFNESS

    if( R.request & DMATPROP ) 
    {
      R.Dmat.resize(4, 3, 3);     // Already initialized to zero
      R.DDmat.resize(4, 4, 3, 3); // Already initialized to zero 
      Matrix3d Pa = exp(_b1*(I1bar-3.0))*(I3oneThird*F - (I1bar/3.0)*FinvT);
      Matrix3d Pb = 2.0*(I4fbar-1.0)*exp(_b2*pow(I4fbar-1.0, 2.0))*(I3oneThird*F*Mff - (I4fbar/3.0)*FinvT);
      Matrix3d Pc = 2.0*(I4sbar-1.0)*exp(_b3*pow(I4sbar-1.0, 2.0))*(I3oneThird*F*Mss - (I4sbar/3.0)*FinvT);
      Matrix3d Pd = I8fsbar*exp(_b4*pow(I8fsbar, 2.0))*(I3oneThird*F*Mfs - (2.0*I8fsbar/3.0)*FinvT);
      
      for (unsigned int i = 0; i<3; i++) {
	for (unsigned int J = 0; J<3; J++) {
	  (R.Dmat).set( 0, i, J, Pa(i,J) );
	  (R.Dmat).set( 1, i, J, Pb(i,J) );
	  (R.Dmat).set( 2, i, J, Pc(i,J) );
	  (R.Dmat).set( 3, i, J, Pd(i,J) );
	} // i
      } // J
    } // DMATPROP

  } // Holzapfel::compute

} // namespace voom


