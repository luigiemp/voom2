#include "FilamentMaterial.h"

namespace voom
{

  void FilamentMaterial::checkConsistency(Filresults & R, const Vector3d & d,
					   const Real h, const Real tol)
  {
    this->compute(R,d);
    Real Wplus = 0.0, Wminus = 0.0, fnum = 0.0, error = 0.0;
    Vector3d fplus;
    fplus << 0.0 , 0.0 , 0.0;
    Vector3d fmin = fplus;
    Vector3d fan = R.f;
    Vector3d dloc = d;

    //First derivative test
    if( (R.request & ENERGY) && (R.request & FORCE) )  
      {
	for (unsigned int i = 0; i<3; i++) {
	  dloc(i) += h;
	  this->compute(R,dloc);
	  Wplus = R.W;
	  dloc(i) -= 2.0*h;
	  this->compute(R,dloc);
	  Wminus = R.W;
	  dloc(i) += h;
	  fnum = (Wplus - Wminus)/(2.0*h);
	  error += square((fnum-fan(i)));
	}
	if (sqrt(error) < tol) { 
	  cout << "First derivative test PASSED with error = " << sqrt(error) << endl; 
	}
	else {
	  cout << "First derivative test FAILED with error = " << sqrt(error) << endl;   
	}
      }

   
    /*
    // Second derivative test
    error = 0.0;
    if( (R.request & FORCE) && (R.request & STIFFNESS) )
    {
      for (unsigned int i = 0; i<3; i++) {
	for (unsigned int J = 0; J<3; J++) {
	  for (unsigned int k = 0; k<3; k++) {
	    for (unsigned int L = 0; L<3; L++) {
    
                Floc(k,L) += h;
		this->compute(R,Floc,Fiber);
		Pplus = R.P;

		Floc(k,L) -= 2.0*h;
		this->compute(R,Floc,Fiber);
		Pminus = R.P;

		Floc(k,L) += h;
                
                Cnum = (Pplus(i,J) - Pminus(i,J))/(2.0*h);
		error += square(Cnum - Kan.get(i,J,k,L));

	    } // i
	  } // J
	} // k
      } // L

      if (sqrt(error) < tol) {
	cout << "Second derivative test PASSED with error = " << sqrt(error) << endl;
      }
      else {
	cout << "Second derivative test FAILED with error = " << sqrt(error) << endl; 
      }

    } // Second derivative test
    */

  } // checkConsistency


} // namespace voom
