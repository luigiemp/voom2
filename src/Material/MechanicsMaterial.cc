#include "MechanicsMaterial.h"

namespace voom
{

void MechanicsMaterial::checkConsistency(FKresults & R, const Matrix3d & F,
					  const Real h, const Real tol)
  {
    this->compute(R,F);
    Real Wplus = 0.0, Wminus = 0.0, Pnum = 0.0, Cnum = 0.0, error = 0.0;
    Matrix3d Pan = R.P, Floc = F, Pplus, Pminus;
    Pplus << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    Pminus = Pplus;
    FourthOrderTensor Kan = R.K;

    // First derivative test
    if( (R.request & ENERGY) && (R.request & FORCE) )
    {
      for (unsigned int i = 0; i<3; i++) {
	for (unsigned int J = 0; J<3; J++) {

	  Floc(i,J) += h;
	  this->compute(R,Floc);
	  Wplus = R.W;

	  Floc(i,J) -= 2.0*h;
	  this->compute(R,Floc);
	  Wminus = R.W;

	  Floc(i,J) += h;

	  Pnum = (Wplus - Wminus)/(2.0*h);
	  error += square(Pnum-Pan(i,J));
	} // i
      } // J

      if (sqrt(error) < tol) {
	cout << "First derivative test PASSED with error = " << sqrt(error) << endl;
      }
      else {
	cout << "First derivative test FAILED with error = " << sqrt(error) << endl; 
      }

    } // First derivative test

   

    // Second derivative test
    error = 0.0;
    if( (R.request & FORCE) && (R.request & STIFFNESS) )
    {
      for (unsigned int i = 0; i<3; i++) {
	for (unsigned int J = 0; J<3; J++) {
	  for (unsigned int k = 0; k<3; k++) {
	    for (unsigned int L = 0; L<3; L++) {
    
                Floc(k,L) += h;
		this->compute(R,Floc);
		Pplus = R.P;

		Floc(k,L) -= 2.0*h;
		this->compute(R,Floc);
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
           
  } // checkConsistency


} // namespace voom
