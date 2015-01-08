#include "MechanicsLoad.h"

namespace voom
{

  void MechanicsLoad::checkConsistency(SourceLoadResults & R, vector<Real > u,
				       const Real h, const Real tol)
  {
    /*
    this->compute(R,u);
    Real Wan = R.W, Wplus = 0.0, Wminus = 0.0, Fnum = 0.0, error = 0.0;
    vector<Real > Fan = R.F;
    Matrix3d Pplus, Pminus;
    Pplus << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    Pminus = Pplus;
    FourthOrderTensor Kan = R.K;

    // First derivative test
    if( (R.request & ENERGY) && (R.request & FORCE) )
    {
      for (unsigned int i = 0; i < u.size(); i++) {
	u[i] += h;
	this->compute(R,u);
	Wplus = R.W;
	
	u[i] -= 2.0*h;
	this->compute(R,u);
	Wminus = R.W;
	
	u[i] += h;

	Fnum = (Wplus - Wminus)/(2.0*h);
	error += square(Fnum-Fan(i));
      } // i
      
      if (sqrt(error) < tol) {
	cout << "First derivative test PASSED with error = " << sqrt(error) << " . Source Energy = << " << Wan << endl;
      }
      else {
	cout << "First derivative test FAILED with error = " << sqrt(error) << " . Source Energy = << " << Wan << endl; 
      }

    } // First derivative test

   

    // Second derivative test
    error = 0.0;
    if( (R.request & FORCE) && (R.request & STIFFNESS) )
    {
      cout << "Second derivative test not implemented for MechanicsLoad. Follower forces not considered yet." << endl; 
    } // Second derivative test
   */ 
  } // checkConsistency
  
  
} // namespace voom
