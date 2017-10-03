#include "PairMaterial.h"

namespace voom
{

  void PairMaterial::checkConsistency(PairMresults & R, Real r,
					   Real h, Real tol)
  {
    this->compute(R, r);
    Real rloc = r, Wplus = 0.0, Wminus = 0.0, 
      Fan = R.F, Fnum = 0.0, Fplus = 0.0, Fminus = 0.0,
      Kan = R.K, Knum = 0.0, Kplus = 0.0, Kminus = 0.0, error = 0.0;

    // First derivative test
    if( (R.request & ENERGY) && (R.request & FORCE) )
    {
      rloc += h;
      this->compute(R, rloc);
      Wplus = R.W;
      
      rloc -= 2.0*h;
      this->compute(R, rloc);
      Wminus = R.W;
      
      rloc += h;
      
      Fnum = (Wplus - Wminus)/(2.0*h);
      error = fabs(Fnum-Fan);
      
      if (error < tol) {
	cout << "First derivative test PASSED with error = " << error << endl;
      }
      else {
	cout << "First derivative test FAILED with error = " << error << endl; 
      }

    } // First derivative test

   

    // Second derivative test
    error = 0.0;
    if( (R.request & FORCE) && (R.request & STIFFNESS) )
    {
      rloc += h;
      this->compute(R, rloc);
      Fplus = R.F;
      
      rloc -= 2.0*h;
      this->compute(R, rloc);
      Fminus = R.F;
      
      rloc += h;
      
      Knum = (Fplus - Fminus)/(2.0*h);
      error = fabs(Knum - Kan);
      
      if (error < tol) {
	cout << "Second derivative test PASSED with error = " << error << endl;
      }
      else {
	cout << "Second derivative test FAILED with error = " << error << endl; 
      }

    } // Second derivative test
           
  } // checkConsistency



} // namespace voom
