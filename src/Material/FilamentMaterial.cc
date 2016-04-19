#include "FilamentMaterial.h"

namespace voom
{

  void FilamentMaterial::checkConsistency(Filresults & R,const Real h, const Real tol)
  {
    
    //this->compute(R,x,X);
    Real Wplus = 0.0, Wminus = 0.0, fnum = 0.0, error = 0.0, Kan = R.k, Knum = 0.0;
    Vector3d fplus;
    fplus << 0.0 , 0.0 , 0.0;
    Vector3d fminus = fplus;
    Vector3d fan = R.f1;
        
    vector<Vector3d> X;
    Vector3d X0;
    X0 << 0.0,0.0,0.0;
    Vector3d X1;
    X1 << 1.0,0.0,0.0;
    X.push_back(X0);
    X.push_back(X1);
    vector<Vector3d> x = X;

    //First derivative test
    if( (R.request & ENERGY) && (R.request & FORCE) )  
      {
	for (unsigned int i = 0; i<3; i++) {
	  x[0](0) += h;
	  this->compute(R,x,X);
	  Wplus = R.W;
	  x[0](0) -= 2.0*h;
	  this->compute(R,x,X);
	  Wminus = R.W;
	  x[0](0) += h;
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

                                                                                                                                                                                                                                          
    // Second derivative test
    error = 0.0;
    if( (R.request & FORCE) && (R.request & STIFFNESS) )
    {   
      Vector3d d = x[1] - x[0];
      for (unsigned int i = 0; i<3; i++){
	x[0](0) += h*d(i)/d.norm();
      }
      
      this->compute(R,x,X);
      fplus = R.f1;
  
      for (unsigned int i = 0; i<3; i++){
        x[0](0) -= 2*h*d(i)/d.norm();
      }
      d = x[1]- x[0];

      this->compute(R,x,X);
      fminus = R.f1;
      for (unsigned int i = 0; i<3; i++){
        x[0](0) += h*d(i)/d.norm();
      }
      d = x[1]- x[0];

      Knum = (fplus.norm()-fminus.norm())/(2.0*h);
      error += square(Knum - Kan);
      
      if (sqrt(error) < tol) {
        cout << "Second derivative test PASSED with error = " << sqrt(error) << endl;
      }
      else {
        cout << "Second derivative test FAILED with error = " << sqrt(error) << endl;
      }

    }//Second derivative Test
    

  } // checkConsistency


} // namespace voom
