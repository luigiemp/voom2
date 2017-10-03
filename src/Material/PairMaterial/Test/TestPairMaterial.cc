#include <vector>
#include <iostream>

#include "Harmonic.h"
#include "Quartic.h"

using namespace voom;

int main()
{
  cout << endl << "Testing PairMaterial class ... " << endl;
  
  {
    cout << ".................................... " << endl << endl;
    cout << endl << "Testing Harmonic material. " << endl;

    srand(time(NULL));
    Real k  = double(rand())/RAND_MAX;
    Real r0 = double(rand())/RAND_MAX;
    Harmonic PairMat( k, r0 );
    
    PairMaterial::PairMresults R;
    R.request = (ENERGY | FORCE | STIFFNESS);
    
    Real r = double(rand())/RAND_MAX;
    
    PairMat.compute(R, r);
    
    cout << "Energy     = " << R.W << endl;
    cout << "Force      = " << R.F << endl;
    cout << "Stiffness  = " << R.K << endl;
 
    PairMat.checkConsistency(R, r);
  }


  {
    cout << ".................................... " << endl << endl;
    cout << endl << "Testing Quartic material. " << endl;

    srand(time(NULL));
    Real k  = double(rand())/RAND_MAX;
    Real r0 = double(rand())/RAND_MAX;
    Quartic PairMat( k, r0 );
    
    PairMaterial::PairMresults R;
    R.request = (ENERGY | FORCE | STIFFNESS);
    
    Real r = double(rand())/RAND_MAX;
    
    PairMat.compute(R, r);
    
    cout << "Energy     = " << R.W << endl;
    cout << "Force      = " << R.F << endl;
    cout << "Stiffness  = " << R.K << endl;
 
    PairMat.checkConsistency(R, r);
  }



  cout << endl << "....................................... " << endl;
  cout << "Test of PairMaterial classes completed " << endl;
  

  
  return 0;
}
