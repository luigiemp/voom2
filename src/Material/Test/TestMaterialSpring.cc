#include <vector>
#include <iostream>

#include "Spring.h"

using namespace voom;

int main()
{
  cout << endl << "Testing mechanics material class ... " << endl;
  cout << ".................................... " << endl << endl;
  
  {
    cout << endl << "Testing Spring material. " << endl;
    Real k = 1.0;
    Vector3d d0;
    Vector3d d;
    d0 << 1.0, 0.0, 0.0;
    d << 1.1,0.0,0.0;
     
    Spring MatSpring(0, k, d0);
     
    FilamentMaterial::Filresults Rf;
    
    Rf.request = 7;

    MatSpring.compute(Rf,d);

    cout << "Energy     = " << Rf.W << endl;    
    cout << "Force  = " << Rf.f << endl;                                                                                                                                                                                             
    cout << "Stiffness = " << Rf.k << endl;   
    
    cout << endl << "Material ID = " << MatSpring.getMatID() << endl << endl;
 
    MatSpring.checkConsistency(Rf,d);

  }


  return 0;
}
