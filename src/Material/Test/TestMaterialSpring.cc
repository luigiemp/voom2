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
    d0[0] = 1.0;
    d0[1] = 0.0;
    d0[2] = 0.0;
    
    Spring MatSpring(0, k, d0);
    
    MechanicsMaterial::Filresults Rf;
    
    
    Rf.request = 15;
    /*
    Matrix3d F;
    F << 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0;
    srand(time(NULL));
    for (unsigned int i = 0; i<3; i++) 
      for (unsigned int J = 0; J<3; J++) 
	F(i,J) += 0.1*(double(rand())/RAND_MAX);
    cout << "determinant(F) = " << F.determinant() << endl;
    

    MatMech.compute(Rm, F);
    
    cout << "Energy     = " << Rm.W << endl;
    cout << "P(2,2)     = " << Rm.P(2,2) << endl;
    cout << "K[0,0,0,0] = " << Rm.K.get(0,0,0,0) << endl;
    // for (unsigned int i = 0; i<3; i++) {
    //   for (unsigned int J = 0; J<3; J++) {
    // 	cout << i << " " << J << " " << (Rm.Dmat).get( 0, i, J ) << " " << (Rm.Dmat).get( 1, i, J ) << endl;
    //   }
    // }
    
    cout << endl << "Material ID = " << MatMech.getMatID() << endl << endl;
 
    MatMech.checkConsistency(Rm,F);
    */
  }


  return 0;
}
