#include <vector>
#include <iostream>

#include "CompNeoHookean.h"
#include "IsotropicDiffusion.h"

using namespace voom;

int main()
{
  cout << endl << "Testing mechanics material class ... " << endl;
  cout << ".................................... " << endl << endl;
  
  CompNeoHookean MatMech(1.0, 3.0);

  MechanicsMaterial::FKresults Rm;
  Rm.request = 7;
  Matrix3d F;
  F << 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0;
  srand(time(NULL));
  for (unsigned int i = 0; i<3; i++) 
    for (unsigned int J = 0; J<3; J++) 
      F(i,J) += 0.1*(double(rand())/RAND_MAX);
  cout << "determinant(F) = " << F.determinant() << endl;


  MatMech.compute(Rm, F);

  cout << "Energy     = " << Rm.W << endl;
  cout << "P(0,0)     = " << Rm.P(0,0) << endl;
  cout << "K[0,0,0,0] = " << Rm.K.get(0,0,0,0) << endl;

  MatMech.checkConsistency(Rm,F);



  cout << endl << "Testing diffusion material class ... " << endl;
          cout << ".................................... " << endl << endl;

  Real k = double(rand())/RAND_MAX;
  IsotropicDiffusion MatDiff(k);
  
  DiffusionMaterial::DiffusionResults Rd;
  MatDiff.compute(Rd);
  cout << "Conductivity = " << k << endl;
  cout << "A = " << Rd.A << endl;
 


  cout << endl << "....................................... " << endl;
          cout << "Test of voom material classes completed " << endl;
  return 0;
}
