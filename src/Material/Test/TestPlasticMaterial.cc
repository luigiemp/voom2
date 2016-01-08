#include <vector>
#include <iostream>

#include "CompNeoHookean.h"
#include "IsotropicDiffusion.h"
#include "PassMyoA.h"
#include "PlasticMaterial.h"
#include "Potential.h"
#include "APForceVelPotential.h"
#include "BlankPotential.h"
#include "HillForceVelPotential.h"
#include "ViscousPotential.h"
#include "BlankViscousPotential.h"
#include "NewtonianViscousPotential.h"

using namespace voom;

int main()
{
  cout << endl << "Testing mechanics Plastic material class ... " << endl;
  cout << ".................................... " << endl << endl;
  
  cout << endl << "Passive Material: CompNeoHookean. " << endl;
  CompNeoHookean PassiveMat(0, 0.4, 0.04);
  cout << endl << "Active Material: CompNeoHookean. " << endl;
  CompNeoHookean ActiveMat(1, 40.0, 4.0);

  // cout << endl << "Creating a Blank Potential Function. " << endl;
  BlankPotential TestPotential;

  cout << endl << "Creating a Hill Force Velocity Potential Function. " << endl;
  // HillForceVelPotential TestPotential(4.4*pow(10,-3), .01*0.59, 25);
  // APForceVelPotential TestPotential(1.0, 250.0);

  BlankViscousPotential ViscPotential;
  // NewtonianViscousPotential ViscPotential(0.1, 0.5);

  PlasticMaterial PlMat(2, &ActiveMat, &PassiveMat, &TestPotential, &ViscPotential);
  vector <Vector3d> dirvec(3, Vector3d::Zero(3,1));
  dirvec[0] << 1, 0, 0;
  dirvec[1] << 0, 1, 0;
  dirvec[2] << 0, 0, 1;

  PlMat.setDirectionVectors(dirvec);
  cout << "Set DirVecs" << endl;

  Vector3d HardParam(0,0,0);
  PlMat.setHardeningParameters(HardParam);
  cout << "Set Hardening Parameters" << endl;
  PlMat.setActiveDeformationGradient(Matrix3d::Identity(3,3));
  cout << "Set Active Deformation Gradient at n" << endl;
  PlMat.setTotalDeformationGradient(Matrix3d::Identity(3,3));
  cout << "Set Total Deformation Gradient at n" << endl;

  PlMat.setTimestep(0.01);
  cout << "Set Timestep to 0.01s" << endl;

  PlMat.setActivationMultiplier(0.0);
  cout << "Set Activation Multiplier to 0.0" << endl;
  
  cout << endl << "Plastic Material Created. " << endl;

  MechanicsMaterial::FKresults Rm;
  Rm.request = 7;
  Matrix3d F;
  F << 1., 0.0, 0.0, 0.0, 1., 0.0, 0.0, 0.0, 1.;
  F << 1.0008, .0015, .0100, .0023, 1.0083, .0008, .0091, .0054, 1.0044;
  /*
  srand(time(NULL));
  for (unsigned int i = 0; i<3; i++) 
    for (unsigned int J = 0; J<3; J++) 
      F(i,J) += 0.1*(double(rand())/RAND_MAX); 
  */
  cout << "determinant(F) = " << F.determinant() << endl;
    

  PlMat.compute(Rm, F);
    
  cout << "Energy     = " << setprecision(15) << Rm.W << endl;
  cout << "P(2,2)     = " << setprecision(15) << Rm.P(2,2) << endl;
  cout << "K[0,0,0,0] = " << setprecision(15) << Rm.K.get(0,0,0,0) << endl;
  
  for (int i = 0; i < 3; i++){
	for (int j = 0; j < 3; j++) {
		for (int k = 0; k < 3; k++) {
			for (int l = 0; l < 3; l++) {
				cout << setprecision(15) << Rm.K.get(i,j,k,l) << endl;
			}
		}
	}
  }
    
  cout << endl << "Material ID = " << PlMat.getMatID() << endl << endl;
 
  PlMat.checkConsistency(Rm,F);
  
  return 0;
}
