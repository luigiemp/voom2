#include <iostream>
#include <fstream>

#include "PassMyoA.h"
#include "CompNeoHookean.h"
#include "Humphrey.h"
#include "Humphrey_Compressible.h"
#include "LinYinActive.h"
#include "LinYinActive_Compressible.h"
#include "PlasticMaterial.h"
#include "HillForceVelPotential.h"
#include "BlankPotential.h"
#include "APForceVelPotential.h"
#include "Potential.h"
#include "ViscousPotential.h"
#include "BlankViscousPotential.h"
#include "NewtonianViscousPotential.h"
#include "FEMesh.h"
#include "MechanicsModel.h"
#include "EigenNRsolver.h"

using namespace voom;

int main(int argc, char** argv)
{
  cout << string(50, '\n'); // Clear Screen

  // Timing
  time_t start, end;
  time(&start);

  // Simulation Time (in ms)
  double simTime = 10;  
  
  // Time Step (in ms)
  double deltaT = 0.05/1000;

  // ActivationFile
  string ActivationFile = "ActivationFunction.dat";

  // writeOutputFlag set to true for output.
  bool writeOutputFlag = false;

  // OutputString
  string outputString = "Output/MaterialPoint/MatPoint.dat";
  ofstream outputFile;
  outputFile.open(outputString.c_str());
  if (!outputFile.is_open())
  {
    cout << "**ERROR: Could not open Output File: " << outputString << endl;
    exit(1);
  }

  // Newton Iteration Tolerance
  double TOL = 1.0e-9;

  // Valve opening threshold
  double Pthres = 17.0;

  // Viscous Parameter
  double eta = 10.0;

  // Force Velocity Potential Parameters
  double APForceVelPotential_a = 1.0;
  double APForceVelPotential_S0 = 100000.;

  // Start with valves open or closed?
  bool valvesOpen = true;

  // string outputString = "SingleElementTest_";


  // ***************************************** //
  // **************Begin Computing************ //
  // ***************************************** //
  cout << endl;
  cout << "\033[1;31mMaterial Point Compute Begin \033[0m" << endl;
  
  
   
  // Initialize Material
  
  vector <Vector3d> el_vectors(3, Vector3d::Zero(3,1));
  el_vectors[0] << 1., 0., 0.;
  el_vectors[1] << 0., 1., 0.;
  el_vectors[2] << 0., 0., 1.;

  // Humphrey PassiveMat(0, 15.98, 55.85, 0.0, -33.27, 30.21);
  Humphrey_Compressible PassiveMat(0, 15.98, 55.85, 0.0, -33.27, 30.21, 3.590, 64.62, el_vectors);
  LinYinActive_Compressible ActiveMat(0, -38.70, 40.83, 25.12, 9.51, 171.18, el_vectors);

  // CompNeoHookean PassiveMat(0, 4.0, 0.4);
  // CompNeoHookean ActiveMat(0, 4.0, 0.4);
  
  // BlankViscousPotential ViscPotential;
  NewtonianViscousPotential ViscPotential(0.0, 0.5);
  // High = 100, Med = 10, Low = 1

  Vector3d HardParam(0.,0.,0.);

  // PLmaterials.push_back(&PassiveMat);
  APForceVelPotential* TestPotential = new APForceVelPotential(APForceVelPotential_a, APForceVelPotential_S0); //1.0, 5.0);
  PlasticMaterial* PlMat = new PlasticMaterial(0, &ActiveMat, &PassiveMat, TestPotential, &ViscPotential);
  PlMat->setDirectionVectors(el_vectors);
  PlMat->setHardeningParameters(HardParam);
  PlMat->setActiveDeformationGradient(Matrix3d::Identity(3,3));
  PlMat->setTotalDeformationGradient(Matrix3d::Identity(3,3));

  PlMat->setTimestep(deltaT);
  PlMat->setActivationMultiplier(0.0);

  // Read in Activation File:
  ifstream myfile;
  myfile.open (ActivationFile.c_str());
  double tempTime = 0.0;
  double tempActivationFactor = 0.0;

  vector <double> Time;
  vector <double> ActivationFactor;

  while(myfile >> tempTime >> tempActivationFactor)
  {
    Time.push_back(tempTime);
    ActivationFactor.push_back(tempActivationFactor);
  }
  myfile.close();
  deltaT = Time[1] - Time[0];

  MechanicsMaterial::FKresults FKres;
  FKres.request = ENERGY | FORCE | STIFFNESS;
  
  // ************************************************************************ //
  // ************************************************************************ //
  // 1. Diastolic Stretching
  cout << "** Begin Uniaxial Diastolic Stretching" << endl;
  PlMat->setTimestep(deltaT/1000);
  PlMat->setActivationMultiplier(0.0);

  int numberOfDiastolicIncrements = 30;
  double Fstart = 0.0;
  double Fend = 0.2;

  double Fincrement = (Fend - Fstart)/(numberOfDiastolicIncrements - 1);

  Matrix3d F;
  Matrix3d ID;
  ID << 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0;

  for (int diasIter = 0; diasIter <= numberOfDiastolicIncrements; diasIter++)
  {
    cout << "Diastolic Iteration: " << diasIter << endl;


    Real stretchValue = diasIter * Fincrement;

    // This is done to make sure that at the last step, Fn = Fn+1 before viscosity is turned on
    if (diasIter == numberOfDiastolicIncrements)
	stretchValue = (diasIter - 1) * Fincrement;


    // Impose a P22 = P33 = 0
    F << ID(0,0) + stretchValue, 0.0, 0.0,
          0.0, 1.0, 0.0, //ID(1,1) + diasIter * Fincrement, 0.0,
          0.0, 0.0, 1.0; // 1.0/((ID(0,0) + diasIter * Fincrement) * (ID(1,1) + diasIter * Fincrement));
    double deltaLambda = 1.0;
    /* while(fabs(deltaLambda) > TOL)
    {
      PlMat->compute(FKres, F, NULL);
      deltaLambda = -1.0 * FKres.P(1,1)/(FKres.K(1,1,1,1) + FKres.K(1,1,2,2));
      F(1,1) += deltaLambda;
      F(2,2) += deltaLambda;
    }
    */

    // cout << "F = \n" << F << endl;
    PlMat->compute(FKres, F);
    Matrix3d ElasticStress = PlMat->getElasticStressOnly();
    Vector3d InternalVariables = PlMat->getCurrentHardeningParameters();

    outputFile << "0.0" << "\t" << F(0,0) << "\t" << F(1,1) << "\t" << F(2,2) << "\t" << ElasticStress(0,0) << "\t" << ElasticStress(1,1) << "\t" << ElasticStress(2,2) << "\t" << InternalVariables[0] << "\t" << InternalVariables[1] << "\t" << InternalVariables[2] << endl;

    PlMat->updateStateVariables();
  }

  // ************************************************************************ //
  // ************************************************************************ //
  // 2. Systole
  cout << "** SYSTOLE" << endl;
  
  Matrix3d Pelastic = PlMat->getElasticStressOnly();
  Vector3d InternalVariables = PlMat->getCurrentHardeningParameters();

  double P11elastic = Pelastic(0,0);
  ViscPotential.setViscosity(eta);

  PlMat->setTimestep(deltaT/1000);

  for (int sysIter = 0; sysIter < ActivationFactor.size(); sysIter++)
  {
    cout << "Systolic Iteration: \t" << sysIter << "\tActivation Factor: \t" << ActivationFactor[sysIter] << endl;
    PlMat->setActivationMultiplier(ActivationFactor[sysIter]);
    
    if (P11elastic >= Pthres && !valvesOpen)
    {
      valvesOpen = true;
      cout << "Valves Open!" << endl;
      ViscPotential.setViscosity(eta);
    }

    if (valvesOpen)
    {
      Real deltaLambda = 1.0;
      // Newton solve stress conditions
      while(abs(deltaLambda) > TOL)
      {
	PlMat->compute(FKres, F);
	deltaLambda = -1.0 * FKres.P(0,0) / FKres.K(0,0,0,0);
	F(0,0) += deltaLambda;
      }
    }
    // cout << "F = \n" << F << endl;
    PlMat->compute(FKres, F);
    Matrix3d ElasticStress = PlMat->getElasticStressOnly();
    InternalVariables = PlMat->getCurrentHardeningParameters();

    outputFile << Time[sysIter] << "\t" << F(0,0) << "\t" << F(1,1) << "\t" << F(2,2) << "\t" << ElasticStress(0,0) << "\t" << ElasticStress(1,1) << "\t" << ElasticStress(2,2) << "\t" << InternalVariables[0] << "\t" << InternalVariables[1] << "\t" << InternalVariables[2] << endl;

    PlMat->updateStateVariables();
    P11elastic = ElasticStress(0,0);
  }
  

  // Write Output
  
  outputFile.close();

  return 0;

}
