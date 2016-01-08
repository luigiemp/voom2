#include "PassMyoA.h"
#include "CompNeoHookean.h"
#include "PlasticMaterial.h"
#include "HillForceVelPotential.h"
#include "BlankPotential.h"
#include "APForceVelPotential.h"
#include "Potential.h"
#include "ViscousPotential.h"
#include "BlankViscousPotential.h"
#include "NewtonianViscousPotential.h"
#include "FEMesh.h"
#include "EigenEllipticResult.h"
#include "MechanicsModel.h"
#include "EigenNRsolver.h"
#include "LBFGSB.h"

using namespace voom;

int main(int argc, char** argv)
{
  cout << string(50, '\n'); // Clear Screen
  
  // Timing
  time_t start, end;
  time(&start);

  // Conduction Velocity
  double cv = 0.05;

  // OutputString
  string outputString = "Cube_PV_HighViscosity";



  // Initialize Mesh
  // Assumptions to use this main as is: strip has a face at z=0; tetrahedral mesh
  FEMesh Cube("Cube1.node", "Cube1.ele");
  FEMesh surfMesh("Cube1.node", "Cube1Surf.ele");
  Real xmax = 1.0;
  // FEMesh Cube("Strip36.node", "Strip36.ele");
  // FEMesh surfMesh("Strip36.node", "Strip36Surf.ele");
  // FEMesh Cube("Strip144.node", "Strip144.ele");
  // FEMesh surfMesh("Strip144.node", "Strip144Surf.ele");
  // Real xmax = 6.0; 

  cout << endl;
  cout << "\033[1;32mNumber Of Nodes \t : \033[0m" << Cube.getNumberOfNodes() << endl;
  cout << "\033[1;32mNumber Of Element \t : \033[0m" << Cube.getNumberOfElements() << endl;
  cout << "\033[1;32mMesh Dimension \t\t : \033[0m" << Cube.getDimension() << endl << endl;
  

    
  // Initialize Material
  uint NumMat =  Cube.getNumberOfElements();
  
  vector<MechanicsMaterial * > PLmaterials;
  PLmaterials.reserve(NumMat);

  CompNeoHookean PassiveMat(0, 4.0, 0.4);
  CompNeoHookean ActiveMat(0, 4.0, 0.4);
  APForceVelPotential TestPotential(1.0, 500.0);
  // HillForceVelPotential TestPotential(4.4*pow(10,-3), .01*0.59, 25);
  // BlankPotential TestPotential;

  // BlankViscousPotential ViscPotential;
  NewtonianViscousPotential ViscPotential(0.05, 0.5);

  vector <Vector3d> dirvec(3, Vector3d::Zero(3,1));
  dirvec[0] << 1., 0., 0.;
  dirvec[1] << 0., 1., 0.;
  dirvec[2] << 0., 0., 1.;

  Vector3d HardParam(0.,0.,0.);

  // Read in Activation File:
  ifstream myfile;
  myfile.open ("ActivationFunction_LargeTimestep.dat");
  vector <double> Time(795, 0.0);
  vector <double> ActivationFactor(795, 0.0);
  for (int i = 0; i < 795; i++)
  {
    myfile >> Time[i];
    myfile >> ActivationFactor[i];
  }
  myfile.close();



  for (int k = 0; k < NumMat; k++) {
    // PLmaterials.push_back(&PassiveMat);
    PlasticMaterial* PlMat = new PlasticMaterial(k, &ActiveMat, &PassiveMat, &TestPotential, &ViscPotential);
    PlMat->setDirectionVectors(dirvec);
    PlMat->setHardeningParameters(HardParam);
    PlMat->setActiveDeformationGradient(Matrix3d::Identity(3,3));
    PlMat->setTotalDeformationGradient(Matrix3d::Identity(3,3));

    PlMat->setTimestep(0.01);
    PlMat->setActivationMultiplier(0.0);

    PLmaterials.push_back(PlMat);
  }

  // Initialize Model
  int NodeDoF = 3;
  int PressureFlag = 0;
  Real Pressure = 0.0;
  int NodalForcesFlag = 0;
  vector<int > ForcesID;
  vector<Real > Forces;
  MechanicsModel myModel(&Cube, PLmaterials, NodeDoF, PressureFlag, &surfMesh,
			 NodalForcesFlag);
  myModel.updatePressure(Pressure);
  myModel.updateNodalForces(&ForcesID, &Forces);
 
  // Initialize Result
  uint myRequest;
  uint PbDoF = (Cube.getNumberOfNodes())*myModel.getDoFperNode();
  EigenEllipticResult myResults(PbDoF, NumMat*2);

  // Run Consistency check

  Real perturbationFactor = 0.1;
  myRequest = 7; // Check both Forces and Stiffness
  Real myH = 1e-6;
  Real myTol = 1e-7;

  // Before checking consistency, the perturbed deformation state must be 
  // set to the current deformation state.

  // myModel.checkConsistency(myResults, perturbationFactor, myRequest, myH, myTol);
  // myModel.checkDmat(myResults, perturbationFactor, myH, myTol);





    // Print initial configuration
    myModel.writeOutputVTK("CubeSmall_", 0);

    // Check on applied pressure
    myRequest = 2;
    myResults.setRequest(myRequest);
    myModel.compute(myResults);
    VectorXd Fstart = *(myResults._residual);
    Real pX = 0.0, pY = 0.0, pZ = 0.0;
    for (int i = 0; i < Fstart.size(); i += 3) {
      pX += Fstart(i);
      pY += Fstart(i+1);
      pZ += Fstart(i+2);
    }
    cout << endl << "Pressure = " << pX << " " << pY << " " << pZ << endl << endl;
     
    // EBC
    vector<int > DoFid;
    vector<Real > DoFvalues;
    vector<int > DoFxmax;
    int ind = 0;
    for(int i = 0; i < Cube.getNumberOfNodes(); i++) {
      // All nodes at x = 0
      if ( Cube.getX(i, 0) < 1.0e-12 ) {
	DoFid.push_back(i*3);
	ind++;
	// All nodes at y = 0
	if ( Cube.getX(i, 1) < 1.0e-12 ) {
	  DoFid.push_back(i*3 + 1);
	  ind++;
	}
	// All nodes at z = 0
	if ( Cube.getX(i, 2) < 1.0e-12 ) {
	  DoFid.push_back(i*3 + 2);
	  ind++;
	}
      }

      // All nodes at x = xmax
      if ( Cube.getX(i, 0) > xmax - 1.0e-12 ) {
	DoFid.push_back(i*3);
	DoFxmax.push_back(ind);
	ind++;
      }
    
    }
 
    for(int i = 0; i < DoFid.size(); i++) {
      DoFvalues.push_back( Cube.getX(floor(double(DoFid[i])/3.0) ,DoFid[i]%3) );
      // cout << DoFid[i] << " " <<  DoFvalues[i] << endl;
    }

    // Lin tet cube - applied displacements //
    // for (uint i = 18; i < 27; i++) { 
    //   DoFid[i-5] = i*3 + 2; // Lower top in z direction by 0.1
    //   DoFvalues[i-5] = 1.95;
    // }
    
    // Solver
    Real NRtol = 1.0e-12;
    uint NRmaxIter = 100;
    EigenNRsolver mySolver(&myModel, DoFid, DoFvalues, CHOL, NRtol, NRmaxIter);

    ind = 0;
    int NumPassiveInc = 00;
    Real DeltaX = 0.2*xmax/NumPassiveInc;
    // cout << DoFxmax.size() << endl;
    for (int s = 0; s < NumPassiveInc; s++) {
      // Impose BC
      for (int i = 0; i < DoFxmax.size(); i++) {
	DoFvalues[DoFxmax[i]] += DeltaX;
      }
      ind++;
      mySolver.solve(DISP);
      for (int k = 0; k < NumMat; k++)
      {
	(PLmaterials[k])->updateStateVariables();
      }
      myModel.writeOutputVTK(outputString, ind);
    }

    for (int s = 0; s < 5; s++) {
      cout << endl << "============================" << endl;
      cout << "Step " << s << endl;
      cout << "Activation Factor: " << ActivationFactor[s] << endl;
      
      if (s == 0)  // Free Contraction = 0, Fixed = 200
      {
	for (int i = 0; i < DoFxmax.size(); i++)
	{
	  DoFid.erase(DoFid.begin() + DoFxmax[i] - i);
	  DoFvalues.erase(DoFvalues.begin() + DoFxmax[i] - i);
	  ViscPotential.setViscosity(0.0E-5);
	  // High  = 1E-4, Low = 1E-5
	}
      }

      for (int k = 0; k < NumMat; k++)
      {
    	(PLmaterials[k])->setTimestep(0.5/1000);
	if(s >= 901 || ActivationFactor[s] < 0.0)
    		(PLmaterials[k])->setActivationMultiplier(0.0);
	else
	    	(PLmaterials[k])->setActivationMultiplier(ActivationFactor[s]);
		// (PLmaterials[k])->setActivationMultiplier(0.0);
      }
      
      ind++;
      mySolver.solve(DISP);
     
      for (int k = 0; k < NumMat; k++)
      {
    	(PLmaterials[k])->updateStateVariables();
      }
      myModel.writeOutputVTK(outputString, ind);
      // myModel.writeField("CubeSmall_", 1);
    }

  // Timing
  time (&end);
  cout << endl << "BVP solved in " << difftime(end,start) << " s" << endl;



  return 0;
}
