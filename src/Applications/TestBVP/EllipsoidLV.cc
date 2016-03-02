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
//#include "EigenResult.h"
#include "MechanicsModel.h"
#include "EigenNRsolver.h"
// #include "LBFGSB.h"

using namespace voom;

int main(int argc, char** argv)
{

  // Timing
  time_t start, end;
  time(&start);

  // Conduction Velocity
  double cv = 0.05;

  // OutputString
  string outputString = "LV_Ellipsoid";



  // Initialize Mesh
  // Assumptions to use this main as is: strip has a face at z=0; tetrahedral mesh
  FEMesh Cube("Small_A.node", "Small_A.ele");
  FEMesh surfMesh("Small_A.node", "Small_A.surfEle");
  string FiberFile = "Small_A.fiber";
  string BCfile = "Small_A.bc";
  ifstream FiberInp(FiberFile.c_str());


  Real xmax = 1.0;
  // Real xmax = 6.0; 

  cout << endl;
  cout << "Number Of Nodes   : " << Cube.getNumberOfNodes() << endl;
  cout << "Number Of Element : " << Cube.getNumberOfElements() << endl;
  cout << "Mesh Dimension    : " << Cube.getDimension() << endl << endl;
  
  // Initialize Material
  uint NumMat =  Cube.getNumberOfElements();
  
  vector<MechanicsMaterial * > PLmaterials;
  PLmaterials.reserve(NumMat);

  CompNeoHookean PassiveMat(0, 0.4, 0.04);
  CompNeoHookean ActiveMat(0, 0.4, 0.04);
  APForceVelPotential TestPotential(1.0, 500.0);
  // HillForceVelPotential TestPotential(4.4*pow(10,-3), .01*0.59, 25);
  // BlankPotential TestPotential;

  BlankViscousPotential ViscPotential;
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
    // Load Fiber Data:
    vector <Vector3d> dirvec(3, Vector3d::Zero(3,1));
    FiberInp >> dirvec[0][0]; FiberInp >> dirvec[0][1]; FiberInp >> dirvec[0][2];
    FiberInp >> dirvec[1][0]; FiberInp >> dirvec[1][1]; FiberInp >> dirvec[1][2];
    FiberInp >> dirvec[2][0]; FiberInp >> dirvec[2][1]; FiberInp >> dirvec[2][2];

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
  EigenResult myResults(PbDoF, NumMat*2);

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
    myModel.writeOutputVTK("LV_Ellipsoid_", 0);

    // Check on applied pressure
    myRequest = 2;
    myResults.setRequest(myRequest);
    myModel.compute(&myResults);
    VectorXd Fstart = *(myResults._residual);
    Real pX = 0.0, pY = 0.0, pZ = 0.0;
    for (int i = 0; i < Fstart.size(); i += 3) {
      pX += Fstart(i);
      pY += Fstart(i+1);
      pZ += Fstart(i+2);
    }
    cout << endl << "Pressure = " << pX << " " << pY << " " << pZ << endl << endl;
     
    // EBC
    cout << "********" << " Setting up EBCs " << "********" << endl;
    int NumBC = 0, node = 0, ind = 0;;
    vector<int > BCnodes;
    vector<int > BCid;
    vector<Real > BCvalues;
    ifstream BCinp(BCfile.c_str());
   
    if (BCinp.is_open())
	cout << "BC File Opened Successfully!" << endl; 
    else
	cout << "ERROR: BC File Failed to Open!" << endl;
   
    BCinp >> NumBC;
    cout << "Number of Nodes with EBC: " << NumBC << endl;
    BCid.reserve(NumBC*3);
    BCvalues.reserve(NumBC*3);
    for(int i = 0; i < NumBC; i++) {
      BCinp >> node;
      BCnodes.push_back(node);
      for (int j = 0; j < 3; j++) {
        BCid.push_back(node*3 + j);
        BCvalues.push_back(Cube.getX(node, j));
        cout << BCid[ind] << " " <<  BCvalues[ind] << endl;
        ind++;
      }
    }
    

    // Lin tet cube - applied displacements //
    // for (uint i = 18; i < 27; i++) { 
    //   DoFid[i-5] = i*3 + 2; // Lower top in z direction by 0.1
    //   DoFvalues[i-5] = 1.95;
    // }
    
    // Solver
    Real NRtol = 1.0e-12;
    uint NRmaxIter = 100;
    EigenNRsolver mySolver(&myModel, BCid, BCvalues, CHOL, NRtol, NRmaxIter);

    ind = 0;
    /*
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
    */
    for (int s = 0; s < 900; s++) {
      cout << "Step " << s << endl;
      // cout << "Activation Factor: " << ActivationFactor[s] << endl;
      
      if (s == 0)  // Free Contraction = 0, Fixed = 200
      {
	  ViscPotential.setViscosity(0.0E-5);
      }

      for (int k = 0; k < NumMat; k++)
      {
    	(PLmaterials[k])->setTimestep(0.5/1000);
	if(s >= 795 || ActivationFactor[s] < 0.0)
    		(PLmaterials[k])->setActivationMultiplier(0.0);
	else
	    	(PLmaterials[k])->setActivationMultiplier(ActivationFactor[s] * 0.1);
      }
      
      ind++;
      mySolver.solve(DISP);
     
      for (int k = 0; k < NumMat; k++)
      {
    	(PLmaterials[k])->updateStateVariables();
      }
      myModel.writeOutputVTK(outputString, ind);
      // myModel.writeField("CubeSmall_", 1);
      
      // Check EBCs after Solve:
      // Get the Field from model
      vector <double > tempField(Cube.getNumberOfNodes() * 3, 0.0);
      myModel.getField(tempField);
      for (int node_iter = 0; node_iter < NumBC; node_iter++)
      {
        for (int dim_iter = 0; dim_iter < 3; dim_iter++)
          cout << BCnodes[node_iter] * 3 + dim_iter << "\t" << tempField[BCnodes[node_iter] * 3 + dim_iter] - Cube.getX(BCnodes[node_iter], dim_iter) << endl;
      }
    }

  // Timing
  time (&end);
  cout << endl << "BVP solved in " << difftime(end,start) << " s" << endl;

  return 0;
}
