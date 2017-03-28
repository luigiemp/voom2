#include "PassMyoA.h"
#include "CompNeoHookean.h"
#include "Humphrey_Compressible.h"
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

  // Initialize Mesh
  // FEMesh Ellipsoid("../Input/Swine2/Swine2.node", "../Input/Swine2/Swine2.ele");
  // FEMesh surfMesh("../Input/Swine2/Swine2.node", "../Input/Swine2/Swine2.endoEle");
  // string FiberFile = "../Input/Swine2/Swine2.fiber";
  // string BCfile = "../Input/Swine2/Swine2.bc";
  // string outputString = "Output/Swine2/Swine_";

  FEMesh Ellipsoid("../Input/Swine2/Swine2.node", "../Input/Swine2/Swine2.ele");
  FEMesh surfMesh("../Input/Swine2/Swine2.node", "../Input/Swine2/Swine2.endoEle");
  string FiberFile = "../Input/Swine2/Swine2.fiber";
  string BCfile = "../Input/Swine2/Swine2.bc";
  string outputString = "Output/Swine2/SwineStiff_";

  // FEMesh Ellipsoid("../Input/Swine2/Swine2expanded.node", "../Input/Swine2/Swine2expanded.ele");
  // FEMesh surfMesh("../Input/Swine2/Swine2expanded.node", "../Input/Swine2/Swine2expanded.endoEle");
  // string FiberFile = "../Input/Swine2/Swine2expanded.fiber";
  // string BCfile = "../Input/Swine2/Swine2expanded.bc";
  // string outputString = "Output/Swine2/SwineExp_";
 
  ifstream FiberInp(FiberFile.c_str());

  cout << endl;
  cout << "Number Of Nodes   : " << Ellipsoid.getNumberOfNodes() << endl;
  cout << "Number Of Element : " << Ellipsoid.getNumberOfElements() << endl;
  cout << "Mesh Dimension    : " << Ellipsoid.getDimension() << endl << endl;
  
  // Initialize Material
  int NumQP  = 4;
  int NumMat =  Ellipsoid.getNumberOfElements()*NumQP;
  
  vector<MechanicsMaterial * > materials;
  materials.reserve(NumMat);

  for (int k = 0; k < NumMat; k++) {
    // Load Fiber Data:
    vector <Vector3d> dirvec(3, Vector3d::Zero(3,1));
    FiberInp >> dirvec[0][0]; FiberInp >> dirvec[0][1]; FiberInp >> dirvec[0][2];
    FiberInp >> dirvec[1][0]; FiberInp >> dirvec[1][1]; FiberInp >> dirvec[1][2];
    FiberInp >> dirvec[2][0]; FiberInp >> dirvec[2][1]; FiberInp >> dirvec[2][2];

    // PLmaterials.push_back(&PassiveMat);
    // CompNeoHookean* PassMat = new CompNeoHookean(k, 2.0, 2.0);
    // Humphrey_Compressible* PassMat = new Humphrey_Compressible(k, 15.98, 55.85, 0.0, -33.27, 30.21, 3.590, 64.62, dirvec);
    Humphrey_Compressible* PassMat = new Humphrey_Compressible(k, 2.0*15.98, 2.0*55.85, 0.0, -2.0*33.27, 2.0*30.21, 2.0*3.590, 2.0*64.62, dirvec);
    // PlMat->setDirectionVectors(dirvec);
    
    materials.push_back(PassMat);
  }

  cout << "here" << endl;

  // Initialize Model
  int NodeDoF = 3;
  int PressureFlag = 1;
  Real Pressure = -0.75;
  int NodalForcesFlag = 0;
  vector<int > ForcesID;
  vector<Real > Forces;
  MechanicsModel myModel(&Ellipsoid, materials, NodeDoF, PressureFlag, &surfMesh, NodalForcesFlag);
  myModel.updatePressure(Pressure);
  myModel.updateNodalForces(&ForcesID, &Forces);
 
  // Initialize Result
  uint myRequest;
  uint PbDoF = (Ellipsoid.getNumberOfNodes())*myModel.getDoFperNode();
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
  // myModel.writeOutputVTK("Output/Swine2/SwineExp_", 0);

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
  BCnodes.reserve(NumBC*3);
  BCid.reserve(NumBC*3);
  BCvalues.reserve(NumBC*3);
  for(int i = 0; i < NumBC; i++) {
    BCinp >> node;
    BCnodes.push_back(node);
    for (int j = 0; j < 3; j++) {
      BCid.push_back(node*3 + j);
      BCvalues.push_back(Ellipsoid.getX(node, j));
      // cout << BCid[ind] << " " <<  BCvalues[ind] << endl;
      ind++;
    }
  }
    


  // Solver
  Real NRtol = 1.0e-8;
  uint NRmaxIter = 100;
  EigenNRsolver mySolver(&myModel, BCid, BCvalues, CHOL, NRtol, NRmaxIter);

  for (int s = 0; s < 21; s++) {
    cout << "Step " << s << endl;
    myModel.updatePressure(Pressure*s);

    mySolver.solve(DISP); 
    
    myModel.writeOutputVTK(outputString, s);
  }

  // Timing
  time (&end);
  cout << endl << "BVP solved in " << difftime(end,start) << " s" << endl;

  return 0;
}
