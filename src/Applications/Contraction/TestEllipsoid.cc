#include "PassMyoA.h"
#include "CompNeoHookean.h"
#include "FEMesh.h"
#include "EigenEllipticResult.h"
#include "MechanicsModel.h"
#include "EigenNRsolver.h"
#include "LBFGSB.h"

using namespace voom;

int main(int argc, char** argv)
{

  // Timing
  time_t start, end;
  time(&start);



  // Initialize Mesh
  // Assumptions to use this main as is: strip has a face at z=0; tetrahedral mesh
  FEMesh Heart("Mesh/HeartEllipseSmall.node", "Mesh/HeartEllipseSmall.ele");
  FEMesh surfMesh("Mesh/HeartEllipseSmall.node", "Mesh/HeartEllipseSmall.surf");
  string BCfile = "Mesh/HeartEllipseSmall.bc";
 
  cout << endl;
  cout << "Number Of Nodes   : " << Heart.getNumberOfNodes() << endl;
  cout << "Number Of Element : " << Heart.getNumberOfElements() << endl;
  cout << "Mesh Dimension    : " << Heart.getDimension() << endl << endl;

  cout << "Number Of Nodes   : " << surfMesh.getNumberOfNodes() << endl;
  cout << "Number Of Element : " << surfMesh.getNumberOfElements() << endl;
  

    
  // Initialize Material
  uint NumMat =  Heart.getNumberOfElements();
  vector<MechanicsMaterial * > materials;
  materials.reserve(NumMat);
  vector<Vector3d > Fibers; 
  Fibers.reserve(NumMat);

  for (int k = 0; k < NumMat; k++) {
    CompNeoHookean* Mat = new CompNeoHookean(k, 1.0, 1.0);
    materials.push_back(Mat);

    Vector3d N; N << 1.0, 0.0, 0.0;
    Fibers.push_back(N);
  }

  // Load fibers into elements
  Heart.setFibers(Fibers);

  // Initialize Model
  int NodeDoF = 3;
  int PressureFlag = 1;
  int NodalForcesFlag = 0;
  MechanicsModel myModel(&Heart, materials, NodeDoF, PressureFlag, &surfMesh,
			 NodalForcesFlag);
 
  // Initialize Result
  uint myRequest;
  uint PbDoF = (Heart.getNumberOfNodes())*myModel.getDoFperNode();
  EigenEllipticResult myResults(PbDoF, NumMat*2);



  // Run Consistency check
    Real perturbationFactor = 0.1;
    myRequest = 7; // Check both Forces and Stiffness
    Real myH = 1e-6;
    Real myTol = 1e-7;
    // myModel.checkConsistency(myResults, perturbationFactor, myRequest, myH, myTol);
    // myModel.checkDmat(myResults, perturbationFactor, myH, myTol);
  
  // Check on applied pressure
  myModel.updatePressure(1.0);
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
  int NumBC = 0, node = 0, ind = 0;
  vector<int > BCid;
  vector<Real > BCvalues;
  ifstream BCinp(BCfile.c_str());

  BCinp >> NumBC;
  cout << "NumBC = " << NumBC << endl;
  BCid.reserve(NumBC*3);
  BCvalues.reserve(NumBC*3);
  for(int i = 0; i < NumBC; i++) {
    BCinp >> node;
    // cout << "node = " << node << endl;
    for (int j = 0; j < 3; j++) {
      BCid.push_back(node*3 + j);
      BCvalues.push_back(Heart.getX(node, j));
      // cout << BCid[ind] << " " <<  BCvalues[ind] << endl;
      ind++;
    }
  }

  // Print initial configuration
  // myModel.writeOutputVTK("Output/CubeSmall_", 0);
  myModel.writeOutputVTK("Output/Heart_", 0);


  // Loop through pressure
  int NumPsteps = 10;
  Real DeltaP = 0.01;

  // Solver
  Real NRtol = 1.0e-12;
  uint NRmaxIter = 100;
  EigenNRsolver mySolver(&myModel, BCid, BCvalues, CHOL, NRtol, NRmaxIter);  

  for (int s = 1; s <= NumPsteps; s++) {
    myModel.updatePressure(Real(s)*DeltaP);
    mySolver.solve(DISP); 

    myModel.writeOutputVTK("Output/Heart_", s);
  }

 

  // Timing
  time (&end);
  cout << endl << "BVP solved in " << difftime(end,start) << " s" << endl;



  return 0;
}
