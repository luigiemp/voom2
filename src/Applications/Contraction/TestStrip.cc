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
  // Assumptions to use this main as is: strip has a face at x=0; tetrahedral mesh
  // FEMesh Cube("Mesh/Cube6.node", "Mesh/Cube6.ele");
  // FEMesh surfMesh("Mesh/Cube6.node", "Mesh/Cube6Surf.ele");
  FEMesh Cube("Mesh/Strip36.node", "Mesh/Strip36.ele");
  FEMesh surfMesh("Mesh/Strip36.node", "Mesh/Strip36Surf.ele");
  // FEMesh Cube("Mesh/Strip144.node", "Mesh/Strip144.ele");
  // FEMesh surfMesh("Mesh/Strip144.node", "Mesh/Strip144Surf.ele");

 
  cout << endl;
  cout << "Number Of Nodes   : " << Cube.getNumberOfNodes() << endl;
  cout << "Number Of Element : " << Cube.getNumberOfElements() << endl;
  cout << "Mesh Dimension    : " << Cube.getDimension() << endl << endl;
  
    
  // Initialize Material
  uint NumMat =  Cube.getNumberOfElements();
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
  Cube.setFibers(Fibers);

  // Initialize Model
  int NodeDoF = 3;
  int PressureFlag = 1;
  Real Pressure = 1.54;
  int NodalForcesFlag = 1;
  vector<int > ForcesID;
  vector<Real > Forces;
  MechanicsModel myModel(&Cube, materials, NodeDoF, PressureFlag, &surfMesh,
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
  myModel.checkConsistency(myResults, perturbationFactor, myRequest, myH, myTol);
  myModel.checkDmat(myResults, perturbationFactor, myH, myTol);
  
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
  for(int i = 0; i < Cube.getNumberOfNodes(); i++) {
    // All nodes at x = 0
    if ( Cube.getX(i, 0) < 1.0e-12 ) {
      DoFid.push_back(i*3);
      // All nodes at y = 0
      if ( Cube.getX(i, 1) < 1.0e-12 ) {
	DoFid.push_back(i*3 + 1);
      }
      // All nodes at z = 0
      if ( Cube.getX(i, 2) < 1.0e-12 ) {
	DoFid.push_back(i*3 + 2);
      }
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

  // Print initial configuration
  // myModel.writeOutputVTK("Output/CubeSmall_", 0);
  myModel.writeOutputVTK("Output/Strip36_", 0);
    
  // Solver
  Real NRtol = 1.0e-12;
  uint NRmaxIter = 100;
  EigenNRsolver mySolver(&myModel, DoFid, DoFvalues, CHOL, NRtol, NRmaxIter);
  mySolver.solve(DISP); 

  // myModel.writeOutputVTK("Output/CubeSmall_", 1);
  myModel.writeOutputVTK("Output/Strip36_", 1);
    

 
 
   
  // Compute Residual
  myRequest = 2;
  myResults.setRequest(myRequest);
  myModel.compute(myResults);
  VectorXd F = *(myResults._residual);
  Real TotRes = 0.0;
  for (int i = 0; i < DoFid.size(); i ++) {
    cout << "DoF = " <<  DoFid[i] << " F = " << -F(DoFid[i]) << endl;
    TotRes -= F(DoFid[i]);
  }

  // Timing
  time (&end);
  cout << endl << "BVP solved in " << difftime(end,start) << " s" << endl;



  return 0;
}
