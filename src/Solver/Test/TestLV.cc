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
  FEMesh LVmesh("../../Mesh/Test/CoarseLV.node", "../../Mesh/Test/CoarseLV.ele");
  FEMesh LVsurf("../../Mesh/Test/CoarseLV.node", "../../Mesh/Test/CoarseLV.surf");
  string BCfile = "CoarseLV.BaseSurfBC";
 
  cout << endl;
  cout << "Number Of Nodes   : " << LVmesh.getNumberOfNodes() << endl;
  cout << "Number Of Element : " << LVmesh.getNumberOfElements() << endl;
  cout << "Mesh Dimension    : " << LVmesh.getDimension() << endl << endl;
  
  // Initialize Material
  uint NumMat =  LVmesh.getNumberOfElements();
  vector<MechanicsMaterial * > materials;
  materials.reserve(NumMat);
  for (int k = 0; k < NumMat; k++) {
    // // PassMyoA* Mat = new PassMyoA(30.48, 7.25, 0.1, 1.60, 2.82);
    // PassMyoA* Mat = new PassMyoA(10.48, 7.25, 100.0, 1.0, 2.60, 2.82);
    // // Vector3d N; N << sqrt(2.0)/2.0, sqrt(2.0)/2.0, 0.0;
    // Vector3d N; N << 0.0, 0.0, 1.0;
    // Mat->setN(N);
    // // vector<Real > temp  = Mat->getMaterialParameters();
    // // cout << temp[0] << " " << temp[1] << " " <<  temp[2] << endl;
    // materials.push_back(Mat);
    materials.push_back(new CompNeoHookean(k, 10.0, 10.0) );
  }

  // Initialize Model
  int NodeDoF = 3;
  int PressureFlag = 1;
  Real Pressure = 10.0;
  MechanicsModel myModel(&LVmesh, materials, NodeDoF, PressureFlag, Pressure, &LVsurf);
 
  // Initialize Result
  uint PbDoF = (LVmesh.getNumberOfNodes())*myModel.getDoFperNode();
  EigenEllipticResult myResults(PbDoF, NumMat);
 
  // Print initial configuration
  myModel.writeOutputVTK("LVpassiveFilling_", 0);

  // Check on applied pressure
  int myRequest = 2;
  myResults.setRequest(myRequest);
  myModel.compute(myResults);
  VectorXd F = *(myResults._residual);
  Real pX = 0.0, pY = 0.0, pZ = 0.0;
  for (int i = 0; i < F.size(); i += 3) {
    pX += F(i);
    pY += F(i+1);
    pZ += F(i+2);
  }
  cout << endl << "Pressure = " << pX << " " << pY << " " << pZ << endl << endl;


  // EBC
  int NumBC = 0, node = 0, ind = 0;;
  vector<int > DoFid;
  vector<Real > DoFvalues;
  ifstream BCinp(BCfile.c_str());

  BCinp >> NumBC;
  DoFid.reserve(NumBC*3);
  DoFvalues.reserve(NumBC*3);
  for(int i = 0; i < NumBC; i++) {
    BCinp >> node;
    for (int j = 0; j < 3; j++) {
      DoFid.push_back(node*3 + j);
      DoFvalues.push_back(LVmesh.getX(node, j));
      // cout << DoFid[ind] << " " <<  DoFvalues[ind] << endl;
      ind++;
    }
  }
  
  // Solver
  Real NRtol = 1.0e-8;
  uint NRmaxIter = 100;
  EigenNRsolver mySolver(&myModel, DoFid, DoFvalues, CHOL, NRtol, NRmaxIter);
  mySolver.solve(DISP); 

  // int m = 5; 
  // double factr = 1.0e+1;
  // double pgtol = 1.0e-5;
  // int iprint = 50;
  // int maxIterations = 10000;
  
  // vector<int > nbd(PbDoF, 0);
  // vector<double > lowerBound(PbDoF, 0.0);
  // vector<double > upperBound(PbDoF, 0.0);
  // for(int i = 0; i < DoFid.size(); i++) {
  //   nbd[DoFid[i]] = 2;
  //   lowerBound[DoFid[i]] = DoFvalues[i];
  //   upperBound[DoFid[i]] = DoFvalues[i];
  // }
  // // cout << "PbDoF = " << PbDoF << endl;
  // // for (int j = 0; j < PbDoF; j++) {
  // //   cout << nbd[j] << " " << lowerBound[j] << " " << upperBound[j] << endl;
  // // }
  
  // LBFGSB mySolver(& myModel, & myResults, m, factr, pgtol, iprint, maxIterations);
  // mySolver.setBounds(nbd, lowerBound, upperBound);
  // mySolver.solve();
  
  // Print final configuration
  myModel.writeOutputVTK("LVpassiveFilling_", 1);



  // Timing
  time (&end);
  cout << endl << "BVP solved in " << difftime(end,start) << " s" << endl;



  return 0;
}
