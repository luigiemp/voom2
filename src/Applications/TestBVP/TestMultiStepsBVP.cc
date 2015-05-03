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
  FEMesh Cube("../../Mesh/Test/CubeQuad.node", "../../Mesh/Test/CubeQuad.ele");
  FEMesh surfMesh("../../Mesh/Test/CubeQuad.node", "../../Mesh/Test/SurfCubeQuad.ele");
 
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
    // PassMyoA* Mat = new PassMyoA(k, 39.02, 7.31, 100.0, 0.01, 2.85, 2.82);
    PassMyoA* Mat = new PassMyoA(k, 35.19, 7.06, 100.0, 0.025, 2.87, 2.82);
    Vector3d N; N << 1.0/sqrt(3.0), 1.0/sqrt(3.0), 1.0/sqrt(3.0);
    materials.push_back(Mat);
    // materials.push_back(new CompNeoHookean(k, 10.0, 3.0) );

    Fibers.push_back(N);
  }

  // Load fibers into elements
  Cube.setFibers(Fibers);

  // Initialize Model
  int NodeDoF = 3;
  int PressureFlag = 1;
  Real Pressure = 1.0;
  int NodalForcesFlag = 1;
  vector<int > ForcesID;
  vector<Real > Forces;
  MechanicsModel myModel(&Cube, materials, NodeDoF, PressureFlag, Pressure, &surfMesh,
			 NodalForcesFlag, &ForcesID, &Forces);
 
  // Initialize Result
  uint PbDoF = (Cube.getNumberOfNodes())*myModel.getDoFperNode();
  EigenEllipticResult myResults(PbDoF, NumMat*2);

  // Run Consistency check
  Real perturbationFactor = 0.1;
  uint myRequest = 7; // Check both Forces and Stiffness
  Real myH = 1e-6;
  Real myTol = 1e-7;
  myModel.checkConsistency(myResults, perturbationFactor, myRequest, myH, myTol);
  myModel.checkDmat(myResults, perturbationFactor, myH, myTol);
  
    // Solve for displacement first
    // Print initial configuration
    myModel.writeOutputVTK("CubeQuad_", 0);

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
    
    // Quad tet cube applied pressure //
    vector<int > DoFid(13,0);
    vector<Real > DoFvalues(13, 0.0);
    DoFid[0] = 0;  DoFid[1] = 1;  DoFid[2] = 2;  
    DoFid[3] = 5;  DoFid[4] = 8;  DoFid[5] = 11;  
    DoFid[6] = 26; DoFid[7] = 29; DoFid[8] = 35;          
    DoFid[9] = 38; DoFid[10] = 47;
    DoFid[11] = 3; DoFid[12] = 10;
 
    for(int i = 0; i < DoFid.size(); i++) {
      DoFvalues[i] = Cube.getX(floor(double(DoFid[i])/3.0) ,DoFid[i]%3);
    }

   
    // Loop through pressure
    int NumPsteps = 10;
    Real DeltaP = 0.2;
    // Solver
    Real NRtol = 1.0e-12;
    uint NRmaxIter = 100;
    EigenNRsolver mySolver(&myModel, DoFid, DoFvalues, CHOL, NRtol, NRmaxIter);
    
    for (int i = 1; i <= NumPsteps; i++) {
      myModel.updatePressure(Real(i)*DeltaP);
      mySolver.solve(DISP); 

      myModel.writeOutputVTK("CubeQuad_", i);
      myModel.writeField("CubeQuad_", i);
    }


    // int m = 5; 
    // double factr = 1.0e+1;
    // double pgtol = 1.0e-5;
    // int iprint = 0;
    // int maxIterations = -1;
    
    // vector<int > nbd(PbDoF, 0);
    // vector<double > lowerBound(PbDoF, 0.0);
    // vector<double > upperBound(PbDoF, 0.0);
    // for(int i = 0; i < DoFid.size(); i++) {
    //   nbd[DoFid[i]] = 2;
    //   lowerBound[DoFid[i]] = DoFvalues[i];
    //   upperBound[DoFid[i]] = DoFvalues[i];
    // }
    // cout << "PbDoF = " << PbDoF << endl;
    // for (int j = 0; j < PbDoF; j++) {
    //   cout << nbd[j] << " " << lowerBound[j] << " " << upperBound[j] << endl;
    // }
    
    // LBFGSB mySolver(& myModel, & myResults, m, factr, pgtol, iprint, maxIterations);
    // mySolver.setBounds(nbd, lowerBound, upperBound);
    // mySolver.solve();
    
    // // Print final configuration
    // myModel.writeOutputVTK("CubeQuad_", 1);
    
    // // Print _field
    // cout << endl;
    // myModel.printField();
    // cout << endl;
    
    // Compute Residual
    myRequest = 2;
    myResults.setRequest(myRequest);
    myModel.compute(myResults);
    VectorXd F = *(myResults._residual);
    for (int i = 0; i < DoFid.size(); i ++) {
      ForcesID.push_back( DoFid[i] );
      Forces.push_back( -F(DoFid[i]) );
    }
 // End of displacement solution

  // // Next recompute material properties using EMFO
  // {
  //   ifstream BCinp(DispFile.c_str());
  //   int DoF;

  //   BCinp >> DoF;
  //   vector<int > DoFid(DoF, 0);
  //   vector<Real > DoFvalues(DoF, 0);
  //   for(int i = 0; i < DoF; i++) {
  //     DoFid[i] = i;
  //     BCinp >> DoFvalues[i];
  //     // cout <<  DoFid[i] << " " <<  DoFvalues[i] << endl;
  //   }
   
  //   // Change Material properties
  //   // Find unique material parameters
  //   set<MechanicsMaterial *> UNIQUEmaterials;
  //   for (uint i = 0; i < materials.size(); i++) 
  //     UNIQUEmaterials.insert(materials[i]);

  //   srand (time(NULL));
  //   for (set<MechanicsMaterial *>::iterator MatIt = UNIQUEmaterials.begin(); MatIt != UNIQUEmaterials.end(); MatIt++) {
  //     int MatID = (*MatIt)->getMatID();
  //     vector<Real > MatProp = (*MatIt)->getMaterialParameters();
  //     int NumPropPerMat = MatProp.size();
  //     for (int m = 0; m < NumPropPerMat; m++) {
  //   	MatProp[m] *= double(rand())/double(RAND_MAX);
  // 	cout << MatProp[m] << " ";
  //     }
  //     cout << endl;
  //     (*MatIt)->setMaterialParameters(MatProp);
  //   }

  //   myModel.writeOutputVTK("CubeQuad_", 2);

  //   // Solve for correct material properties
  //   Real NRtol = 1.0e-12;
  //   uint NRmaxIter = 10;
  //   EigenNRsolver mySolver(&myModel, DoFid, DoFvalues, CHOL, NRtol, NRmaxIter);
  //   mySolver.solve(MAT); 
 
  //   // Print found material properties
  //   for (set<MechanicsMaterial *>::iterator MatIt = UNIQUEmaterials.begin(); MatIt != UNIQUEmaterials.end(); MatIt++) {
  //     int MatID = (*MatIt)->getMatID();
  //     vector<Real > MatProp = (*MatIt)->getMaterialParameters();
  //     int NumPropPerMat = MatProp.size();
  //     for (int m = 0; m < NumPropPerMat; m++) {
  // 	cout << MatProp[m] << " ";
  //     }
  //     cout << endl;
  //   }

  //   myModel.writeOutputVTK("CubeQuad_", 3);

  // } // End of material properties solution


  // Timing
  time (&end);
  cout << endl << "BVP solved in " << difftime(end,start) << " s" << endl;



  return 0;
}
